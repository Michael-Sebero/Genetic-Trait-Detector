#!/usr/bin/env python3
"""Offline trait and variant report for consumer raw genotype data.

Reads 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA, Illumina final
reports, generic rsid/chromosome/position/genotype tables and single-sample VCF,
plain or gzip/bzip2/zip compressed. Every database allele is stored on the GRCh37
plus strand. Non-palindromic calls reported on the opposite strand are complemented;
anything that fits neither orientation is reported as an allele conflict.

The text report lists every database SNP found in the file, split into traits present
and traits not present; every result title gets a colour no other title uses.

  genetic-trait-detector.py genome.zip
  genetic-trait-detector.py genome.zip --snp rs1006737
  genetic-trait-detector.py genome.txt --format json > report.json
  genetic-trait-detector.py genome.txt --pgs PGS000001.txt.gz
  genetic-trait-detector.py --list
  genetic-trait-detector.py --selftest
"""
from __future__ import annotations

import argparse
import bz2
import colorsys
import gzip
import io
import json
import math
import os
import re
import shutil
import sys
import tempfile
import textwrap
import zipfile
from dataclasses import dataclass, field

CATEGORIES = {
    "appearance": ("Appearance & pigmentation", "33"),
    "senses": ("Senses & physiology", "36"),
    "diet": ("Diet & metabolism", "32"),
    "substances": ("Alcohol & nicotine", "35"),
    "drugs": ("Drug response", "34"),
    "health": ("Health", "31"),
    "neuro": ("Neuro & psychiatric", "95"),
    "longevity": ("Longevity", "94"),
}

EVIDENCE = {
    "A": "clinical grade: guideline-backed or established pathogenic/protective",
    "B": "strong association, large effect",
    "C": "replicated association, small effect",
    "D": "functional variant, trait associations inconsistent",
}
EVIDENCE_SHORT = "A clinical grade | B strong, large effect | C replicated, small effect | D weak"

COMPLEMENT = str.maketrans("ACGT", "TGCA")
BASES = frozenset("ACGTDI")
PAR_X = (2781479, 154931044)


@dataclass(frozen=True)
class Variant:
    rsid: str
    gene: str
    cat: str
    ref: str
    alt: str
    ev: str
    text: tuple
    src: str
    chrom: str = ""
    pos: int = 0
    hidden: bool = False
    levels: str = "000"
    rare: bool = False
    aliases: tuple = ()
    trait: str = ""
    has: str = "011"

    @property
    def palindromic(self) -> bool:
        return {self.ref, self.alt} in ({"A", "T"}, {"C", "G"})


def _v(rsid, gene, cat, ref, alt, ev, text, src, loc="", hidden=False, levels="000", rare=False, aliases=(),
       trait="", has="011"):
    chrom, _, pos = loc.partition(":")
    return Variant(rsid, gene, cat, ref, alt, ev, tuple(text), src, chrom, int(pos or 0),
                   hidden, levels, rare, tuple(aliases), trait, has)


MC1R_SRC = "Valverde 1995 Nat Genet 11:328; Sulem 2007 Nat Genet 39:1443"
IRIS_SRC = "Walsh 2011 Forensic Sci Int Genet 5:170"

DB = [
    _v("rs12913832", "HERC2", "appearance", "A", "G", "B",
       ("brown-eye genotype", "brown or intermediate (green/hazel) eyes most likely", "blue/light eyes expected"),
       "Sturm 2008 Am J Hum Genet 82:424; Eiberg 2008 Hum Genet 123:177", "15:28365618",
       trait="Blue/light eyes (HERC2)", has="001"),
    _v("rs1800407", "OCA2", "appearance", "C", "T", "C",
       ("no R419Q", "R419Q: shifts prediction from brown toward blue/intermediate", "R419Q homozygous"),
       IRIS_SRC, "15:28230318", hidden=True, trait="OCA2 R419Q", has="011"),
    _v("rs12896399", "SLC24A4", "appearance", "G", "T", "C",
       ("no light-pigmentation allele", "one allele associated with light eyes and hair",
        "two alleles associated with light eyes and hair"),
       "Sulem 2007 Nat Genet 39:1443", "14:92773663", hidden=True,
       trait="Light eye/hair allele (SLC24A4)", has="011"),
    _v("rs16891982", "SLC45A2", "appearance", "C", "G", "B",
       ("374Leu/Leu: ancestral darker-pigmentation genotype, common outside Europe",
        "374Leu/Phe: one light-pigmentation allele",
        "374Phe/Phe: light-pigmentation genotype typical of Europeans"),
       "Graf 2005 Hum Mutat 25:278; " + IRIS_SRC, "5:33951693",
       trait="Light-pigmentation allele (SLC45A2 374Phe)", has="011"),
    _v("rs1393350", "TYR", "appearance", "G", "A", "C",
       ("no TYR light-eye allele", "blue vs green eyes OR 1.52 per A; more sun sensitivity",
        "blue vs green eyes OR 1.52 per A (two copies); more sun sensitivity"),
       "Sulem 2007 Nat Genet 39:1443", "11:89011046", hidden=True, trait="Light-eye allele (TYR)", has="011"),
    _v("rs12203592", "IRF4", "appearance", "C", "T", "C",
       ("no IRF4 effect allele", "more freckling and sun sensitivity",
        "more freckling and sun sensitivity (two alleles)"),
       "Han 2008 PLoS Genet 4:e1000074; Praetorius 2013 Cell 155:1022", "6:396321",
       trait="Freckling and sun sensitivity (IRF4)", has="011"),
    _v("rs1805006", "MC1R", "appearance", "C", "A", "B", ("no D84E", "D84E (R allele)", "D84E homozygous"),
       MC1R_SRC, "16:89985918", hidden=True, trait="MC1R D84E", has="011"),
    _v("rs11547464", "MC1R", "appearance", "G", "A", "B", ("no R142H", "R142H (R allele)", "R142H homozygous"),
       MC1R_SRC, "16:89986091", hidden=True, trait="MC1R R142H", has="011"),
    _v("rs1805007", "MC1R", "appearance", "C", "T", "B", ("no R151C", "R151C (R allele)", "R151C homozygous"),
       MC1R_SRC, "16:89986117", hidden=True, trait="MC1R R151C", has="011"),
    _v("rs1805008", "MC1R", "appearance", "C", "T", "B", ("no R160W", "R160W (R allele)", "R160W homozygous"),
       MC1R_SRC, "16:89986144", hidden=True, trait="MC1R R160W", has="011"),
    _v("rs1805009", "MC1R", "appearance", "G", "C", "B", ("no D294H", "D294H (R allele)", "D294H homozygous"),
       MC1R_SRC, "16:89986546", hidden=True, aliases=("i3002507",), trait="MC1R D294H", has="011"),
    _v("rs1805005", "MC1R", "appearance", "G", "T", "C", ("no V60L", "V60L (r allele)", "V60L homozygous"),
       MC1R_SRC, "16:89985844", hidden=True, trait="MC1R V60L", has="011"),
    _v("rs2228479", "MC1R", "appearance", "G", "A", "C", ("no V92M", "V92M (r allele)", "V92M homozygous"),
       MC1R_SRC, "16:89985940", hidden=True, trait="MC1R V92M", has="011"),
    _v("rs885479", "MC1R", "appearance", "G", "A", "C", ("no R163Q", "R163Q (r allele)", "R163Q homozygous"),
       MC1R_SRC, "16:89986154", hidden=True, trait="MC1R R163Q", has="011"),
    _v("rs12821256", "KITLG", "appearance", "T", "C", "C",
       ("no KITLG blond allele", "one blond-hair-associated allele", "two blond-hair-associated alleles"),
       "Sulem 2007 Nat Genet 39:1443; Guenther 2014 Nat Genet 46:748", "12:89328335",
       trait="Blond-hair allele (KITLG)", has="011"),
    _v("rs35264875", "TPCN2", "appearance", "A", "T", "C",
       ("no M484L", "M484L: blond-hair-associated allele", "M484L homozygous: two blond-hair-associated alleles"),
       "Sulem 2008 Nat Genet 40:835", "11:68846399", trait="Blond-hair allele (TPCN2 M484L)", has="011"),
    _v("rs3829241", "TPCN2", "appearance", "G", "A", "C",
       ("no G734E", "G734E: blond-hair-associated allele", "G734E homozygous: two blond-hair-associated alleles"),
       "Sulem 2008 Nat Genet 40:835", "11:68855363", trait="Blond-hair allele (TPCN2 G734E)", has="011"),
    _v("rs11803731", "TCHH", "appearance", "T", "A", "C",
       ("no straight-hair allele (curlier hair more likely)", "one straight-hair allele",
        "two straight-hair alleles; most common European genotype (T/T in gene orientation)"),
       "Medland 2009 Am J Hum Genet 85:750", "1:152083325", trait="Straight-hair allele (TCHH)", has="011"),
    _v("rs3827760", "EDAR", "appearance", "A", "G", "B",
       ("370Val/Val: ancestral EDAR, typical in Europeans and Africans",
        "370Val/Ala: one hair-thickness allele",
        "370Ala/Ala: thick straight scalp hair; common in East Asians and Native Americans"),
       "Fujimoto 2008 Hum Mol Genet 17:835; Kamberov 2013 Cell 152:691", "2:109513601",
       trait="Thick straight hair allele (EDAR 370Ala)", has="011"),
    _v("rs6152", "AR", "appearance", "G", "A", "C",
       ("common allele; higher odds of male-pattern baldness (G allele OR 2.68)",
        "one lower-risk A allele",
        "lower-risk A allele; A carriers were 24% of non-bald vs 7-13% of bald men"),
       "Hillmer 2005 Am J Hum Genet 77:140; Hayes 2005 Cancer Epidemiol Biomarkers Prev 14:993; "
       "Zhuo 2012 PMID 21981665",
       "X:0", trait="Male-pattern baldness risk allele (AR)", has="110"),

    _v("rs17822931", "ABCC11", "senses", "C", "T", "B",
       ("wet earwax", "wet earwax (dry-type carrier)", "dry earwax and reduced axillary odour"),
       "Yoshiura 2006 Nat Genet 38:324; Martin 2010 J Invest Dermatol 130:529", "16:48258198",
       trait="Dry earwax (ABCC11)", has="001"),
    _v("rs72921001", "OR6A2", "senses", "A", "C", "C",
       ("lowest odds of cilantro tasting soapy (A allele OR 0.81)", "intermediate odds of cilantro tasting soapy",
        "highest odds of cilantro tasting soapy; locus explains ~0.5% of variance"),
       "Eriksson 2012 Flavour 1:22", "11:6889648", trait="Soapy-cilantro allele (OR6A2)", has="011"),
    _v("rs4481887", "OR2M7", "senses", "G", "A", "C",
       ("highest odds of asparagus-metabolite anosmia",
        "likely able to smell asparagus metabolites (dominant OR 0.48 for anosmia)",
        "most likely able to smell asparagus metabolites"),
       "Eriksson 2010 PLoS Genet 6:e1000993", "1:248496863",
       trait="Can smell asparagus metabolites (OR2M7)", has="011"),
    _v("rs10427255", "ZEB2", "senses", "T", "C", "C",
       ("baseline odds of photic sneeze reflex", "higher odds of photic sneeze reflex (OR 1.32 per C)",
        "highest odds of photic sneeze reflex"),
       "Eriksson 2010 PLoS Genet 6:e1000993", "2:146125523", trait="Photic sneeze reflex allele (ZEB2)", has="011"),
    _v("rs1815739", "ACTN3", "senses", "C", "T", "B",
       ("RR: alpha-actinin-3 present in fast fibres; enriched in sprint/power athletes", "RX: alpha-actinin-3 present",
        "XX: alpha-actinin-3 deficient (~18% of Europeans); under-represented in elite sprinters"),
       "Yang 2003 Am J Hum Genet 73:627", "11:66328095", trait="Alpha-actinin-3 deficiency (ACTN3 XX)", has="001"),
    _v("rs713598", "TAS2R38", "senses", "C", "G", "B", ("Ala49 (AVI)", "Ala49/Pro49", "Pro49 (PAV)"),
       "Kim 2003 Science 299:1221", "7:141673345", hidden=True, trait="TAS2R38 Pro49", has="011"),
    _v("rs1726866", "TAS2R38", "senses", "A", "G", "B", ("Val262 (AVI)", "Val262/Ala262", "Ala262 (PAV)"),
       "Kim 2003 Science 299:1221", "7:141672705", hidden=True, trait="TAS2R38 Ala262", has="011"),
    _v("rs10246939", "TAS2R38", "senses", "T", "C", "B", ("Ile296 (AVI)", "Ile296/Val296", "Val296 (PAV)"),
       "Kim 2003 Science 299:1221", "7:141672604", hidden=True, trait="TAS2R38 Val296", has="011"),

    _v("rs4988235", "MCM6/LCT", "diet", "G", "A", "B",
       ("lactase non-persistence likely (adult lactose malabsorption)", "lactase persistence (dominant)",
        "lactase persistence"),
       "Enattah 2002 Nat Genet 30:233; Tishkoff 2007 Nat Genet 39:31 (non-European persistence alleles not assessed)",
       "2:136608646", levels="100", trait="Lactase persistence (MCM6/LCT)", has="011"),
    _v("rs762551", "CYP1A2", "diet", "C", "A", "C",
       ("slower caffeine metabolism (*1F absent)", "slower caffeine metabolism",
        "rapid caffeine metabolism (*1F/*1F, highly inducible)"),
       "Sachse 1999 Br J Clin Pharmacol 47:445; Cornelis 2006 JAMA 295:1135", "15:75041917",
       trait="Fast caffeine metabolism (CYP1A2 *1F/*1F)", has="001"),
    _v("rs601338", "FUT2", "diet", "G", "A", "B",
       ("secretor", "secretor (non-secretor carrier)",
        "non-secretor (W154X/W154X): largely resistant to symptomatic infection by common norovirus strains"),
       "Lindesmith 2003 Nat Med 9:548; Thorven 2005 J Virol 79:15351", "19:49206674",
       trait="Non-secretor (FUT2)", has="001"),
    _v("rs9939609", "FTO", "diet", "T", "A", "C",
       ("no FTO risk allele", "one obesity-risk allele", "~3 kg heavier on average; obesity OR 1.67 vs T/T"),
       "Frayling 2007 Science 316:889", "16:53820527", trait="Obesity-risk allele (FTO)", has="011"),
    _v("rs17782313", "MC4R", "diet", "T", "C", "C",
       ("no MC4R risk allele", "one BMI-raising allele", "two BMI-raising alleles"),
       "Loos 2008 Nat Genet 40:768", "18:57851097", trait="BMI-raising allele (MC4R)", has="011"),
    _v("rs7903146", "TCF7L2", "diet", "C", "T", "C",
       ("no TCF7L2 risk allele", "type 2 diabetes relative risk 1.45", "type 2 diabetes relative risk 2.41"),
       "Grant 2006 Nat Genet 38:320", "10:114758349", levels="012",
       trait="Type 2 diabetes risk allele (TCF7L2)", has="011"),
    _v("rs2282679", "GC", "diet", "T", "G", "C",
       ("no GC low-vitamin-D allele", "modestly lower serum 25(OH)D",
        "lower serum 25(OH)D, higher odds of insufficiency"),
       "Wang 2010 Lancet 376:180", "4:72608383", trait="Lower vitamin D allele (GC)", has="011"),
    _v("rs662799", "APOA5", "diet", "A", "G", "C",
       ("no -1131C allele", "higher plasma triglycerides (-1131C)", "higher plasma triglycerides (-1131C/C)"),
       "Pennacchio 2001 Science 294:169", "11:116663707", trait="Higher triglyceride allele (APOA5)", has="011"),
    _v("rs1801133", "MTHFR", "diet", "G", "A", "B", ("677CC", "677CT", "677TT"),
       "Frosst 1995 Nat Genet 10:111", "1:11856378", hidden=True, trait="MTHFR 677T", has="011"),
    _v("rs1801131", "MTHFR", "diet", "T", "G", "B", ("1298AA", "1298AC", "1298CC"),
       "van der Put 1998 Am J Hum Genet 62:1044", "1:11854476", hidden=True, trait="MTHFR 1298C", has="011"),

    _v("rs671", "ALDH2", "substances", "G", "A", "B",
       ("normal ALDH2 activity", "ALDH2*1/*2: alcohol flushing; drinking raises oesophageal cancer risk",
        "ALDH2*2/*2: near-absent ALDH2 activity, severe flushing"),
       "Brooks 2009 PLoS Med 6:e50", "12:112241766", levels="022", trait="Alcohol flush (ALDH2*2)", has="011"),
    _v("rs1229984", "ADH1B", "substances", "C", "T", "B",
       ("ADH1B*1/*1 (typical in Europeans)", "one fast ADH1B*2 allele: lower risk of alcohol dependence",
        "ADH1B*2/*2: lower risk of alcohol dependence"),
       "Walters 2018 Nat Neurosci 21:1656", "4:100239319", trait="Fast alcohol metabolism (ADH1B*2)", has="011"),
    _v("rs16969968", "CHRNA5", "substances", "G", "A", "C",
       ("no D398N risk allele", "heavier smoking if a smoker; higher lung cancer risk",
        "heaviest smoking quantity if a smoker; higher lung cancer risk"),
       "Thorgeirsson 2008 Nature 452:638; Saccone 2010 PLoS Genet 6:e1001053", "15:78882925", levels="011",
       trait="Heavy-smoking risk allele (CHRNA5)", has="011"),

    _v("rs4244285", "CYP2C19", "drugs", "G", "A", "A", ("no *2", "*2 (no function)", "*2/*2"),
       "Lee 2022 Clin Pharmacol Ther 112:959", "10:96541616", hidden=True, trait="CYP2C19*2", has="011"),
    _v("rs4986893", "CYP2C19", "drugs", "G", "A", "A", ("no *3", "*3 (no function)", "*3/*3"),
       "Lee 2022 Clin Pharmacol Ther 112:959", "10:96540410", hidden=True, trait="CYP2C19*3", has="011"),
    _v("rs12248560", "CYP2C19", "drugs", "C", "T", "A", ("no *17", "*17 (increased function)", "*17/*17"),
       "Lee 2022 Clin Pharmacol Ther 112:959", "10:96521657", hidden=True, trait="CYP2C19*17", has="011"),
    _v("rs1799853", "CYP2C9", "drugs", "C", "T", "A", ("no *2", "*2 (decreased function)", "*2/*2"),
       "Johnson 2017 Clin Pharmacol Ther 102:397", "10:96702047", hidden=True, trait="CYP2C9*2", has="011"),
    _v("rs1057910", "CYP2C9", "drugs", "A", "C", "A", ("no *3", "*3 (no function)", "*3/*3"),
       "Johnson 2017 Clin Pharmacol Ther 102:397", "10:96741053", hidden=True, trait="CYP2C9*3", has="011"),
    _v("rs9923231", "VKORC1", "drugs", "C", "T", "A",
       ("-1639GG: usual warfarin sensitivity", "-1639GA: increased warfarin sensitivity",
        "-1639AA: high warfarin sensitivity"),
       "Rieder 2005 N Engl J Med 352:2285; Johnson 2017 Clin Pharmacol Ther 102:397 (CPIC)", "16:31107689",
       levels="011", trait="Warfarin sensitivity (VKORC1 -1639A)", has="011"),
    _v("rs4149056", "SLCO1B1", "drugs", "T", "C", "A",
       ("normal SLCO1B1 function", "decreased function: myopathy OR 4.5 on simvastatin 80 mg/day",
        "poor function: myopathy OR 16.9 on simvastatin 80 mg/day"),
       "SEARCH 2008 N Engl J Med 359:789; Cooper-DeHoff 2022 Clin Pharmacol Ther 111:1007", "12:21331549",
       levels="012", trait="Statin myopathy risk (SLCO1B1)", has="011"),
    _v("rs776746", "CYP3A5", "drugs", "T", "C", "A",
       ("*1/*1: CYP3A5 expressor; CPIC raises the tacrolimus starting dose 1.5-2x",
        "*1/*3: CYP3A5 expressor; CPIC raises the tacrolimus starting dose 1.5-2x",
        "*3/*3: non-expressor (typical in Europeans); standard tacrolimus dosing"),
       "Birdwell 2015 Clin Pharmacol Ther 98:19", "7:99270539", levels="110",
       trait="CYP3A5 expressor (tacrolimus dosing)", has="110"),
    _v("rs3892097", "CYP2D6", "drugs", "C", "T", "A",
       ("no *4 (gene deletions/duplications cannot be detected on chips)",
        "one *4 no-function allele: intermediate metabolizer if the other allele is normal",
        "*4/*4: poor metabolizer (codeine, tramadol, tamoxifen, many antidepressants)"),
       "Caudle 2020 Clin Transl Sci 13:116", "22:42524947", levels="012",
       trait="CYP2D6*4 no-function allele", has="011"),
    _v("rs2395029", "HCP5", "drugs", "T", "G", "A",
       ("HLA-B*57:01 tag absent (tag validated in Europeans)",
        "HLA-B*57:01 tag present: abacavir hypersensitivity risk, confirm by HLA typing",
        "HLA-B*57:01 tag present (two copies): abacavir hypersensitivity risk, confirm by HLA typing"),
       "Mallal 2008 N Engl J Med 358:568; Colombo 2008 J Infect Dis 198:864", "6:31431780", levels="022",
       trait="HLA-B*57:01 tag (abacavir hypersensitivity)", has="011"),
    _v("rs3918290", "DPYD", "drugs", "C", "T", "A", ("no *2A", "*2A (no function)", "*2A/*2A"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", "1:97915614", hidden=True, rare=True, trait="DPYD*2A", has="011"),
    _v("rs55886062", "DPYD", "drugs", "A", "C", "A", ("no *13", "*13 (no function)", "*13/*13"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", "1:97981343", hidden=True, rare=True, trait="DPYD*13", has="011"),
    _v("rs67376798", "DPYD", "drugs", "T", "A", "A",
       ("no c.2846A>T", "c.2846A>T (decreased function)", "c.2846A>T homozygous"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", "1:97547947", hidden=True, rare=True,
       trait="DPYD c.2846A>T", has="011"),
    _v("rs56038477", "DPYD", "drugs", "C", "T", "A", ("no HapB3", "HapB3 (decreased function)", "HapB3 homozygous"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", hidden=True, trait="DPYD HapB3", has="011"),
    _v("rs1800462", "TPMT", "drugs", "C", "G", "A", ("no *2", "*2 (no function)", "*2/*2"),
       "Relling 2019 Clin Pharmacol Ther 105:1095", "6:18143955", hidden=True, rare=True, trait="TPMT*2", has="011"),
    _v("rs1800460", "TPMT", "drugs", "C", "T", "A", ("no 460A", "460A (*3B, or *3A with 719G)", "460A homozygous"),
       "Relling 2019 Clin Pharmacol Ther 105:1095", "6:18139228", hidden=True,
       trait="TPMT 460A (*3B or *3A)", has="011"),
    _v("rs1142345", "TPMT", "drugs", "T", "C", "A", ("no 719G", "719G (*3C, or *3A with 460A)", "719G homozygous"),
       "Relling 2019 Clin Pharmacol Ther 105:1095", "6:18130918", hidden=True,
       trait="TPMT 719G (*3C or *3A)", has="011"),
    _v("rs116855232", "NUDT15", "drugs", "C", "T", "A",
       ("no NUDT15*3: normal thiopurine tolerance at this site",
        "NUDT15*3 carrier: intermediate metabolizer; reduced thiopurine starting dose",
        "NUDT15*3/*3: poor metabolizer; thiopurines need drastic dose reduction"),
       "Relling 2019 Clin Pharmacol Ther 105:1095 (CPIC)", levels="022",
       trait="Reduced NUDT15 activity (thiopurines)", has="011"),

    _v("rs429358", "APOE", "health", "T", "C", "A", ("112Cys/Cys", "112Cys/Arg", "112Arg/Arg"),
       "Corder 1993 Science 261:921", "19:45411941", hidden=True,
       trait="APOE e4-defining allele (rs429358 C)", has="011"),
    _v("rs7412", "APOE", "health", "C", "T", "A", ("158Arg/Arg", "158Arg/Cys", "158Cys/Cys"),
       "Corder 1993 Science 261:921", "19:45412079", hidden=True,
       trait="APOE e2-defining allele (rs7412 T)", has="011"),
    _v("rs6025", "F5", "health", "C", "T", "A",
       ("no factor V Leiden", "factor V Leiden heterozygote: venous thrombosis risk ~7-fold",
        "factor V Leiden homozygote: venous thrombosis risk ~80-fold"),
       "Bertina 1994 Nature 369:64; Rosendaal 1995 Blood 85:1504", "1:169519049", levels="022",
       trait="Factor V Leiden (F5)", has="011"),
    _v("rs1799963", "F2", "health", "G", "A", "A",
       ("no prothrombin G20210A", "G20210A heterozygote: venous thrombosis risk 2.8-fold",
        "G20210A homozygote: strongly increased venous thrombosis risk"),
       "Poort 1996 Blood 88:3698", "11:46761055", levels="022", aliases=("i3002432",),
       trait="Prothrombin G20210A (F2)", has="011"),
    _v("rs1800562", "HFE", "health", "G", "A", "A", ("no C282Y", "C282Y", "C282Y/C282Y"),
       "Feder 1996 Nat Genet 13:399", "6:26093141", hidden=True, trait="HFE C282Y", has="011"),
    _v("rs1799945", "HFE", "health", "C", "G", "A", ("no H63D", "H63D", "H63D/H63D"),
       "Feder 1996 Nat Genet 13:399", "6:26091179", hidden=True, trait="HFE H63D", has="011"),
    _v("rs333", "CCR5", "health", "I", "D", "B",
       ("no CCR5-delta32", "CCR5-delta32 carrier; slower HIV-1 progression reported",
        "no functional CCR5: strong resistance to CCR5-tropic (not CXCR4-tropic) HIV-1"),
       "Samson 1996 Nature 382:722; Liu 1996 Cell 86:367; Dean 1996 Science 273:1856", "3:46414947",
       aliases=("i3003626",), trait="CCR5-delta32", has="011"),
    _v("rs63750847", "APP", "health", "C", "T", "B",
       ("no APP A673T", "A673T: protective against Alzheimer's disease and cognitive decline",
        "A673T homozygous: protective"),
       "Jonsson 2012 Nature 488:96", "21:27269932", rare=True, trait="APP A673T (Alzheimer's-protective)", has="011"),
    _v("rs2187668", "HLA-DQA1", "health", "C", "T", "B", ("no DQ2.5 tag", "DQ2.5 tag", "DQ2.5 tag x2"),
       "Monsuur 2008 PLoS One 3:e2270", "6:32605884", hidden=True, trait="HLA-DQ2.5 tag", has="011"),
    _v("rs7454108", "HLA-DQB1", "health", "T", "C", "B", ("no DQ8 tag", "DQ8 tag", "DQ8 tag x2"),
       "Monsuur 2008 PLoS One 3:e2270", "6:32681483", hidden=True, trait="HLA-DQ8 tag", has="011"),
    _v("rs3184504", "SH2B3", "health", "C", "T", "C",
       ("no R262W risk allele", "autoimmune risk allele (celiac disease, type 1 diabetes)",
        "two autoimmune risk alleles"),
       "Hunt 2008 Nat Genet 40:395", "12:111884608", trait="Autoimmune risk allele (SH2B3 R262W)", has="011"),

    _v("rs4680", "COMT", "neuro", "G", "A", "D",
       ("Val/Val: higher COMT activity", "Val/Met: intermediate COMT activity",
        "Met/Met: lower COMT activity; behavioural associations inconsistent"),
       "Lachman 1996 Pharmacogenetics 6:243; Chen 2004 Am J Hum Genet 75:807", "22:19951271",
       trait="Lower COMT activity (Met allele)", has="011"),
    _v("rs6265", "BDNF", "neuro", "C", "T", "D",
       ("Val/Val", "Val/Met: reduced activity-dependent BDNF secretion; no depression association in large samples",
        "Met/Met: reduced activity-dependent BDNF secretion; no depression association in large samples"),
       "Egan 2003 Cell 112:257; Border 2019 Am J Psychiatry 176:376", "11:27679916",
       trait="BDNF Met allele", has="011"),
    _v("rs1006737", "CACNA1C", "neuro", "G", "A", "C",
       ("no CACNA1C risk allele", "bipolar disorder OR 1.18 per allele (negligible absolute effect)",
        "bipolar disorder OR 1.18 per allele, two copies (negligible absolute effect)"),
       "Ferreira 2008 Nat Genet 40:1056", "12:2345295", trait="Bipolar-associated allele (CACNA1C)", has="011"),
    _v("rs10994336", "ANK3", "neuro", "C", "T", "C",
       ("no ANK3 risk allele", "bipolar disorder OR 1.45 per allele (small absolute effect)",
        "bipolar disorder OR 1.45 per allele, two copies (small absolute effect)"),
       "Ferreira 2008 Nat Genet 40:1056", "10:62179812", trait="Bipolar-associated allele (ANK3)", has="011"),

    _v("rs2802292", "FOXO3", "longevity", "T", "G", "C",
       ("no FOXO3 longevity-associated allele", "one longevity-associated allele", "two longevity-associated alleles"),
       "Willcox 2008 Proc Natl Acad Sci USA 105:13987; Flachsbart 2009 Proc Natl Acad Sci USA 106:2700", "6:108908518",
       trait="Longevity-associated allele (FOXO3)", has="011"),
]

DB_BY_ID = {v.rsid: v for v in DB}


class ParseError(Exception):
    pass


@dataclass
class Record:
    chrom: str
    pos: int
    alleles: tuple | None
    raw: str


@dataclass
class Genome:
    path: str
    fmt: str = "unknown"
    build: str = ""
    build_note: str = ""
    sample: str = ""
    markers: int = 0
    nocalls: int = 0
    x_single: int = 0
    x_diploid: int = 0
    x_het: int = 0
    sex: str = "unknown"
    strand: str = "plus"
    by_id: dict = field(default_factory=dict)
    by_pos: dict = field(default_factory=dict)


def iter_text(path):
    """Yield text lines from plain, gzip, bzip2 or zip input (largest zip member)."""
    with open(path, "rb") as fh:
        magic = fh.read(4)
    if magic[:2] == b"\x1f\x8b":
        with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as fh:
            yield from fh
    elif magic[:3] == b"BZh":
        with bz2.open(path, "rt", encoding="utf-8-sig", errors="replace") as fh:
            yield from fh
    elif magic == b"PK\x03\x04":
        with zipfile.ZipFile(path) as zf:
            members = [m for m in zf.infolist() if not m.is_dir()
                       and not m.filename.startswith("__MACOSX") and not os.path.basename(m.filename).startswith(".")]
            if not members:
                raise ParseError("zip archive contains no files")
            member = max(members, key=lambda m: m.file_size)
            with zf.open(member) as raw:
                stream = gzip.GzipFile(fileobj=raw) if member.filename.endswith(".gz") else raw
                yield from io.TextIOWrapper(stream, encoding="utf-8-sig", errors="replace")
    else:
        with open(path, encoding="utf-8-sig", errors="replace") as fh:
            yield from fh


CHROM_ALIASES = {"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT", "PAR": "XY"}


def norm_chrom(c: str) -> str:
    c = c.strip().upper()
    if c.startswith("CHR"):
        c = c[3:]
    return CHROM_ALIASES.get(c, c)


def parse_gt(s: str):
    """Normalise a genotype string to a sorted allele tuple; None for no-calls."""
    s = s.strip().upper().replace("/", "").replace("|", "").replace(" ", "")
    if len(s) in (1, 2) and all(c in BASES for c in s):
        return tuple(sorted(s))
    return None


def _ncol(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


ID_COLS = ("rsid", "snpname", "name", "markername", "snp", "snpid", "id")
CHR_COLS = ("chromosome", "chrom", "chr")
POS_COLS = ("position", "pos", "basepairposition", "physicalposition")
GT_COLS = ("genotype", "result", "gt", "call")
A1_COLS = ("allele1plus", "allele1forward", "allele1")
A2_COLS = ("allele2plus", "allele2forward", "allele2")


def _split_header(s: str):
    if "\t" in s:
        return [f.strip().strip('"') for f in s.split("\t")]
    if "," in s:
        return [f.strip().strip('"') for f in s.split(",")]
    return s.split()


KNOWN_COLS = frozenset(ID_COLS + CHR_COLS + POS_COLS + GT_COLS + A1_COLS + A2_COLS)


def _schema_from_header(s: str):
    cols = [_ncol(f) for f in _split_header(s)]
    if len(cols) < 3 or not (cols[0] in ID_COLS or sum(c in KNOWN_COLS for c in cols) * 2 >= len(cols)):
        return None

    def find(names):
        for n in names:
            if n in cols:
                return cols.index(n)
        return None

    sc = {"id": find(ID_COLS), "chrom": find(CHR_COLS), "pos": find(POS_COLS),
          "gt": find(GT_COLS), "a1": find(A1_COLS), "a2": find(A2_COLS)}
    if sc["id"] is None or (sc["gt"] is None and (sc["a1"] is None or sc["a2"] is None)):
        return None
    mode = "csv" if ("," in s and "\t" not in s) else "ws"
    return sc, mode


def _schema_from_data(s: str):
    mode = "csv" if ("," in s and "\t" not in s) else "ws"
    f = _fields(s, mode)
    if not f or not f[0]:
        return None
    if len(f) == 4:
        return {"id": 0, "chrom": 1, "pos": 2, "gt": 3, "a1": None, "a2": None}, mode
    if len(f) == 5:
        return {"id": 0, "chrom": 1, "pos": 2, "gt": None, "a1": 3, "a2": 4}, mode
    return None


def _fields(s: str, mode: str):
    if mode == "csv":
        return [f.strip() for f in s.replace('"', "").split(",")]
    return s.split()


def _pick(f, i):
    return f[i] if i is not None and i < len(f) else ""


BUILD_PATTERNS = (
    ("37", re.compile(r"build\s*37|grch37|hg19|\bb37\b|249250621")),
    ("36", re.compile(r"build\s*36|ncbi36|hg18|\bb36\b")),
    ("38", re.compile(r"build\s*38|grch38|hg38|\bb38\b|248956422")),
)
VENDORS = (("23andme", "23andMe"), ("ancestrydna", "AncestryDNA"), ("myheritage", "MyHeritage"),
           ("famfinder", "FamilyTreeDNA"), ("familytreedna", "FamilyTreeDNA"), ("living dna", "LivingDNA"),
           ("dna.land", "DNA.Land"), ("genes for good", "Genes for Good"))


def _meta(genome: Genome, meta: list):
    text = "\n".join(meta).lower()
    for build, pat in BUILD_PATTERNS:
        if pat.search(text):
            genome.build, genome.build_note = build, "header"
            break
    for key, name in VENDORS:
        if key in text:
            genome.fmt = name
            break


def _store(genome: Genome, ids, chrom, pos, alleles, raw, want_ids, want_pos):
    genome.markers += 1
    if alleles is None:
        genome.nocalls += 1
    elif chrom == "X" and PAR_X[0] < pos < PAR_X[1]:
        if len(alleles) == 1:
            genome.x_single += 1
        else:
            genome.x_diploid += 1
            genome.x_het += alleles[0] != alleles[1]
    rec = None
    for vid in ids:
        if vid in want_ids:
            rec = rec or Record(chrom, pos, alleles, raw)
            old = genome.by_id.get(vid)
            if old is None or (old.alleles is None and alleles is not None):
                genome.by_id[vid] = rec
    if pos and (chrom, pos) in want_pos:
        genome.by_pos.setdefault((chrom, pos), []).append(rec or Record(chrom, pos, alleles, raw))


def _vcf_alleles(ref: str, alt: str, fmt: str, sample: str):
    keys = fmt.split(":")
    if "GT" not in keys:
        return None
    vals = sample.split(":")
    k = keys.index("GT")
    gt = vals[k] if k < len(vals) else "."
    alleles = [ref.upper()] + alt.upper().split(",")
    out = []
    for i in re.split(r"[/|]", gt):
        if not i.isdigit() or int(i) >= len(alleles):
            return None
        a = alleles[int(i)]
        if a in (".", "*", "") or a.startswith("<"):
            return None
        out.append(a)
    if all(len(a) == 1 and a in "ACGT" for a in out):
        return tuple(sorted(out))
    lens = {len(a) for a in alleles if a not in (".", "*") and not a.startswith("<")}
    if len(lens) == 2:
        hi = max(lens)
        return tuple(sorted("I" if len(a) == hi else "D" for a in out))
    return None


def load_genome(path: str, want_ids: set, want_pos: set) -> Genome:
    genome = Genome(path=path)
    lines = iter_text(path)
    meta, schema, mode, skip_section, vcf = [], None, "ws", False, None
    for n, line in enumerate(lines):
        s = line.strip()
        if not s:
            continue
        if n == 0 and s.startswith("##fileformat=VCF"):
            vcf = True
        if vcf:
            if s.startswith("##"):
                if len(meta) < 5000:
                    meta.append(s)
                continue
            if s.startswith("#CHROM"):
                cols = s.lstrip("#").split("\t")
                genome.sample = cols[9] if len(cols) > 9 else ""
                continue
            f = s.split("\t")
            if len(f) < 10:
                continue
            ids = [x for x in f[2].split(";") if x and x != "."]
            pos = int(f[1]) if f[1].isdigit() else 0
            alleles = _vcf_alleles(f[3], f[4], f[8], f[9])
            _store(genome, ids, norm_chrom(f[0]), pos, alleles, "/".join(alleles or ("-",)), want_ids, want_pos)
            continue
        if s.startswith("[") and s.endswith("]"):
            skip_section = s.lower() != "[data]"
            continue
        if skip_section:
            continue
        if s.startswith("#"):
            if len(meta) < 500:
                meta.append(s)
            if schema is None:
                found = _schema_from_header(s.lstrip("#").strip())
                if found:
                    schema, mode = found
            continue
        if schema is None:
            found = _schema_from_header(s)
            if found:
                schema, mode = found
                continue
            found = _schema_from_data(s)
            if not found:
                raise ParseError(f"unrecognised layout near line {n + 1}: {s[:80]!r}")
            schema, mode = found
        f = _fields(s, mode)
        vid = _pick(f, schema["id"])
        if not vid:
            continue
        p = _pick(f, schema["pos"])
        raw = _pick(f, schema["gt"]) if schema["gt"] is not None else _pick(f, schema["a1"]) + _pick(f, schema["a2"])
        _store(genome, (vid,), norm_chrom(_pick(f, schema["chrom"])), int(p) if p.isdigit() else 0,
               parse_gt(raw), raw, want_ids, want_pos)
    if genome.markers == 0:
        raise ParseError("no genotype records found")
    _meta(genome, meta)
    if vcf:
        genome.fmt = "VCF" + (f" (sample {genome.sample})" if genome.sample else "")
    elif genome.fmt == "unknown" and schema is not None:
        genome.fmt = "generic table (" + ("allele columns" if schema["gt"] is None else "genotype column") + ")"
    return genome


def finalize(genome: Genome, db, forced_build: str = ""):
    checked = agree = 0
    for v in db:
        rec = genome.by_id.get(v.rsid)
        if v.pos and rec is not None and rec.pos:
            checked += 1
            agree += rec.pos == v.pos
    if forced_build:
        genome.build, genome.build_note = forced_build, "--build"
    elif checked >= 5 and agree / checked >= 0.9:
        note = f"{agree}/{checked} database positions match GRCh37"
        if genome.build in ("", "37"):
            genome.build, genome.build_note = "37", (genome.build_note + "; " if genome.build_note else "") + note
        else:
            genome.build_note += f"; but {note}"
    elif checked >= 5 and genome.build == "37":
        genome.build_note += f"; only {agree}/{checked} database positions match GRCh37"
    diploid = genome.x_diploid
    if genome.x_single + diploid >= 50:
        if genome.x_single / (genome.x_single + diploid) >= 0.9:
            genome.sex = "male (hemizygous X calls)"
        elif diploid and genome.x_het / diploid < 0.03:
            genome.sex = "male (no X heterozygosity)"
        elif diploid and genome.x_het / diploid > 0.08:
            genome.sex = "female (X heterozygosity {:.1%})".format(genome.x_het / diploid)
    return genome


@dataclass
class Call:
    v: Variant
    status: str = "absent"
    alleles: tuple = ()
    raw: str = ""
    flipped: bool = False
    via: str = ""
    unverified: bool = False

    def count(self, allele: str):
        if self.status != "ok":
            return None
        k = self.alleles.count(allele)
        return k * 2 if len(self.alleles) == 1 else k

    @property
    def dosage(self):
        return self.count(self.v.alt)

    @property
    def genotype(self) -> str:
        if self.status == "ok":
            return "/".join(self.alleles)
        return {"nocall": "--", "absent": "n/a"}.get(self.status, self.raw or "?")


def orient(v: Variant, alleles):
    """Return (alleles on the plus strand, flipped) or None if the call fits neither strand."""
    allowed = {v.ref, v.alt}
    if set(alleles) <= allowed:
        return alleles, False
    if not v.palindromic:
        flipped = tuple(sorted(a.translate(COMPLEMENT) for a in alleles))
        if set(flipped) <= allowed:
            return flipped, True
    return None


def call_variant(v: Variant, genome: Genome) -> Call:
    rec, via = None, ""
    for key in (v.rsid,) + v.aliases:
        r = genome.by_id.get(key)
        if r is not None and (rec is None or rec.alleles is None):
            rec, via = r, "rsid" if key == v.rsid else key
    if (rec is None or rec.alleles is None) and genome.build == "37" and v.pos:
        hits = [r for r in genome.by_pos.get((v.chrom, v.pos), []) if r.alleles and orient(v, r.alleles)]
        if hits and len({orient(v, r.alleles)[0] for r in hits}) == 1:
            rec, via = hits[0], "position"
    if rec is None:
        return Call(v)
    if rec.alleles is None:
        return Call(v, "nocall", raw=rec.raw, via=via)
    oriented = orient(v, rec.alleles)
    if oriented is None:
        return Call(v, "conflict", (), rec.raw, False, via)
    return Call(v, "ok", oriented[0], rec.raw, oriented[1], via)


class CallSet(dict):
    def n(self, rsid: str, allele: str | None = None):
        c = self.get(rsid)
        if c is None:
            return None
        return c.count(allele or c.v.alt)

    def missing(self, rsids):
        return [r for r in rsids if self.n(r) is None]


@dataclass
class Finding:
    cat: str
    name: str
    title: str
    has: bool | None
    basis: str
    text: str
    ev: str
    src: str
    level: int = 0
    data: dict = field(default_factory=dict)


def _basis(cs, rsids):
    def one(r):
        c = cs[r]
        if c.status == "ok":
            return f"{r} {c.genotype}"
        return f"{r} " + {"nocall": "(no-call)", "absent": "(not in your file)"}.get(c.status, f"({c.raw!r} conflict)")
    found = all(cs[r].status == "ok" for r in rsids)
    return ", ".join(one(r) for r in rsids) + (" in your file" if found else "")


def _na(cat, name, title, basis, ev, src, rsids):
    return Finding(cat, name, title, None, basis, "not assessable: " + ", ".join(rsids) + " absent or no-call",
                   ev, src, -1)


def _untested(miss):
    return f" (not genotyped: {', '.join(miss)})" if miss else ""


APOE_SRC = "Corder 1993 Science 261:921; Farrer 1997 JAMA 278:1349"
APOE_TABLE = {(0, 0): "e3/e3", (0, 1): "e2/e3", (0, 2): "e2/e2", (1, 0): "e3/e4", (1, 1): "e2/e4",
              (2, 0): "e4/e4", (1, 2): "e1/e2", (2, 1): "e1/e4", (2, 2): "e1/e1"}
APOE_TEXT = {
    "e2/e2": ("Alzheimer's disease OR 0.6 vs e3/e3", 0),
    "e2/e3": ("Alzheimer's disease OR 0.6 vs e3/e3", 0),
    "e3/e3": ("most common genotype, reference risk", 0),
    "e2/e4": ("Alzheimer's disease OR 2.6 vs e3/e3 (rare e1/e3 gives the same genotypes)", 1),
    "e3/e4": ("Alzheimer's disease OR 3.2 vs e3/e3", 1),
    "e4/e4": ("Alzheimer's disease OR 14.9 vs e3/e3", 2),
}


def f_apoe(cs):
    ids = ("rs429358", "rs7412")
    title, basis = "APOE e4 (Alzheimer's risk allele)", _basis(cs, ids)
    a, b = cs.n("rs429358"), cs.n("rs7412")
    if a is None or b is None:
        return [_na("health", "APOE", title, basis, "A", APOE_SRC, cs.missing(ids))]
    geno = APOE_TABLE[(a, b)]
    text, level = APOE_TEXT.get(geno, ("implies the rare e1 haplotype; a genotyping error is more likely, verify", 1))
    has = None if geno in ("e1/e2", "e1/e1") else "e4" in geno
    out = [Finding("health", "APOE", title, has, basis, f"{geno}: {text} (Farrer 1997, Caucasian clinical series)",
                   "A", APOE_SRC, level, {"genotype": geno})]
    if "e1" not in geno:
        if "e4" in geno:
            lt = f"{geno}: e4 is under-represented among long-lived individuals"
        elif "e2" in geno:
            lt = f"{geno}: e2 is over-represented among long-lived individuals"
        else:
            lt = "e3/e3: neutral at the strongest longevity locus"
        out.append(Finding("longevity", "APOE longevity", "Longevity-favouring APOE (e2, no e4)",
                           "e2" in geno and "e4" not in geno, basis, lt, "B", "Deelen 2019 Nat Commun 10:3669", 0,
                           {"genotype": geno}))
    return out


IRISPLEX = (
    ("rs12913832", "A", -5.385897557, -2.302942706),
    ("rs1800407", "T", 1.330922472, 0.9785432025),
    ("rs12896399", "T", 0.7636299425, 0.2541024725),
    ("rs16891982", "C", -1.529280117, -0.9342328224),
    ("rs1393350", "A", 0.4340394174, 0.2087410752),
    ("rs12203592", "T", 0.653488393, 0.6457702022),
)
IRISPLEX_ALPHA = (2.575575265, 0.306416526)
IRISPLEX_SRC = "Walsh 2011 Forensic Sci Int Genet 5:170; Liu 2009 Curr Biol 19:R192 (HIrisPlex-S webtool parameters)"
EYE_TITLES = {"blue": "Blue eyes (IrisPlex)", "intermediate": "Green/hazel eyes (IrisPlex)",
              "brown": "Brown eyes (IrisPlex)"}


def irisplex(counts: dict) -> dict:
    lb = IRISPLEX_ALPHA[0] + sum(bb * counts[r] for r, _, bb, _ in IRISPLEX)
    li = IRISPLEX_ALPHA[1] + sum(bi * counts[r] for r, _, _, bi in IRISPLEX)
    eb, ei = math.exp(lb), math.exp(li)
    d = 1.0 + eb + ei
    return {"blue": eb / d, "intermediate": ei / d, "brown": 1.0 / d}


def f_irisplex(cs):
    counts = {r: cs.n(r, a) for r, a, _, _ in IRISPLEX}
    miss = [r for r, k in counts.items() if k is None]
    basis = f"IrisPlex model, {len(IRISPLEX) - len(miss)} of {len(IRISPLEX)} SNPs in your file"
    if miss:
        return [_na("appearance", "IrisPlex", "Eye colour (IrisPlex)", basis, "B", IRISPLEX_SRC, miss)]
    p = irisplex(counts)
    best = max(p, key=p.get)
    called = p[best] >= 0.7
    text = f"blue {p['blue']:.3f}, intermediate {p['intermediate']:.3f}, brown {p['brown']:.3f}"
    text += f" -> {best}" if called else " -> inconclusive (no category reaches 0.7)"
    title = EYE_TITLES[best] if called else "Eye colour (IrisPlex)"
    return [Finding("appearance", "IrisPlex", title, True if called else None, basis, text, "B", IRISPLEX_SRC, 0,
                    {"p": p, "counts": counts})]


MC1R_R = (("rs1805006", "D84E"), ("rs11547464", "R142H"), ("rs1805007", "R151C"), ("rs1805008", "R160W"),
          ("rs1805009", "D294H"))
MC1R_r = (("rs1805005", "V60L"), ("rs2228479", "V92M"), ("rs885479", "R163Q"))


def f_mc1r(cs):
    src = ("Flanagan 2000 Hum Mol Genet 9:2531; Duffy 2010 J Invest Dermatol 130:520; "
           "Raimondi 2008 Int J Cancer 122:2753")
    carried, big, small, miss = [], 0, 0, []
    for group, strong in ((MC1R_R, True), (MC1R_r, False)):
        for rsid, name in group:
            k = cs.n(rsid)
            if k is None:
                miss.append(name)
                continue
            if k:
                carried.append(name if k == 1 else f"{name} x2")
            big, small = (big + k, small) if strong else (big, small + k)
    total = len(MC1R_R) + len(MC1R_r)
    basis = f"{total - len(miss)} of {total} MC1R coding variants in your file"
    if len(miss) == total:
        return [_na("appearance", "MC1R", "Red hair genotype (MC1R)", basis, "B", src,
                    [r for r, _ in MC1R_R + MC1R_r])]
    if big >= 2:
        title, text, level = "Red hair genotype (MC1R R/R)", "red hair likely (recessive; MC1R variants almost " \
            "never share a haplotype, so two R alleles are in trans); fair skin, freckling, higher melanoma risk", 1
    elif big == 1 and small:
        title, text, level = "MC1R R/r red-hair alleles", "red hair possible; more freckling and sun " \
            "sensitivity, higher melanoma risk", 1
    elif big == 1:
        title, text, level = "MC1R red-hair allele carrier (R/+)", "red hair unlikely; more freckling and sun " \
            "sensitivity, higher melanoma risk", 0
    elif small:
        title, text, level = "MC1R weak red-hair alleles (r)", "small pigmentation effect", 0
    else:
        title, text, level = "MC1R red-hair alleles", "no tested MC1R variant", 0
    if carried:
        text += " [" + ", ".join(carried) + "]"
    return [Finding("appearance", "MC1R", title, big + small > 0, basis, text + _untested(miss), "B", src, level,
                    {"R": big, "r": small, "carried": carried})]


def f_tas2r38(cs):
    ids = ("rs713598", "rs1726866", "rs10246939")
    title, basis, src = "PTC/PROP bitter taster (TAS2R38)", _basis(cs, ids), "Kim 2003 Science 299:1221"
    k = [cs.n(r) for r in ids]
    if None in k:
        return [_na("senses", "TAS2R38", title, basis, "B", src, cs.missing(ids))]
    if k[0] == k[1] == k[2]:
        has = k[0] > 0
        text = {2: "PAV/PAV: taster, strongest bitterness perception", 1: "PAV/AVI: taster",
                0: "AVI/AVI: non-taster"}[k[0]]
    else:
        has, text = None, "includes a rare haplotype (AAV, AVV, PVI...); taster status uncertain"
    return [Finding("senses", "TAS2R38", title, has, basis, text, "B", src, 0, {"pav_counts": k})]


def f_hfe(cs):
    src = "Feder 1996 Nat Genet 13:399; Allen 2008 N Engl J Med 358:221"
    ids = ("rs1800562", "rs1799945")
    basis = _basis(cs, ids)
    c, h = cs.n("rs1800562"), cs.n("rs1799945")
    if c is None:
        return [_na("health", "HFE", "HFE haemochromatosis variants", basis, "A", src, ["rs1800562"])]
    h0 = h or 0
    if c == 2:
        title, text, level = "Haemochromatosis genotype (HFE C282Y/C282Y)", "iron-overload disease in 28.4% of " \
            "male and 1.2% of female homozygotes; ferritin and transferrin saturation are the follow-up tests", 2
    elif c == 1 and h0 == 1:
        title, text, level = "HFE compound heterozygote (C282Y/H63D)", "mild iron loading possible; overload " \
            "disease uncommon", 1
    elif c == 1:
        title, text, level = "HFE C282Y carrier", "no iron overload expected", 0
    elif h0 == 2:
        title, text, level = "HFE H63D/H63D", "clinical iron overload uncommon", 0
    elif h0 == 1:
        title, text, level = "HFE H63D carrier", "no iron overload expected", 0
    else:
        title, text, level = "HFE haemochromatosis variants", "no C282Y or H63D", 0
    if c + h0 > 2:
        text += "; unusual combination, verify genotyping"
    text += "" if h is not None else " (H63D not genotyped)"
    return [Finding("health", "HFE", title, c + h0 > 0, basis, text, "A", src, level, {"C282Y": c, "H63D": h})]


def f_mthfr(cs):
    src = "Frosst 1995 Nat Genet 10:111; van der Put 1998 Am J Hum Genet 62:1044; Hickey 2013 Genet Med 15:153"
    ids = ("rs1801133", "rs1801131")
    title, basis = "Reduced MTHFR activity", _basis(cs, ids)
    a, b = cs.n("rs1801133"), cs.n("rs1801131")
    if a is None and b is None:
        return [_na("diet", "MTHFR", title, basis, "B", src, list(ids))]
    a0, b0 = a or 0, b or 0
    table = {(2, 0): "677TT: thermolabile enzyme with reduced activity; homocysteine rises mainly when folate is low",
             (1, 1): "677CT/1298AC compound heterozygote: mildly reduced activity",
             (1, 0): "677CT: mildly reduced activity", (0, 2): "1298CC: mildly reduced activity",
             (0, 1): "1298AC: minimal effect", (0, 0): "677CC/1298AA: reference genotype"}
    text = table.get((a0, b0), "unusual combination (677T and 1298C are rarely in cis); verify genotyping")
    text += "; ACMG advises against MTHFR testing for thrombophilia"
    miss = [n for n, k in (("C677T", a), ("A1298C", b)) if k is None]
    return [Finding("diet", "MTHFR", title, a0 >= 1 or b0 == 2, basis, text + _untested(miss), "B", src, 0,
                    {"677T": a, "1298C": b})]


def f_celiac(cs):
    src = "Monsuur 2008 PLoS One 3:e2270"
    ids = ("rs2187668", "rs7454108")
    title, basis = "Celiac-permissive HLA-DQ (DQ2.5/DQ8)", _basis(cs, ids)
    dq25, dq8 = cs.n("rs2187668"), cs.n("rs7454108")
    if dq25 is None and dq8 is None:
        return [_na("health", "HLA-DQ", title, basis, "B", src, list(ids))]
    tags = [n for n, k in (("DQ2.5", dq25), ("DQ8", dq8)) if k]
    miss = [n for n, k in (("DQ2.5", dq25), ("DQ8", dq8)) if k is None]
    if tags:
        has = True
        text = " and ".join(tags) + " tag present: permissive for celiac disease; most carriers never develop it"
    else:
        has, text = (None if miss else False), "no DQ2.5/DQ8 tag: celiac disease unlikely (DQ2.2 not assessed)"
    return [Finding("health", "HLA-DQ", title, has, basis, text + _untested(miss), "B", src, 0,
                    {"DQ2.5": dq25, "DQ8": dq8})]


def _diplotype(alleles):
    alleles = list(alleles)
    while len(alleles) < 2:
        alleles.insert(0, "*1")
    return "/".join(alleles)


def f_cyp2c19(cs):
    src = "Lee 2022 Clin Pharmacol Ther 112:959 (CPIC)"
    ids = ("rs4244285", "rs4986893", "rs12248560")
    basis = _basis(cs, ids)
    n2, n3, n17 = cs.n("rs4244285"), cs.n("rs4986893"), cs.n("rs12248560")
    if n2 is None or n17 is None:
        return [_na("drugs", "CYP2C19", "Altered CYP2C19 function", basis, "A", src,
                    cs.missing(("rs4244285", "rs12248560")))]
    n3 = n3 or 0
    lof, gof = n2 + n3, n17
    if lof + gof > 2:
        return [Finding("drugs", "CYP2C19", "Altered CYP2C19 function", None, basis,
                        "more than two variant alleles; cannot assign a diplotype", "A", src, 1)]
    dip = _diplotype(["*2"] * n2 + ["*3"] * n3 + ["*17"] * n17)
    if lof == 2:
        phen, level = "poor metabolizer", 2
    elif lof == 1:
        phen, level = "intermediate metabolizer", 1
    elif gof == 2:
        phen, level = "ultrarapid metabolizer", 1
    elif gof == 1:
        phen, level = "rapid metabolizer", 0
    else:
        phen, level = "normal metabolizer", 0
    title = ("Reduced CYP2C19 function" if lof else "Increased CYP2C19 function" if gof
             else "Altered CYP2C19 function")
    text = f"{dip}: {phen}; relevant to clopidogrel, PPIs, citalopram/escitalopram, voriconazole"
    text += _untested(["*3"] if cs.n("rs4986893") is None else [])
    return [Finding("drugs", "CYP2C19", title, bool(lof or gof), basis, text, "A", src, level,
                    {"diplotype": dip, "phenotype": phen})]


def f_cyp2c9(cs):
    src = "Johnson 2017 Clin Pharmacol Ther 102:397 (CPIC)"
    ids = ("rs1799853", "rs1057910")
    title, basis = "Reduced CYP2C9 function", _basis(cs, ids)
    n2, n3 = cs.n("rs1799853"), cs.n("rs1057910")
    if n2 is None or n3 is None:
        return [_na("drugs", "CYP2C9", title, basis, "A", src, cs.missing(ids))]
    if n2 + n3 > 2:
        return [Finding("drugs", "CYP2C9", title, True, basis, "more than two variant alleles; verify", "A", src, 1)]
    score = 2 - 0.5 * n2 - 1.0 * n3
    phen = "normal" if score == 2 else "intermediate" if score >= 1 else "poor"
    dip = _diplotype(["*2"] * n2 + ["*3"] * n3)
    text = f"{dip}: {phen} metabolizer (activity score {score:g}); relevant to warfarin, phenytoin, several NSAIDs"
    return [Finding("drugs", "CYP2C9", title, score < 2, basis, text, "A", src,
                    {"normal": 0, "intermediate": 1, "poor": 2}[phen], {"diplotype": dip, "activity_score": score})]


DPYD = (("rs3918290", "*2A", 0.0), ("rs55886062", "*13", 0.0), ("rs67376798", "c.2846A>T", 0.5),
        ("rs56038477", "HapB3", 0.5))


def f_dpyd(cs):
    src = "Amstutz 2018 Clin Pharmacol Ther 103:210 (CPIC)"
    title = "Reduced DPYD activity (fluoropyrimidine toxicity)"
    found, miss, total, score = [], [], 0, 2.0
    for rsid, name, value in DPYD:
        k = cs.n(rsid)
        if k is None:
            miss.append(name)
            continue
        if k:
            found.append(name if k == 1 else f"{name} x2")
            total += k
            score -= k * (1 - value)
    basis = f"{len(DPYD) - len(miss)} of {len(DPYD)} CPIC DPYD variants in your file"
    if len(miss) == len(DPYD):
        return [_na("drugs", "DPYD", title, basis, "A", src, [r for r, _, _ in DPYD])]
    if total > 2:
        return [Finding("drugs", "DPYD", title, True, basis, "more than two variant alleles; cannot score", "A",
                        src, 2)]
    if score >= 2:
        text, level = "no tested decreased/no-function variant: normal activity assumed", 0
    elif score >= 1:
        text, level = f"activity score {score:g}: intermediate metabolizer; CPIC advises a reduced " \
                      f"fluoropyrimidine (5-FU, capecitabine) starting dose", 2
    else:
        text, level = f"activity score {score:g}: poor metabolizer; CPIC advises avoiding fluoropyrimidines", 2
    if found:
        text += " [" + ", ".join(found) + "; rare variant, confirm clinically]"
    return [Finding("drugs", "DPYD", title, score < 2, basis, text + _untested(miss), "A", src, level,
                    {"activity_score": score, "variants": found})]


def f_tpmt(cs):
    src = "Relling 2019 Clin Pharmacol Ther 105:1095 (CPIC)"
    ids = ("rs1800462", "rs1800460", "rs1142345")
    title, basis = "Reduced TPMT activity (thiopurines)", _basis(cs, ids)
    t2, t3b, t3c = cs.n("rs1800462"), cs.n("rs1800460"), cs.n("rs1142345")
    if t3b is None or t3c is None:
        return [_na("drugs", "TPMT", title, basis, "A", src, cs.missing(ids[1:]))]
    a3 = min(t3b, t3c)
    nb, nc, n2 = t3b - a3, t3c - a3, t2 or 0
    total = n2 + a3 + nb + nc
    dip = _diplotype(["*2"] * n2 + ["*3A"] * a3 + ["*3B"] * nb + ["*3C"] * nc)
    if total > 2:
        text, level = "more than two variant alleles; verify", 2
    elif total == 0:
        text, level = "*1/*1: normal metabolizer", 0
    else:
        phen = "intermediate" if total == 1 else "poor"
        note = "; *3B/*3C in trans (poor) is rare but not excluded" if (a3 == 1 and total == 1) else ""
        text, level = f"{dip}: {phen} metabolizer{note}", 2
    text += "; relevant to azathioprine, mercaptopurine, thioguanine"
    text += "" if t2 is not None else " (*2 not genotyped)"
    return [Finding("drugs", "TPMT", title, total > 0, basis, text, "A", src, level, {"tpmt": dip})]


FINDERS = (f_irisplex, f_mc1r, f_tas2r38, f_mthfr, f_cyp2c19, f_cyp2c9, f_dpyd, f_tpmt, f_apoe, f_hfe, f_celiac)


def resolve_strand(genome: Genome, cs: CallSet):
    """Complement A/T and C/G calls when the file is consistently minus-strand; flag them when mixed."""
    plain = [c for c in cs.values() if c.status == "ok" and not c.v.palindromic]
    flips = sum(c.flipped for c in plain)
    if not flips:
        return
    frac = flips / len(plain)
    if flips >= 10 and frac >= 0.9:
        genome.strand = "minus"
    elif flips >= 3 and frac >= 0.1:
        genome.strand = "mixed"
    else:
        return
    for c in cs.values():
        if c.status == "ok" and c.v.palindromic:
            if genome.strand == "minus":
                c.alleles, c.flipped = tuple(sorted(a.translate(COMPLEMENT) for a in c.alleles)), True
            else:
                c.unverified = True


def analyse(genome: Genome, db=DB):
    cs = CallSet((v.rsid, call_variant(v, genome)) for v in db)
    resolve_strand(genome, cs)
    findings = [f for finder in FINDERS for f in finder(cs)]
    return cs, findings


def _col(cols: dict, *names):
    for n in names:
        if n in cols:
            return cols[n]
    return None


class PGS:
    """PGS Catalog scoring file (format 1.0/2.0, original or harmonized)."""

    def __init__(self, path: str):
        self.path = path
        self.meta, self.cols = {}, {}
        it = iter_text(path)
        for line in it:
            s = line.rstrip("\r\n")
            if not s:
                continue
            if s.startswith("#"):
                key, sep, val = s.lstrip("#").partition("=")
                if sep:
                    self.meta[key.strip()] = val.strip()
                continue
            self.cols = {c.strip(): i for i, c in enumerate(s.split("\t"))}
            break
        it.close()
        c = self.cols
        self.i_rs = [i for i in (_col(c, "hm_rsID"), _col(c, "rsID", "rsid", "SNP")) if i is not None]
        self.i_ea = _col(c, "effect_allele")
        self.i_oa = _col(c, "other_allele", "hm_inferOtherAllele", "reference_allele")
        self.i_w = _col(c, "effect_weight")
        self.i_or = _col(c, "OR", "HR")
        self.i_af = _col(c, "allelefrequency_effect")
        self.i_dos = [_col(c, f"dosage_{k}_weight") for k in range(3)]
        hm = _col(c, "hm_pos") is not None
        self.i_chr = _col(c, "hm_chr") if hm else _col(c, "chr_name")
        self.i_pos = _col(c, "hm_pos") if hm else _col(c, "chr_position")
        build = self.meta.get("HmPOS_build" if hm else "genome_build", "")
        self.build = next((b for b in ("36", "37", "38") if b in build), "")
        wt = self.meta.get("weight_type", "").strip().lower()
        self.log_weights = self.i_w is not None and wt in ("or", "hr", "or/hr", "odds ratio", "hazard ratio")
        if self.i_ea is None or (self.i_w is None and self.i_or is None and None in self.i_dos):
            raise ParseError(f"{path}: not a PGS Catalog scoring file (needs effect_allele and effect_weight)")

    def rows(self):
        header_seen = False
        for line in iter_text(self.path):
            s = line.rstrip("\r\n")
            if not s or s.startswith("#"):
                continue
            if not header_seen:
                header_seen = True
                continue
            yield s.split("\t")

    @staticmethod
    def _cell(f, i):
        return f[i].strip() if i is not None and i < len(f) else ""

    def rsid(self, f):
        for i in self.i_rs:
            v = self._cell(f, i)
            if v and v != "NA":
                return v
        return ""

    def locus(self, f):
        p = self._cell(f, self.i_pos)
        return (norm_chrom(self._cell(f, self.i_chr)), int(p)) if p.isdigit() else None

    def wants(self):
        ids, pos = set(), set()
        for f in self.rows():
            r = self.rsid(f)
            if r:
                ids.add(r)
            loc = self.locus(f)
            if loc:
                pos.add(loc)
        return ids, pos

    def score(self, genome: Genome) -> dict:
        use_pos = bool(self.build) and self.build == genome.build
        res = dict(id=self.meta.get("pgs_id", os.path.basename(self.path)), name=self.meta.get("pgs_name", ""),
                   trait=self.meta.get("trait_reported", ""), variants=0, matched=0, missing=0, flipped=0,
                   conflicts=0, skipped=0, score=0.0)
        mean = var = imputed = 0.0
        have_freq = True
        for f in self.rows():
            res["variants"] += 1
            eff = self._cell(f, self.i_ea).upper()
            oth = self._cell(f, self.i_oa).upper()
            dos = [self._cell(f, i) for i in self.i_dos]
            if any(dos):
                w = None
                dw = [float(x) if x else 0.0 for x in dos]
                have_freq = False
            else:
                raw_w = self._cell(f, self.i_w) if self.i_w is not None else self._cell(f, self.i_or)
                try:
                    w = float(raw_w)
                    if self.log_weights or self.i_w is None:
                        w = math.log(w)
                except (ValueError, TypeError):
                    res["skipped"] += 1
                    continue
            af = self._cell(f, self.i_af)
            try:
                p = float(af)
            except ValueError:
                p, have_freq = None, False
            if w is not None and p is not None:
                mean += 2 * p * w
                var += 2 * p * (1 - p) * w * w
            if len(eff) != 1 or eff not in "ACGT":
                res["skipped"] += 1
                if w is not None and p is not None:
                    imputed += 2 * p * w
                continue
            rec = genome.by_id.get(self.rsid(f))
            if (rec is None or rec.alleles is None) and use_pos:
                loc = self.locus(f)
                hits = [r for r in genome.by_pos.get(loc, []) if r.alleles is not None] if loc else []
                rec = hits[0] if len(hits) == 1 else rec
            obs = rec.alleles if rec is not None else None
            if obs is None or not set(obs) <= set("ACGT"):
                res["missing"] += 1
                if w is not None and p is not None:
                    imputed += 2 * p * w
                continue
            allowed = {eff, oth} if len(oth) == 1 and oth in "ACGT" else None
            if allowed is None or set(obs) <= allowed:
                al = obs
            elif {eff, oth} not in ({"A", "T"}, {"C", "G"}) and set(a.translate(COMPLEMENT) for a in obs) <= allowed:
                al = tuple(a.translate(COMPLEMENT) for a in obs)
                res["flipped"] += 1
            else:
                res["conflicts"] += 1
                if w is not None and p is not None:
                    imputed += 2 * p * w
                continue
            k = al.count(eff) * (2 if len(al) == 1 else 1)
            contrib = dw[k] if w is None else w * k
            res["score"] += contrib
            imputed += contrib
            res["matched"] += 1
        if have_freq and res["variants"] and var > 0:
            z = (imputed - mean) / math.sqrt(var)
            res.update(imputed_score=imputed, expected=mean, sd=math.sqrt(var), z=z,
                       percentile=50 * (1 + math.erf(z / math.sqrt(2))))
        return res


EV_ORDER = "ABCD"
LEVEL_STYLE = {-1: "2", 0: "", 1: "33", 2: "31;1"}


class Style:
    def __init__(self, on: bool):
        self.on = on

    def __call__(self, codes: str, s: str) -> str:
        return f"\033[{codes}m{s}\033[0m" if self.on and codes else s


def _ev_ok(ev: str, minimum: str) -> bool:
    return EV_ORDER.index(ev) <= EV_ORDER.index(minimum)


EV_WORDS = {"A": "clinical grade", "B": "strong", "C": "small effect", "D": "weak"}


def call_text(c: Call):
    """Return (interpretation, level, technical notes) for one variant call."""
    v, notes = c.v, []
    if c.status == "ok":
        d = c.dosage
        text, level = v.text[d], int(v.levels[d])
        if c.flipped:
            notes.append("strand-flipped")
        if c.unverified:
            notes.append("A/T or C/G SNP in a mixed-strand file, orientation unverifiable")
            level = max(level, 1)
        if c.via == "position":
            notes.append("matched by GRCh37 position")
        elif c.via not in ("", "rsid"):
            notes.append(f"reported as {c.via}")
        if v.rare and d:
            text += "; rare variant, confirm clinically"
            level = max(level, 1)
        return text, level, notes
    if c.status == "conflict":
        return f"genotype {c.raw!r} fits neither strand of {v.ref}/{v.alt}; check file build/format", 1, notes
    return ("no call" if c.status == "nocall" else "not on this chip"), -1, notes


def trait_allele(v: Variant) -> str:
    return v.ref if v.has == "110" else v.alt


def copies(c: Call) -> str:
    """What the file holds for this SNP and how many copies of the trait allele that is."""
    v = c.v
    if c.status == "absent":
        return f"{v.rsid} not in your file"
    if c.status == "nocall":
        return f"{v.rsid} in your file as a no-call ({c.raw or '--'})"
    if c.status == "conflict":
        return f"{v.rsid} {c.raw!r} in your file; alleles do not match {v.ref}/{v.alt}"
    allele = trait_allele(v)
    k = c.alleles.count(allele)
    note = f"{k} {'copy' if k == 1 else 'copies'} of trait allele {allele}"
    if len(c.alleles) == 1:
        note += " (hemizygous)"
    if v.has == "001":
        note += ", 2 needed"
    return f"{v.rsid} {c.genotype} in your file; {note}"


def sentence(text: str) -> str:
    if re.match(r"[a-z]{2,}\b", text):
        text = text[0].upper() + text[1:]
    return text if text.endswith((".", "?", "!")) else text + "."


def _wrap(text: str, width: int, indent: int):
    return textwrap.wrap(text, max(width - indent, 30)) or [""]


def _luminance(rgb) -> float:
    def lin(c):
        c /= 255
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (lin(x) for x in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


class Palette:
    """Gives every result title a colour no other title in the report uses. Hues step by the golden ratio so
    neighbouring titles contrast; relative luminance is held between 0.18 and 0.55 so titles stay readable on
    dark and light backgrounds. 24-bit when the terminal advertises it, else the xterm 256-colour cube."""

    GOLDEN = 0.6180339887498949
    LEVELS = (0, 95, 135, 175, 215, 255)

    def __init__(self, truecolor: bool):
        self.truecolor, self.used, self.n = truecolor, set(), 0
        cube = [(r, g, b) for r in range(6) for g in range(6) for b in range(6)
                if max(r, g, b) - min(r, g, b) >= 2
                and 0.18 <= _luminance((self.LEVELS[r], self.LEVELS[g], self.LEVELS[b])) <= 0.55]
        self.cube = sorted(cube, key=lambda c: colorsys.rgb_to_hls(*(x / 5 for x in c))[0])

    @staticmethod
    def _hls(hue: float, light: float):
        return tuple(round(x * 255) for x in colorsys.hls_to_rgb(hue, light, 0.75))

    def _rgb(self, hue: float, target: float):
        lo, hi = 0.0, 1.0
        for _ in range(30):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if _luminance(self._hls(hue, mid)) < target else (lo, mid)
        return self._hls(hue, (lo + hi) / 2)

    def take(self) -> str:
        while True:
            self.n += 1
            hue = (self.n * self.GOLDEN) % 1.0
            if not self.truecolor and len(self.used) < len(self.cube):
                start = int(hue * len(self.cube))
                c = next(self.cube[(start + j) % len(self.cube)] for j in range(len(self.cube))
                         if self.cube[(start + j) % len(self.cube)] not in self.used)
                self.used.add(c)
                return f"38;5;{16 + 36 * c[0] + 6 * c[1] + c[2]}"
            rgb = self._rgb(hue, (0.26, 0.34, 0.42)[self.n % 3])
            if rgb not in self.used:
                self.used.add(rgb)
                return "38;2;{};{};{}".format(*rgb)


def truecolor_supported() -> bool:
    return os.environ.get("COLORTERM", "").lower() in ("truecolor", "24bit") or "WT_SESSION" in os.environ


@dataclass
class Entry:
    title: str
    has: bool | None
    cat: str
    text: str
    level: int
    meta: str
    src: str
    in_file: bool


def entries(cs, findings, category=None, min_evidence="D", show_all=False):
    out = []
    for cat in CATEGORIES:
        if category and cat not in category:
            continue
        for f in findings:
            if f.cat == cat and _ev_ok(f.ev, min_evidence):
                meta = f"{f.basis}; evidence {f.ev} ({EV_WORDS[f.ev]})"
                out.append(Entry(f.title, f.has, cat, f.text, f.level, meta, f.src, f.level >= 0))
        for c in cs.values():
            v = c.v
            if v.cat != cat or not _ev_ok(v.ev, min_evidence) or (v.hidden and not show_all and c.status != "conflict"):
                continue
            text, level, notes = call_text(c)
            meta = "; ".join([copies(c), f"evidence {v.ev} ({EV_WORDS[v.ev]})"] + notes)
            has = v.has[c.dosage] == "1" if c.status == "ok" else None
            out.append(Entry(v.trait, has, cat, text, level, meta, v.src, c.status in ("ok", "conflict")))
    return out


SECTIONS = (("TRAITS PRESENT", lambda e: e.in_file and e.has is True),
            ("FOUND IN YOUR FILE, TRAIT NOT PRESENT", lambda e: e.in_file and e.has is False),
            ("FOUND IN YOUR FILE, INCONCLUSIVE", lambda e: e.in_file and e.has is None),
            ("NOT IN YOUR FILE OR NO-CALL", lambda e: not e.in_file))


def render_text(genome, cs, findings, pgs_results, args, out):
    st = Style(args.color)
    palette = Palette(getattr(args, "truecolor", False))
    width = min(shutil.get_terminal_size((100, 24)).columns, 100) if args.color else 100
    calls = list(cs.values())
    n = {s: sum(c.status == s for c in calls) for s in ("ok", "nocall", "absent", "conflict")}
    flips = sum(c.flipped for c in calls)
    build = {"36": "NCBI36", "37": "GRCh37", "38": "GRCh38"}.get(genome.build, "build unknown")
    items = entries(cs, findings, args.category, args.min_evidence, args.all)
    groups = [(label, [e for e in items if test(e)]) for label, test in SECTIONS]
    shown = groups if args.all else [groups[0], groups[2]]
    rate = genome.nocalls / genome.markers if genome.markers else 0.0
    sex = genome.sex.split(" (")[0]
    w = out.write
    w(st("1", f"genetic-trait-detector") + "\n")
    w(f"file      {genome.path}\n")
    w(f"data      {genome.fmt}, {build}, {'sex unknown' if sex == 'unknown' else sex}, {genome.markers:,} markers "
      f"({rate:.1%} no-call)\n")
    w(f"coverage  {n['ok']} of {len(calls)} database SNPs found in your file, {flips} strand-flipped, "
      f"{n['conflict']} conflicts\n")
    present, absent, unclear, missing = (g for _, g in groups)
    summary = f"{len(present)} trait present" + (f", {len(unclear)} inconclusive" if unclear else "")
    rest = f"{len(absent)} not present, {len(missing)} not in your file"
    w(f"results   {summary}; {rest}" + ("" if args.all else " (-a lists them)") + "\n")
    if genome.strand == "minus":
        w(st("33", "warning   file reports the minus strand; all calls, including A/T and C/G SNPs, were "
                   "complemented") + "\n")
    elif genome.strand == "mixed":
        w(st("33", "warning   file mixes strands; A/T and C/G SNPs cannot be oriented and are flagged") + "\n")
    if "only" in genome.build_note or "but" in genome.build_note:
        w(st("33", f"warning   build check: {genome.build_note}") + "\n")
    if genome.build != "37":
        w(f"note      {build}: rsID matching only, no position fallback\n")
    for label, group in shown:
        if not group:
            continue
        w("\n" + st("1", f"{label} ({len(group)})") + "\n")
        for cat, (cat_label, _) in CATEGORIES.items():
            block = [e for e in group if e.cat == cat]
            if not block:
                continue
            w("\n" + st("1", cat_label) + "\n")
            for e in block:
                w("\n  " + st("1;" + palette.take(), e.title) + "\n")
                for line in _wrap(sentence(e.text), width, 4):
                    w("    " + st(LEVEL_STYLE[e.level], line) + "\n")
                for line in _wrap(e.meta, width, 4):
                    w("    " + st("2", line) + "\n")
                if not args.no_sources:
                    for line in _wrap(e.src, width, 4):
                        w("    " + st("2", line) + "\n")
    for r in pgs_results:
        w("\n" + st("1", f"POLYGENIC SCORE {r['id']}") + (f" ({r['name']}; {r['trait']})" if r["trait"] else "") + "\n")
        w(f"  {r['variants']:,} variants: {r['matched']:,} matched, {r['missing']:,} missing, {r['flipped']:,} "
          f"strand-flipped, {r['conflicts']:,} conflicts, {r['skipped']:,} skipped (non-SNV or unusable weight)\n")
        w(f"  raw score {r['score']:.6g}\n")
        if "z" in r:
            w(f"  mean-imputed score {r['imputed_score']:.6g} vs expected {r['expected']:.6g} (SD {r['sd']:.4g}): "
              f"z {r['z']:+.2f}, percentile ~{r['percentile']:.0f} (HWE, independent SNPs)\n")
        else:
            w("  no percentile: scoring file lacks effect-allele frequencies for every variant\n")
    w("\nevidence  " + EVIDENCE_SHORT + "\n")


def lookup(genome, cs, rsids, args, out):
    """Show exactly what the file holds for each requested rsID and how the database reads it."""
    st = Style(args.color)
    for rsid in rsids:
        rec = genome.by_id.get(rsid)
        out.write("\n" + st("1", rsid) + "\n")
        if rec is None:
            out.write("  not in your file\n")
        else:
            gt = "/".join(rec.alleles) if rec.alleles else "no-call"
            out.write(f"  in your file: chromosome {rec.chrom or '?'}, position {rec.pos or '?'}, raw {rec.raw!r}, "
                      f"parsed {gt}\n")
        c = cs.get(rsid)
        if c is None:
            out.write("  not in the trait database\n")
            continue
        text, level, notes = call_text(c)
        v = c.v
        verdict = {True: "trait present", False: "trait not present"}.get(
            v.has[c.dosage] == "1" if c.status == "ok" else None, "not determinable")
        out.write(f"  database: {v.trait}; {v.gene} {v.ref}>{v.alt}, trait allele {trait_allele(v)}\n")
        out.write(f"  result: {copies(c)} -> {verdict}\n")
        out.write(f"  {sentence(text)}" + (f" [{'; '.join(notes)}]" if notes else "") + "\n")
        out.write("  " + st("2", v.src) + "\n")


def _variant_json(c: Call):
    text, level, notes = call_text(c)
    v = c.v
    return {"rsid": v.rsid, "gene": v.gene, "category": v.cat, "trait": v.trait,
            "has": (v.has[c.dosage] == "1") if c.status == "ok" else None, "trait_allele": trait_allele(v),
            "trait_allele_copies": c.alleles.count(trait_allele(v)) if c.status == "ok" else None,
            "notes": notes, "ref": v.ref, "effect_allele": v.alt,
            "status": c.status, "genotype": c.genotype if c.status == "ok" else None, "raw": c.raw,
            "dosage": c.dosage, "flipped": c.flipped, "matched_by": c.via or None, "evidence": v.ev,
            "level": level, "interpretation": text, "component": v.hidden, "source": v.src,
            "snpedia": f"https://www.snpedia.com/index.php/{v.rsid}"}


def render_json(genome, cs, findings, pgs_results, args, out):
    items = entries(cs, findings)
    doc = {"tool": "genetic-trait-detector",
           "file": {"path": genome.path, "format": genome.fmt, "build": genome.build or None,
                    "build_note": genome.build_note, "sex": genome.sex, "strand": genome.strand,
                    "markers": genome.markers, "nocalls": genome.nocalls},
           "summary": {"present": [e.title for e in items if e.has is True],
                       "absent": [e.title for e in items if e.has is False],
                       "not_determinable": [e.title for e in items if e.has is None]},
           "variants": [_variant_json(c) for c in cs.values()],
           "findings": [{"name": f.name, "title": f.title, "has": f.has, "category": f.cat, "basis": f.basis,
                         "result": f.text, "evidence": f.ev, "level": f.level, "source": f.src, "data": f.data}
                        for f in findings],
           "pgs": pgs_results, "evidence_scale": EVIDENCE}
    json.dump(doc, out, indent=2)
    out.write("\n")


def render_tsv(genome, cs, findings, pgs_results, args, out):
    cols = ("kind", "category", "id", "gene", "trait", "has", "genotype", "effect_allele", "dosage", "status",
            "flipped", "evidence", "interpretation", "source")
    out.write("\t".join(cols) + "\n")
    tri = {True: "yes", False: "no", None: ""}
    for c in cs.values():
        j = _variant_json(c)
        row = ("variant", j["category"], j["rsid"], j["gene"], j["trait"], tri[j["has"]], j["genotype"] or "",
               j["effect_allele"], "" if j["dosage"] is None else j["dosage"], j["status"], int(j["flipped"]),
               j["evidence"], j["interpretation"], j["source"])
        out.write("\t".join(map(str, row)) + "\n")
    for f in findings:
        row = ("finding", f.cat, f.name, "", f.title, tri[f.has], f.basis, "", "", "ok" if f.has is not None else "na",
               0, f.ev, f.text, f.src)
        out.write("\t".join(map(str, row)) + "\n")
    for r in pgs_results:
        row = ("pgs", "", r["id"], "", r["trait"], "", "", "", "", "ok", 0, "", f"score {r['score']:.6g}; "
               f"matched {r['matched']}/{r['variants']}" + (f"; z {r['z']:+.3f}" if "z" in r else ""), "")
        out.write("\t".join(map(str, row)) + "\n")


def list_db(args, out=sys.stdout):
    st = Style(args.color)
    for cat, (label, color) in CATEGORIES.items():
        if args.category and cat not in args.category:
            continue
        w = [v for v in DB if v.cat == cat and _ev_ok(v.ev, args.min_evidence)]
        if not w:
            continue
        out.write(st(color + ";1", label) + "\n")
        for v in w:
            loc = f"{v.chrom}:{v.pos}" if v.pos else (v.chrom or "-")
            out.write(f"  {v.rsid:<11} {v.gene:<9} {v.ref}>{v.alt}  [{v.ev}] {loc:<13} {v.trait}"
                      f"{' (component)' if v.hidden else ''}  {st('2', v.src)}\n")
    out.write("\ncompound calls: " + ", ".join(f.__name__[2:] for f in FINDERS) + "\n")
    out.write("alleles: GRCh37 plus strand, ref>effect\n")


def build_parser():
    p = argparse.ArgumentParser(prog="genetic-trait-detector",
                                description="Offline trait and variant report for consumer raw DNA files.")
    p.add_argument("file", nargs="?", help="raw genotype file (.txt/.csv/.tsv/.vcf, optionally .gz/.bz2/.zip)")
    p.add_argument("-f", "--format", choices=("text", "json", "tsv"), default="text")
    p.add_argument("-c", "--category", action="append", choices=list(CATEGORIES),
                   help="limit output to a category (repeatable)")
    p.add_argument("-e", "--min-evidence", choices=tuple(EV_ORDER), default="D",
                   help="hide entries weaker than this evidence tier (A strongest)")
    p.add_argument("-a", "--all", action="store_true",
                   help="also list traits not present, SNPs not in your file, and component SNPs")
    p.add_argument("--snp", action="append", default=[], metavar="RSID",
                   help="show exactly what the file holds for an rsID and how it is interpreted (repeatable)")
    p.add_argument("--pgs", action="append", default=[], metavar="FILE",
                   help="apply a PGS Catalog scoring file (repeatable)")
    p.add_argument("--build", choices=("36", "37", "38"), help="override the detected genome build")
    p.add_argument("--no-color", action="store_true", help="disable colour (NO_COLOR is also honoured)")
    p.add_argument("--no-sources", action="store_true", help="omit literature sources in text output")
    p.add_argument("--list", action="store_true", help="print the variant database and exit")
    p.add_argument("--selftest", action="store_true", help="run built-in tests and exit")
    p.add_argument("--version", action="version", version=f"%(prog)s")
    return p


def wanted(pgs_list):
    ids = {v.rsid for v in DB} | {a for v in DB for a in v.aliases}
    pos = {(v.chrom, v.pos) for v in DB if v.pos}
    for s in pgs_list:
        i, p = s.wants()
        ids |= i
        pos |= p
    return ids, pos


def run(path, args, out):
    pgs_list = [PGS(p) for p in args.pgs]
    ids, pos = wanted(pgs_list)
    snps = [s.strip() for s in getattr(args, "snp", []) if s.strip()]
    genome = finalize(load_genome(path, ids | set(snps), pos), DB, args.build or "")
    cs, findings = analyse(genome)
    if snps:
        lookup(genome, cs, snps, args, out)
        return genome, cs, findings, []
    pgs_results = [s.score(genome) for s in pgs_list]
    {"text": render_text, "json": render_json, "tsv": render_tsv}[args.format](
        genome, cs, findings, pgs_results, args, out)
    return genome, cs, findings, pgs_results


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        sys.stdout.reconfigure(errors="replace")
    except (AttributeError, ValueError):
        pass
    if args.selftest:
        return selftest()
    args.color = (args.format == "text" and not args.no_color and not os.environ.get("NO_COLOR")
                  and sys.stdout.isatty())
    args.truecolor = truecolor_supported()
    if args.color and os.name == "nt":
        try:
            import colorama
            colorama.just_fix_windows_console()
        except Exception:
            os.system("")
    if args.list:
        try:
            list_db(args)
        except BrokenPipeError:
            os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0
    path = args.file
    if not path:
        if not sys.stdin.isatty():
            print("error: no input file given", file=sys.stderr)
            return 2
        path = input("Enter DNA data file location: ").strip().strip("'\"")
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        print(f"error: file not found: {path}", file=sys.stderr)
        return 2
    try:
        run(path, args, sys.stdout)
    except BrokenPipeError:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 1
    except (OSError, EOFError, ParseError, zipfile.BadZipFile, ValueError) as e:
        print(f"error: {path}: {e}", file=sys.stderr)
        return 2
    return 0


ST_TRUTH = {
    "rs12913832": "AG", "rs1800407": "CC", "rs12896399": "GG", "rs16891982": "GG", "rs1393350": "AG",
    "rs12203592": "CT", "rs429358": "CT", "rs7412": "CC", "rs11547464": "AG", "rs1805007": "CT",
    "rs713598": "CG", "rs1726866": "AG", "rs10246939": "CT", "rs4244285": "AG", "rs4986893": "GG",
    "rs12248560": "CT", "rs1800460": "CT", "rs1142345": "CT", "rs1800462": "CC", "rs3918290": "CT",
    "rs55886062": "AA", "rs67376798": "TT", "rs6025": "CT", "rs1800562": "AA", "rs1799945": "CC",
    "rs17822931": "TT",
}
ST_CONFLICT = {"rs1815739": "AC"}
ST_MINUS = {"rs4988235": "CT"}
ST_ALIAS = {"i3003626": ("3", 46414947, "DI"), "i3002432": ("11", 46761055, "AG")}
ST_POSONLY = ("vg0001", "1", 11856378, "AA")


def _st_rows(sex="male"):
    rows = []
    for table in (ST_TRUTH, ST_CONFLICT, ST_MINUS):
        for rsid, gt in table.items():
            v = DB_BY_ID[rsid]
            rows.append((rsid, v.chrom, v.pos, gt))
    rows += [(k, c, p, g) for k, (c, p, g) in ST_ALIAS.items()]
    rows.append(ST_POSONLY)
    rows.append(("rs6152", "X", 66765627, "A" if sex == "male" else "AG"))
    for i in range(120):
        a = "ACGT"[i % 4]
        gt = a if sex == "male" else (a + "ACGT"[(i + 1) % 4] if i % 3 == 0 else a + a)
        rows.append((f"rs9{i:07d}", "X", 3000000 + i * 1000, gt))
    return rows


def _st_write(fmt, path, sex="male"):
    rows = _st_rows(sex)
    if fmt == "23andme":
        head = ["# This data file generated by 23andMe", "# We are using reference human assembly build 37",
                "# rsid\tchromosome\tposition\tgenotype"]
        body = [f"{r}\t{c}\t{p}\t{g}" for r, c, p, g in rows]
        text = "\r\n".join(head + body) + "\r\n"
    elif fmt == "ancestry":
        head = ["#AncestryDNA raw data download", "#coordinates use human reference build 37.1",
                "rsid\tchromosome\tposition\tallele1\tallele2"]
        body = [f"{r}\t{'23' if c == 'X' else c}\t{p}\t{g[0]}\t{g[-1]}" for r, c, p, g in rows]
        text = "\n".join(head + body) + "\n"
    elif fmt == "myheritage":
        head = ["# MyHeritage DNA raw data.", "RSID,CHROMOSOME,POSITION,RESULT"]
        body = [f'"{r}","{c}","{p}","{g}"' for r, c, p, g in rows]
        text = "\n".join(head + body) + "\n"
    elif fmt == "famfinder":
        head = ["# famfinder, https://www.familytreedna.com", "# name,chromosome,position,allele1,allele2"]
        body = [f"{r},{c},{p},{g[0]},{g[-1]}" for r, c, p, g in rows]
        text = "\n".join(head + body) + "\n"
    elif fmt == "vcf":
        head = ["##fileformat=VCFv4.2", "##reference=GRCh37", "##contig=<ID=1,length=249250621>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1"]
        body = []
        for r, c, p, g in rows:
            if r in ST_CONFLICT or r in ST_ALIAS:
                continue
            if r in ST_MINUS:
                g = "".join(sorted(a.translate(COMPLEMENT) for a in g))
            if r == "vg0001":
                body.append(f"1\t{p}\t.\tG\tA\t.\tPASS\t.\tGT\t1/1")
                continue
            if r.startswith("rs9") and len(r) == 10:
                ref, alt = g[0], "ACGT"[("ACGT".index(g[0]) + 1) % 4]
            else:
                ref, alt = DB_BY_ID[r].ref, DB_BY_ID[r].alt
            gt = "/".join(str(int(a == alt)) for a in g) if len(g) == 2 else str(int(g == alt))
            body.append(f"{c}\t{p}\t{r}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}")
        body.append("3\t46414943\trs333\tTCAGTATCAATTCTGGAAGAATTTCCAGACAT\tT\t.\tPASS\t.\tGT\t0/1")
        body.append("11\t46761055\trs1799963\tG\tA\t.\tPASS\t.\tGT\t0/1")
        text = "\n".join(head + body) + "\n"
    else:
        raise ValueError(fmt)
    with open(path, "w", newline="") as fh:
        fh.write(text)
    return path


ST_PGS = """###PGS CATALOG SCORING FILE - selftest
#format_version=2.0
#pgs_id=PGS_SELFTEST
#pgs_name=selftest
#trait_reported=test trait
#genome_build=GRCh37
#weight_type=beta
rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tallelefrequency_effect
rs12913832\t15\t28365618\tG\tA\t0.5\t0.5
rs4988235\t2\t136608646\tA\tG\t1.0\t0.7
rs9999999999\t1\t1\tT\tC\t2.0\t0.1
rs17822931\t16\t48258198\tC\tT\t-1.0\t0.8
"""

ST_LAYOUTS = {
    "23andme_crlf": "# rsid\tchromosome\tposition\tgenotype\r\nrs5\t1\t105\t--\r\nrs6\t1\t106\tGC\r\n"
                    "rs7\t1\t107\tTC\r\nrs8\t1\t108\tAT\r\n",
    "ancestry_multisep": "rsid\tchromosome\tposition\tallele1\tallele2\nrs5\t1\t105\t0\t0\nrs6\t\t1\t106\tG\tC\n"
                         "rs7\t1\t\t 107\tT\tC\nrs8\t1\t108\tA\tT\n",
    "myheritage_quotes": 'RSID,CHROMOSOME,POSITION,RESULT\n"rs5"",""1"",""105"",""--"\n"rs6"",""1"",""106"",""GC"\n'
                         '"rs7"",""1"",""107"",""TC"\n"rs8"",""1"",""108"",""AT"\n',
    "illumina_report": "[Header]\nContent\t\tX\n[Data]\nSample Name\tSNP Name\tChr\tPosition\tAllele1 - Forward\t"
                       "Allele2 - Forward\n1\trs5\t1\t105\t-\t-\n1\trs6\t1\t106\tG\tC\n1\trs7\t1\t107\tT\tC\n"
                       "1\trs8\t1\t108\tA\tT\n",
    "slash_genotypes": "# MARKERNAME\tCHROM\tPOS\tGT\nrs6\tchr1\t106\tG/C\nrs7\tchr1\t107\tT/C\nrs8\tchr1\t108\tA/T\n",
    "no_header": "rs5\t1\t105\t--\nrs6\t1\t106\tGC\nrs7\t1\t107\tTC\nrs8\t1\t108\tAT\n",
}


def selftest() -> int:
    import contextlib
    passed = []

    def check(name, cond, detail=""):
        passed.append(bool(cond))
        print(("PASS  " if cond else "FAIL  ") + name + ("" if cond or not detail else f"  ({detail})"))

    cases = {"AG": ("A", "G"), "ga": ("A", "G"), "--": None, "00": None, "A": ("A",), "DI": ("D", "I"),
             "A/G": ("A", "G"), "N": None, "": None, "AGT": None, "0 0": None, "T|C": ("C", "T")}
    bad = {k: parse_gt(k) for k, want in cases.items() if parse_gt(k) != want}
    check("genotype normalisation", not bad, bad)
    chroms = {"chr23": "X", "26": "MT", "M": "MT", "chrX": "X", "25": "XY", "7": "7"}
    check("chromosome normalisation", all(norm_chrom(k) == v for k, v in chroms.items()))
    vcf = [(("A", "AGC", "GT", "0/1"), ("D", "I")), (("A", ".", "GT", "0/0"), ("A", "A")),
           ((".", "C", "GT", "1/1"), ("C", "C")), (("C", ".", "GT", "./."), None),
           (("G", "T,C", "GT:DP", "1/2:9"), ("C", "T")), (("A", "G", "GT", "1"), ("G",)),
           (("A", "<DEL>", "GT", "0/1"), None)]
    bad = [a for a, want in vcf if _vcf_alleles(*a) != want]
    check("VCF genotype decoding incl. indel D/I and haploid", not bad, bad)

    counts = {"rs12913832": 1, "rs1800407": 0, "rs12896399": 0, "rs16891982": 0, "rs1393350": 1, "rs12203592": 1}
    p = irisplex(counts)
    got = tuple(round(p[k], 3) for k in ("blue", "intermediate", "brown"))
    check("IrisPlex reproduces published webtool output 0.119/0.213/0.668", got == (0.119, 0.213, 0.668), got)

    def fake(pairs):
        cs = CallSet()
        for rsid, gt in pairs.items():
            v = DB_BY_ID[rsid]
            cs[rsid] = Call(v, "ok", tuple(sorted(gt)), gt)
        return cs

    apoe_want = {("TT", "CC"): "e3/e3", ("TT", "CT"): "e2/e3", ("TT", "TT"): "e2/e2", ("CT", "CC"): "e3/e4",
                 ("CT", "CT"): "e2/e4", ("CC", "CC"): "e4/e4"}
    bad = [k for k, want in apoe_want.items()
           if f_apoe(fake({"rs429358": k[0], "rs7412": k[1]}))[0].data.get("genotype") != want]
    check("APOE haplotype table", not bad, bad)
    cyp_want = {("GG", "CC"): "normal", ("AG", "CC"): "intermediate", ("AA", "CC"): "poor",
                ("GG", "CT"): "rapid", ("GG", "TT"): "ultrarapid", ("AG", "CT"): "intermediate"}
    bad = [k for k, want in cyp_want.items()
           if f_cyp2c19(fake({"rs4244285": k[0], "rs4986893": "GG", "rs12248560": k[1]}))[0].data.get(
               "phenotype") != f"{want} metabolizer"]
    check("CYP2C19 phenotype table (CPIC)", not bad, bad)
    dp = f_dpyd(fake({"rs3918290": "CT", "rs55886062": "AA", "rs67376798": "AT", "rs56038477": "CC"}))[0]
    check("DPYD activity score *2A + c.2846A>T = 0.5 (poor)", dp.data.get("activity_score") == 0.5, dp.text)
    tp = f_tpmt(fake({"rs1800462": "CC", "rs1800460": "TT", "rs1142345": "CC"}))[0]
    check("TPMT *3A/*3A poor metabolizer", tp.data.get("tpmt") == "*3A/*3A" and "poor" in tp.text, tp.text)
    ta = f_tas2r38(fake({"rs713598": "GG", "rs1726866": "AG", "rs10246939": "CC"}))[0]
    check("TAS2R38 inconsistent haplotypes flagged", "rare haplotype" in ta.text, ta.text)
    pal = call_variant(DB_BY_ID["rs16891982"], Genome("x", by_id={"rs16891982": Record("5", 1, ("A", "T"), "AT")}))
    check("palindromic C/G SNP is never strand-flipped", pal.status == "conflict", pal.status)

    tmp = tempfile.mkdtemp(prefix="gtd-selftest-")
    try:
        pgs_path = os.path.join(tmp, "score.txt")
        with open(pgs_path, "w") as fh:
            fh.write(ST_PGS)
        pgs = PGS(pgs_path)
        base = _st_write("23andme", os.path.join(tmp, "g23.txt"))
        gz = os.path.join(tmp, "g23.txt.gz")
        with open(base, "rb") as src, gzip.open(gz, "wb") as dst:
            dst.write(src.read())
        anc = _st_write("ancestry", os.path.join(tmp, "anc.txt"))
        zp = os.path.join(tmp, "anc.zip")
        with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.write(anc, "AncestryDNA.txt")
        files = {"23andMe CRLF": base, "23andMe gzip": gz, "AncestryDNA": anc, "AncestryDNA zip": zp,
                 "MyHeritage CSV": _st_write("myheritage", os.path.join(tmp, "mh.csv")),
                 "FTDNA famfinder": _st_write("famfinder", os.path.join(tmp, "ff.csv")),
                 "VCF": _st_write("vcf", os.path.join(tmp, "g.vcf"))}
        ids, pos = wanted([pgs])
        for label, path in files.items():
            g = finalize(load_genome(path, ids, pos), DB)
            cs, fs = analyse(g)
            f = {x.name: x for x in fs if x.cat != "longevity"}
            vcf_like = label == "VCF"
            problems = []
            want = [
                ("build 37", g.build == "37"),
                ("sex male", g.sex.startswith("male")),
                ("APOE e3/e4", f["APOE"].data.get("genotype") == "e3/e4"),
                ("IrisPlex", tuple(round(f["IrisPlex"].data["p"][k], 3)
                                   for k in ("blue", "intermediate", "brown")) == (0.119, 0.213, 0.668)),
                ("MC1R R=2", f["MC1R"].data.get("R") == 2),
                ("TAS2R38 PAV/AVI", f["TAS2R38"].text.startswith("PAV/AVI")),
                ("CYP2C19 *2/*17", f["CYP2C19"].data.get("diplotype") == "*2/*17"),
                ("TPMT *1/*3A", f["TPMT"].data.get("tpmt") == "*1/*3A"),
                ("DPYD AS 1", f["DPYD"].data.get("activity_score") == 1.0),
                ("HFE C282Y/C282Y", f["HFE"].level == 2),
                ("MTHFR 677TT via position", f["MTHFR"].data.get("677T") == 2
                 and cs["rs1801133"].via == "position"),
                ("LCT dosage 1", cs["rs4988235"].dosage == 1 and cs["rs4988235"].flipped != vcf_like),
                ("AR hemizygous A", cs["rs6152"].dosage == 2),
                ("CCR5 delta32 carrier", cs["rs333"].dosage == 1),
                ("F2 via alias", cs["rs1799963"].dosage == 1 and cs["rs1799963"].via in ("i3002432", "rsid")),
                ("ACTN3 conflict", cs["rs1815739"].status == ("absent" if vcf_like else "conflict")),
            ]
            problems = [n for n, ok in want if not ok]
            check(f"{label}: parse and interpret", not problems, ", ".join(problems))
            if label == "23andMe CRLF":
                r = pgs.score(g)
                ok = (abs(r["score"] - 1.5) < 1e-9 and (r["matched"], r["missing"], r["flipped"]) == (3, 1, 1)
                      and abs(r.get("z", 0) - 0.9532) < 1e-3)
                check("PGS scoring: flip, mean imputation, HWE z-score", ok, r)
        plus = {x.name + x.cat: x.text for x in analyse(finalize(load_genome(base, ids, pos), DB))[1]}
        minus = os.path.join(tmp, "minus.txt")
        with open(base) as src, open(minus, "w") as dst:
            for line in src:
                f = line.rstrip("\n").split("\t")
                if not line.startswith("#") and len(f) == 4:
                    f[3] = f[3].translate(COMPLEMENT)
                dst.write("\t".join(f) + "\n")
        g = finalize(load_genome(minus, ids, pos), DB)
        cs, fs = analyse(g)
        diff = [x.name for x in fs if plus.get(x.name + x.cat) != x.text]
        check("minus-strand file: palindromic SNPs complemented, findings unchanged",
              g.strand == "minus" and not diff and cs["rs16891982"].genotype == "G/G", (g.strand, diff))
        fem = _st_write("23andme", os.path.join(tmp, "fem.txt"), sex="female")
        g = finalize(load_genome(fem, ids, pos), DB)
        check("sex inference female", g.sex.startswith("female"), g.sex)
        for name, text in ST_LAYOUTS.items():
            path = os.path.join(tmp, name + ".txt")
            with open(path, "w", newline="") as fh:
                fh.write(text)
            g = load_genome(path, {"rs5", "rs6", "rs7", "rs8"}, set())
            got = {k: g.by_id[k].alleles for k in ("rs6", "rs7", "rs8")}
            nc = g.by_id.get("rs5")
            ok = got == {"rs6": ("C", "G"), "rs7": ("C", "T"), "rs8": ("A", "T")} and (nc is None or nc.alleles is None)
            check(f"layout {name}", ok, got)
        g = finalize(load_genome(base, ids, pos), DB)
        cs, fs = analyse(g)
        items = entries(cs, fs, show_all=True)
        split = {True: {e.title for e in items if e.has is True}, False: {e.title for e in items if e.has is False},
                 None: {e.title for e in items if e.has is None}}
        expect = {True: ("Dry earwax (ABCC11)", "Lactase persistence (MCM6/LCT)", "Red hair genotype (MC1R R/R)",
                         "APOE e4 (Alzheimer's risk allele)", "Reduced CYP2C19 function", "Factor V Leiden (F5)",
                         "Haemochromatosis genotype (HFE C282Y/C282Y)", "PTC/PROP bitter taster (TAS2R38)",
                         "Reduced DPYD activity (fluoropyrimidine toxicity)", "Reduced TPMT activity (thiopurines)"),
                  False: ("Male-pattern baldness risk allele (AR)", "Longevity-favouring APOE (e2, no e4)",
                          "OCA2 R419Q", "HFE H63D"),
                  None: ("Eye colour (IrisPlex)", "Alpha-actinin-3 deficiency (ACTN3 XX)")}
        wrong = [(t, k) for k, titles in expect.items() for t in titles if t not in split[k]]
        check("present/absent/not-determinable classification", not wrong, wrong)
        for truecolor in (False, True):
            pal = Palette(truecolor)
            codes = [pal.take() for _ in range(300)]
            check(f"title palette unique over 300 titles ({'24-bit' if truecolor else '256-colour'})",
                  len(set(codes)) == 300)
        for show_all in (False, True):
            buf = io.StringIO()
            ns = argparse.Namespace(color=True, truecolor=False, all=show_all, category=None, min_evidence="D",
                                    no_sources=True)
            render_text(g, cs, fs, [], ns, buf)
            text = buf.getvalue()
            titles = re.findall(r"\033\[1;(38;[0-9;]+)m(.+?)\033\[0m", text)
            want = len([e for e in entries(cs, fs, show_all=show_all)
                        if show_all or (e.in_file and e.has is not False)])
            present = len([e for e in entries(cs, fs, show_all=show_all) if e.in_file and e.has is True])
            ok = (len(titles) == want and len({c for c, _ in titles}) == len(titles)
                  and f"results   {present} trait present" in text
                  and "FOUND IN YOUR FILE, TRAIT PRESENT" in text
                  and ("TRAIT NOT PRESENT" in text) == show_all and ("NOT IN YOUR FILE" in text) == show_all
                  and "Research use" not in text and "\033[4m" not in text)
            check(f"text report ({'-a' if show_all else 'default'}): sections, bold headings, one colour per title",
                  ok, (len(titles), want))
        for alleles, detected in (("G\tG", False), ("A\tG", True)):
            path = os.path.join(tmp, "cacna1c.txt")
            with open(path, "w") as fh:
                fh.write("rsid\tchromosome\tposition\tallele1\tallele2\n")
                fh.write(f"rs1006737\t12\t2345295\t{alleles}\n")
            gg = finalize(load_genome(path, wanted([])[0], set()), DB)
            cc, ff = analyse(gg)
            outs = {}
            for show_all in (False, True):
                buf = io.StringIO()
                ns = argparse.Namespace(color=False, all=show_all, category=None, min_evidence="D", no_sources=True)
                render_text(gg, cc, ff, [], ns, buf)
                outs[show_all] = buf.getvalue()
            line = f"rs1006737 {alleles.replace(chr(9), '/')} in your file; {int(detected)} " \
                   f"{'copy' if detected else 'copies'} of trait allele A"
            section = "TRAIT PRESENT" if detected else "TRAIT NOT PRESENT"
            ok = (line in outs[True].split(f"FOUND IN YOUR FILE, {section}")[-1]
                  and (line in outs[False]) == detected)
            check(f"rs1006737 {alleles.replace(chr(9), '/')}: listed by default only when A is carried", ok,
                  copies(cc["rs1006737"]))
            buf = io.StringIO()
            lookup(gg, cc, ["rs1006737", "rs0000001"], ns, buf)
            out = buf.getvalue()
            check(f"--snp lookup rs1006737 {alleles.replace(chr(9), '/')}",
                  "position 2345295" in out and ("-> trait present" if detected else "-> trait not present") in out
                  and "rs0000001\n  not in your file" in out, out)
        buf = io.StringIO()
        ns = argparse.Namespace(format="json", pgs=[pgs_path], build=None)
        with contextlib.redirect_stdout(buf):
            run(base, ns, buf)
        doc = json.loads(buf.getvalue())
        check("JSON report round-trip", doc["pgs"][0]["matched"] == 3 and len(doc["variants"]) == len(DB))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"selftest: {sum(passed)}/{len(passed)} passed")
    return 0 if all(passed) else 1


if __name__ == "__main__":
    sys.exit(main())
