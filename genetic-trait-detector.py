#!/usr/bin/env python3
"""Offline trait and variant report for consumer raw genotype data.

Reads 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA, Illumina final
reports, generic rsid/chromosome/position/genotype tables and single-sample VCF,
plain or gzip/bzip2/zip compressed. Every database allele is stored on the GRCh37
plus strand. Non-palindromic calls reported on the opposite strand are complemented;
anything that fits neither orientation is reported as an allele conflict.

  genetic-trait-detector.py genome.zip
  genetic-trait-detector.py genome.txt --format json > report.json
  genetic-trait-detector.py genome.txt --pgs PGS000001.txt.gz
  genetic-trait-detector.py --list
  genetic-trait-detector.py --selftest

Research use only. Chip calls of rare variants are frequently false positives
(Weedon 2021 BMJ 372:n214); confirm medically relevant results with a clinical test.
"""
from __future__ import annotations

import argparse
import bz2
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

__version__ = "2.0.0"

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
EVIDENCE_SHORT = "A clinical grade | B strong, large effect | C replicated, small effect | D functional only"

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

    @property
    def palindromic(self) -> bool:
        return {self.ref, self.alt} in ({"A", "T"}, {"C", "G"})


def _v(rsid, gene, cat, ref, alt, ev, text, src, loc="", hidden=False, levels="000", rare=False, aliases=()):
    chrom, _, pos = loc.partition(":")
    return Variant(rsid, gene, cat, ref, alt, ev, tuple(text), src, chrom, int(pos or 0),
                   hidden, levels, rare, tuple(aliases))


MC1R_SRC = "Valverde 1995 Nat Genet 11:328; Sulem 2007 Nat Genet 39:1443"
IRIS_SRC = "Walsh 2011 Forensic Sci Int Genet 5:170"

DB = [
    _v("rs12913832", "HERC2", "appearance", "A", "G", "B",
       ("brown-eye genotype", "brown or intermediate (green/hazel) eyes most likely", "blue/light eyes expected"),
       "Sturm 2008 Am J Hum Genet 82:424; Eiberg 2008 Hum Genet 123:177", "15:28365618"),
    _v("rs1800407", "OCA2", "appearance", "C", "T", "C",
       ("no R419Q", "R419Q: shifts prediction from brown toward blue/intermediate", "R419Q homozygous"),
       IRIS_SRC, "15:28230318", hidden=True),
    _v("rs12896399", "SLC24A4", "appearance", "G", "T", "C",
       ("no light-pigmentation allele", "one allele associated with light eyes and hair",
        "two alleles associated with light eyes and hair"),
       "Sulem 2007 Nat Genet 39:1443", "14:92773663", hidden=True),
    _v("rs16891982", "SLC45A2", "appearance", "C", "G", "B",
       ("374Leu/Leu: ancestral darker-pigmentation genotype, common outside Europe",
        "374Leu/Phe: one light-pigmentation allele",
        "374Phe/Phe: light-pigmentation genotype typical of Europeans"),
       "Graf 2005 Hum Mutat 25:278; " + IRIS_SRC, "5:33951693"),
    _v("rs1393350", "TYR", "appearance", "G", "A", "C",
       ("no TYR light-eye allele", "blue vs green eyes OR 1.52 per A; more sun sensitivity",
        "blue vs green eyes OR 1.52 per A (two copies); more sun sensitivity"),
       "Sulem 2007 Nat Genet 39:1443", "11:89011046", hidden=True),
    _v("rs12203592", "IRF4", "appearance", "C", "T", "C",
       ("no IRF4 effect allele", "more freckling and sun sensitivity",
        "more freckling and sun sensitivity (two alleles)"),
       "Han 2008 PLoS Genet 4:e1000074; Praetorius 2013 Cell 155:1022", "6:396321"),
    _v("rs1805006", "MC1R", "appearance", "C", "A", "B", ("no D84E", "D84E (R allele)", "D84E homozygous"),
       MC1R_SRC, "16:89985918", hidden=True),
    _v("rs11547464", "MC1R", "appearance", "G", "A", "B", ("no R142H", "R142H (R allele)", "R142H homozygous"),
       MC1R_SRC, "16:89986091", hidden=True),
    _v("rs1805007", "MC1R", "appearance", "C", "T", "B", ("no R151C", "R151C (R allele)", "R151C homozygous"),
       MC1R_SRC, "16:89986117", hidden=True),
    _v("rs1805008", "MC1R", "appearance", "C", "T", "B", ("no R160W", "R160W (R allele)", "R160W homozygous"),
       MC1R_SRC, "16:89986144", hidden=True),
    _v("rs1805009", "MC1R", "appearance", "G", "C", "B", ("no D294H", "D294H (R allele)", "D294H homozygous"),
       MC1R_SRC, "16:89986546", hidden=True, aliases=("i3002507",)),
    _v("rs1805005", "MC1R", "appearance", "G", "T", "C", ("no V60L", "V60L (r allele)", "V60L homozygous"),
       MC1R_SRC, "16:89985844", hidden=True),
    _v("rs2228479", "MC1R", "appearance", "G", "A", "C", ("no V92M", "V92M (r allele)", "V92M homozygous"),
       MC1R_SRC, "16:89985940", hidden=True),
    _v("rs885479", "MC1R", "appearance", "G", "A", "C", ("no R163Q", "R163Q (r allele)", "R163Q homozygous"),
       MC1R_SRC, "16:89986154", hidden=True),
    _v("rs12821256", "KITLG", "appearance", "T", "C", "C",
       ("no KITLG blond allele", "one blond-hair-associated allele", "two blond-hair-associated alleles"),
       "Sulem 2007 Nat Genet 39:1443; Guenther 2014 Nat Genet 46:748", "12:89328335"),
    _v("rs35264875", "TPCN2", "appearance", "A", "T", "C",
       ("no M484L", "M484L: blond-hair-associated allele", "M484L homozygous: two blond-hair-associated alleles"),
       "Sulem 2008 Nat Genet 40:835", "11:68846399"),
    _v("rs3829241", "TPCN2", "appearance", "G", "A", "C",
       ("no G734E", "G734E: blond-hair-associated allele", "G734E homozygous: two blond-hair-associated alleles"),
       "Sulem 2008 Nat Genet 40:835", "11:68855363"),
    _v("rs11803731", "TCHH", "appearance", "T", "A", "C",
       ("no straight-hair allele (curlier hair more likely)", "one straight-hair allele",
        "two straight-hair alleles; most common European genotype (T/T in gene orientation)"),
       "Medland 2009 Am J Hum Genet 85:750", "1:152083325"),
    _v("rs3827760", "EDAR", "appearance", "A", "G", "B",
       ("370Val/Val: ancestral EDAR, typical in Europeans and Africans",
        "370Val/Ala: one hair-thickness allele",
        "370Ala/Ala: thick straight scalp hair; common in East Asians and Native Americans"),
       "Fujimoto 2008 Hum Mol Genet 17:835; Kamberov 2013 Cell 152:691", "2:109513601"),
    _v("rs6152", "AR", "appearance", "G", "A", "C",
       ("common allele; higher odds of male-pattern baldness (G allele OR 2.68)",
        "one lower-risk A allele",
        "lower-risk A allele; A carriers were 24% of non-bald vs 7-13% of bald men"),
       "Hillmer 2005 Am J Hum Genet 77:140; Hayes 2005 Cancer Epidemiol Biomarkers Prev 14:993; "
       "Zhuo 2012 PMID 21981665",
       "X:0"),

    _v("rs17822931", "ABCC11", "senses", "C", "T", "B",
       ("wet earwax", "wet earwax (dry-type carrier)", "dry earwax and reduced axillary odour"),
       "Yoshiura 2006 Nat Genet 38:324; Martin 2010 J Invest Dermatol 130:529", "16:48258198"),
    _v("rs72921001", "OR6A2", "senses", "A", "C", "C",
       ("lowest odds of cilantro tasting soapy (A allele OR 0.81)", "intermediate odds of cilantro tasting soapy",
        "highest odds of cilantro tasting soapy; locus explains ~0.5% of variance"),
       "Eriksson 2012 Flavour 1:22", "11:6889648"),
    _v("rs4481887", "OR2M7", "senses", "G", "A", "C",
       ("highest odds of asparagus-metabolite anosmia",
        "likely able to smell asparagus metabolites (dominant OR 0.48 for anosmia)",
        "most likely able to smell asparagus metabolites"),
       "Eriksson 2010 PLoS Genet 6:e1000993", "1:248496863"),
    _v("rs10427255", "ZEB2", "senses", "T", "C", "C",
       ("baseline odds of photic sneeze reflex", "higher odds of photic sneeze reflex (OR 1.32 per C)",
        "highest odds of photic sneeze reflex"),
       "Eriksson 2010 PLoS Genet 6:e1000993", "2:146125523"),
    _v("rs1815739", "ACTN3", "senses", "C", "T", "B",
       ("RR: alpha-actinin-3 present in fast fibres; enriched in sprint/power athletes", "RX: alpha-actinin-3 present",
        "XX: alpha-actinin-3 deficient (~18% of Europeans); under-represented in elite sprinters"),
       "Yang 2003 Am J Hum Genet 73:627", "11:66328095"),
    _v("rs713598", "TAS2R38", "senses", "C", "G", "B", ("Ala49 (AVI)", "Ala49/Pro49", "Pro49 (PAV)"),
       "Kim 2003 Science 299:1221", "7:141673345", hidden=True),
    _v("rs1726866", "TAS2R38", "senses", "A", "G", "B", ("Val262 (AVI)", "Val262/Ala262", "Ala262 (PAV)"),
       "Kim 2003 Science 299:1221", "7:141672705", hidden=True),
    _v("rs10246939", "TAS2R38", "senses", "T", "C", "B", ("Ile296 (AVI)", "Ile296/Val296", "Val296 (PAV)"),
       "Kim 2003 Science 299:1221", "7:141672604", hidden=True),

    _v("rs4988235", "MCM6/LCT", "diet", "G", "A", "B",
       ("lactase non-persistence likely (adult lactose malabsorption)", "lactase persistence (dominant)",
        "lactase persistence"),
       "Enattah 2002 Nat Genet 30:233; Tishkoff 2007 Nat Genet 39:31 (non-European persistence alleles not assessed)",
       "2:136608646", levels="100"),
    _v("rs762551", "CYP1A2", "diet", "C", "A", "C",
       ("slower caffeine metabolism (*1F absent)", "slower caffeine metabolism",
        "rapid caffeine metabolism (*1F/*1F, highly inducible)"),
       "Sachse 1999 Br J Clin Pharmacol 47:445; Cornelis 2006 JAMA 295:1135", "15:75041917"),
    _v("rs601338", "FUT2", "diet", "G", "A", "B",
       ("secretor", "secretor (non-secretor carrier)",
        "non-secretor (W154X/W154X): largely resistant to symptomatic infection by common norovirus strains"),
       "Lindesmith 2003 Nat Med 9:548; Thorven 2005 J Virol 79:15351", "19:49206674"),
    _v("rs9939609", "FTO", "diet", "T", "A", "C",
       ("no FTO risk allele", "one obesity-risk allele", "~3 kg heavier on average; obesity OR 1.67 vs T/T"),
       "Frayling 2007 Science 316:889", "16:53820527"),
    _v("rs17782313", "MC4R", "diet", "T", "C", "C",
       ("no MC4R risk allele", "one BMI-raising allele", "two BMI-raising alleles"),
       "Loos 2008 Nat Genet 40:768", "18:57851097"),
    _v("rs7903146", "TCF7L2", "diet", "C", "T", "C",
       ("no TCF7L2 risk allele", "type 2 diabetes relative risk 1.45", "type 2 diabetes relative risk 2.41"),
       "Grant 2006 Nat Genet 38:320", "10:114758349", levels="012"),
    _v("rs2282679", "GC", "diet", "T", "G", "C",
       ("no GC low-vitamin-D allele", "modestly lower serum 25(OH)D",
        "lower serum 25(OH)D, higher odds of insufficiency"),
       "Wang 2010 Lancet 376:180", "4:72608383"),
    _v("rs662799", "APOA5", "diet", "A", "G", "C",
       ("no -1131C allele", "higher plasma triglycerides (-1131C)", "higher plasma triglycerides (-1131C/C)"),
       "Pennacchio 2001 Science 294:169", "11:116663707"),
    _v("rs1801133", "MTHFR", "diet", "G", "A", "B", ("677CC", "677CT", "677TT"),
       "Frosst 1995 Nat Genet 10:111", "1:11856378", hidden=True),
    _v("rs1801131", "MTHFR", "diet", "T", "G", "B", ("1298AA", "1298AC", "1298CC"),
       "van der Put 1998 Am J Hum Genet 62:1044", "1:11854476", hidden=True),

    _v("rs671", "ALDH2", "substances", "G", "A", "B",
       ("normal ALDH2 activity", "ALDH2*1/*2: alcohol flushing; drinking raises oesophageal cancer risk",
        "ALDH2*2/*2: near-absent ALDH2 activity, severe flushing"),
       "Brooks 2009 PLoS Med 6:e50", "12:112241766", levels="022"),
    _v("rs1229984", "ADH1B", "substances", "C", "T", "B",
       ("ADH1B*1/*1 (typical in Europeans)", "one fast ADH1B*2 allele: lower risk of alcohol dependence",
        "ADH1B*2/*2: lower risk of alcohol dependence"),
       "Walters 2018 Nat Neurosci 21:1656", "4:100239319"),
    _v("rs16969968", "CHRNA5", "substances", "G", "A", "C",
       ("no D398N risk allele", "heavier smoking if a smoker; higher lung cancer risk",
        "heaviest smoking quantity if a smoker; higher lung cancer risk"),
       "Thorgeirsson 2008 Nature 452:638; Saccone 2010 PLoS Genet 6:e1001053", "15:78882925", levels="011"),

    _v("rs4244285", "CYP2C19", "drugs", "G", "A", "A", ("no *2", "*2 (no function)", "*2/*2"),
       "Lee 2022 Clin Pharmacol Ther 112:959", "10:96541616", hidden=True),
    _v("rs4986893", "CYP2C19", "drugs", "G", "A", "A", ("no *3", "*3 (no function)", "*3/*3"),
       "Lee 2022 Clin Pharmacol Ther 112:959", "10:96540410", hidden=True),
    _v("rs12248560", "CYP2C19", "drugs", "C", "T", "A", ("no *17", "*17 (increased function)", "*17/*17"),
       "Lee 2022 Clin Pharmacol Ther 112:959", "10:96521657", hidden=True),
    _v("rs1799853", "CYP2C9", "drugs", "C", "T", "A", ("no *2", "*2 (decreased function)", "*2/*2"),
       "Johnson 2017 Clin Pharmacol Ther 102:397", "10:96702047", hidden=True),
    _v("rs1057910", "CYP2C9", "drugs", "A", "C", "A", ("no *3", "*3 (no function)", "*3/*3"),
       "Johnson 2017 Clin Pharmacol Ther 102:397", "10:96741053", hidden=True),
    _v("rs9923231", "VKORC1", "drugs", "C", "T", "A", ("-1639GG", "-1639GA", "-1639AA"),
       "Rieder 2005 N Engl J Med 352:2285; Johnson 2017 Clin Pharmacol Ther 102:397", "16:31107689", hidden=True),
    _v("rs4149056", "SLCO1B1", "drugs", "T", "C", "A",
       ("normal SLCO1B1 function", "decreased function: myopathy OR 4.5 on simvastatin 80 mg/day",
        "poor function: myopathy OR 16.9 on simvastatin 80 mg/day"),
       "SEARCH 2008 N Engl J Med 359:789; Cooper-DeHoff 2022 Clin Pharmacol Ther 111:1007", "12:21331549",
       levels="012"),
    _v("rs776746", "CYP3A5", "drugs", "T", "C", "A",
       ("*1/*1: CYP3A5 expressor; CPIC raises the tacrolimus starting dose 1.5-2x",
        "*1/*3: CYP3A5 expressor; CPIC raises the tacrolimus starting dose 1.5-2x",
        "*3/*3: non-expressor (typical in Europeans); standard tacrolimus dosing"),
       "Birdwell 2015 Clin Pharmacol Ther 98:19", "7:99270539", levels="110"),
    _v("rs3892097", "CYP2D6", "drugs", "C", "T", "A",
       ("no *4 (gene deletions/duplications cannot be detected on chips)",
        "one *4 no-function allele: intermediate metabolizer if the other allele is normal",
        "*4/*4: poor metabolizer (codeine, tramadol, tamoxifen, many antidepressants)"),
       "Caudle 2020 Clin Transl Sci 13:116", "22:42524947", levels="012"),
    _v("rs2395029", "HCP5", "drugs", "T", "G", "A",
       ("HLA-B*57:01 tag absent (tag validated in Europeans)",
        "HLA-B*57:01 tag present: abacavir hypersensitivity risk, confirm by HLA typing",
        "HLA-B*57:01 tag present (two copies): abacavir hypersensitivity risk, confirm by HLA typing"),
       "Mallal 2008 N Engl J Med 358:568; Colombo 2008 J Infect Dis 198:864", "6:31431780", levels="022"),
    _v("rs3918290", "DPYD", "drugs", "C", "T", "A", ("no *2A", "*2A (no function)", "*2A/*2A"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", "1:97915614", hidden=True, rare=True),
    _v("rs55886062", "DPYD", "drugs", "A", "C", "A", ("no *13", "*13 (no function)", "*13/*13"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", "1:97981343", hidden=True, rare=True),
    _v("rs67376798", "DPYD", "drugs", "T", "A", "A",
       ("no c.2846A>T", "c.2846A>T (decreased function)", "c.2846A>T homozygous"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", "1:97547947", hidden=True, rare=True),
    _v("rs56038477", "DPYD", "drugs", "C", "T", "A", ("no HapB3", "HapB3 (decreased function)", "HapB3 homozygous"),
       "Amstutz 2018 Clin Pharmacol Ther 103:210", hidden=True),
    _v("rs1800462", "TPMT", "drugs", "C", "G", "A", ("no *2", "*2 (no function)", "*2/*2"),
       "Relling 2019 Clin Pharmacol Ther 105:1095", "6:18143955", hidden=True, rare=True),
    _v("rs1800460", "TPMT", "drugs", "C", "T", "A", ("no 460A", "460A (*3B, or *3A with 719G)", "460A homozygous"),
       "Relling 2019 Clin Pharmacol Ther 105:1095", "6:18139228", hidden=True),
    _v("rs1142345", "TPMT", "drugs", "T", "C", "A", ("no 719G", "719G (*3C, or *3A with 460A)", "719G homozygous"),
       "Relling 2019 Clin Pharmacol Ther 105:1095", "6:18130918", hidden=True),
    _v("rs116855232", "NUDT15", "drugs", "C", "T", "A", ("no *3", "*3 (no function)", "*3/*3"),
       "Relling 2019 Clin Pharmacol Ther 105:1095", hidden=True),

    _v("rs429358", "APOE", "health", "T", "C", "A", ("112Cys/Cys", "112Cys/Arg", "112Arg/Arg"),
       "Corder 1993 Science 261:921", "19:45411941", hidden=True),
    _v("rs7412", "APOE", "health", "C", "T", "A", ("158Arg/Arg", "158Arg/Cys", "158Cys/Cys"),
       "Corder 1993 Science 261:921", "19:45412079", hidden=True),
    _v("rs6025", "F5", "health", "C", "T", "A",
       ("no factor V Leiden", "factor V Leiden heterozygote: venous thrombosis risk ~7-fold",
        "factor V Leiden homozygote: venous thrombosis risk ~80-fold"),
       "Bertina 1994 Nature 369:64; Rosendaal 1995 Blood 85:1504", "1:169519049", levels="022"),
    _v("rs1799963", "F2", "health", "G", "A", "A",
       ("no prothrombin G20210A", "G20210A heterozygote: venous thrombosis risk 2.8-fold",
        "G20210A homozygote: strongly increased venous thrombosis risk"),
       "Poort 1996 Blood 88:3698", "11:46761055", levels="022", aliases=("i3002432",)),
    _v("rs1800562", "HFE", "health", "G", "A", "A", ("no C282Y", "C282Y", "C282Y/C282Y"),
       "Feder 1996 Nat Genet 13:399", "6:26093141", hidden=True),
    _v("rs1799945", "HFE", "health", "C", "G", "A", ("no H63D", "H63D", "H63D/H63D"),
       "Feder 1996 Nat Genet 13:399", "6:26091179", hidden=True),
    _v("rs333", "CCR5", "health", "I", "D", "B",
       ("no CCR5-delta32", "CCR5-delta32 carrier; slower HIV-1 progression reported",
        "no functional CCR5: strong resistance to CCR5-tropic (not CXCR4-tropic) HIV-1"),
       "Samson 1996 Nature 382:722; Liu 1996 Cell 86:367; Dean 1996 Science 273:1856", "3:46414947",
       aliases=("i3003626",)),
    _v("rs63750847", "APP", "health", "C", "T", "B",
       ("no APP A673T", "A673T: protective against Alzheimer's disease and cognitive decline",
        "A673T homozygous: protective"),
       "Jonsson 2012 Nature 488:96", "21:27269932", rare=True),
    _v("rs2187668", "HLA-DQA1", "health", "C", "T", "B", ("no DQ2.5 tag", "DQ2.5 tag", "DQ2.5 tag x2"),
       "Monsuur 2008 PLoS One 3:e2270", "6:32605884", hidden=True),
    _v("rs7454108", "HLA-DQB1", "health", "T", "C", "B", ("no DQ8 tag", "DQ8 tag", "DQ8 tag x2"),
       "Monsuur 2008 PLoS One 3:e2270", "6:32681483", hidden=True),
    _v("rs3184504", "SH2B3", "health", "C", "T", "C",
       ("no R262W risk allele", "autoimmune risk allele (celiac disease, type 1 diabetes)",
        "two autoimmune risk alleles"),
       "Hunt 2008 Nat Genet 40:395", "12:111884608"),

    _v("rs4680", "COMT", "neuro", "G", "A", "D",
       ("Val/Val: higher COMT activity", "Val/Met: intermediate COMT activity",
        "Met/Met: lower COMT activity; behavioural associations inconsistent"),
       "Lachman 1996 Pharmacogenetics 6:243; Chen 2004 Am J Hum Genet 75:807", "22:19951271"),
    _v("rs6265", "BDNF", "neuro", "C", "T", "D",
       ("Val/Val", "Val/Met: reduced activity-dependent BDNF secretion; no depression association in large samples",
        "Met/Met: reduced activity-dependent BDNF secretion; no depression association in large samples"),
       "Egan 2003 Cell 112:257; Border 2019 Am J Psychiatry 176:376", "11:27679916"),
    _v("rs1006737", "CACNA1C", "neuro", "G", "A", "C",
       ("no CACNA1C risk allele", "bipolar disorder OR 1.18 per allele (negligible absolute effect)",
        "bipolar disorder OR 1.18 per allele, two copies (negligible absolute effect)"),
       "Ferreira 2008 Nat Genet 40:1056", "12:2345295"),
    _v("rs10994336", "ANK3", "neuro", "C", "T", "C",
       ("no ANK3 risk allele", "bipolar disorder OR 1.45 per allele (small absolute effect)",
        "bipolar disorder OR 1.45 per allele, two copies (small absolute effect)"),
       "Ferreira 2008 Nat Genet 40:1056", "10:62179812"),

    _v("rs2802292", "FOXO3", "longevity", "T", "G", "C",
       ("no FOXO3 longevity-associated allele", "one longevity-associated allele", "two longevity-associated alleles"),
       "Willcox 2008 Proc Natl Acad Sci USA 105:13987; Flachsbart 2009 Proc Natl Acad Sci USA 106:2700", "6:108908518"),
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
            return "/".join(sorted(self.alleles, key=lambda a: (a != self.v.ref, a)))
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
    text: str
    ev: str
    src: str
    level: int = 0
    data: dict = field(default_factory=dict)


def _na(cat, name, ev, src, rsids):
    return Finding(cat, name, "not assessable: " + ", ".join(rsids) + " absent or no-call", ev, src, -1)


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
    a, b = cs.n("rs429358"), cs.n("rs7412")
    if a is None or b is None:
        return [_na("health", "APOE", "A", APOE_SRC, cs.missing(("rs429358", "rs7412")))]
    geno = APOE_TABLE[(a, b)]
    text, level = APOE_TEXT.get(geno, ("implies the rare e1 haplotype; a genotyping error is more likely, verify", 1))
    out = [Finding("health", "APOE", f"{geno}: {text} (Farrer 1997, Caucasian clinical series)", "A", APOE_SRC,
                   level, {"genotype": geno})]
    if "e1" not in geno:
        if "e4" in geno:
            lt = "e4 carrier: e4 is under-represented among long-lived individuals"
        elif "e2" in geno:
            lt = "e2 carrier without e4: e2 is over-represented among long-lived individuals"
        else:
            lt = "e3/e3: neutral at the strongest longevity locus"
        out.append(Finding("longevity", "APOE", lt, "B", "Deelen 2019 Nat Commun 10:3669", 0, {"genotype": geno}))
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


def irisplex(counts: dict) -> dict:
    lb = IRISPLEX_ALPHA[0] + sum(bb * counts[r] for r, _, bb, _ in IRISPLEX)
    li = IRISPLEX_ALPHA[1] + sum(bi * counts[r] for r, _, _, bi in IRISPLEX)
    eb, ei = math.exp(lb), math.exp(li)
    d = 1.0 + eb + ei
    return {"blue": eb / d, "intermediate": ei / d, "brown": 1.0 / d}


def f_irisplex(cs):
    counts = {r: cs.n(r, a) for r, a, _, _ in IRISPLEX}
    miss = [r for r, k in counts.items() if k is None]
    if miss:
        return [_na("appearance", "IrisPlex eye colour", "B", IRISPLEX_SRC, miss)]
    p = irisplex(counts)
    best = max(p, key=p.get)
    call = best if p[best] >= 0.7 else f"inconclusive below 0.7 (highest: {best})"
    text = f"blue {p['blue']:.3f}, intermediate {p['intermediate']:.3f}, brown {p['brown']:.3f} -> {call}"
    return [Finding("appearance", "IrisPlex eye colour", text, "B", IRISPLEX_SRC, 0, {"p": p, "counts": counts})]


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
    if len(miss) == len(MC1R_R) + len(MC1R_r):
        return [_na("appearance", "MC1R red hair", "B", src, [r for r, _ in MC1R_R + MC1R_r])]
    if big >= 2:
        text, level = "R/R: red hair likely (recessive; MC1R variants almost never share a haplotype, so two R " \
                      "alleles are in trans); fair skin, freckling, higher melanoma risk", 1
    elif big == 1 and small:
        text, level = "R/r: red hair possible; more freckling and sun sensitivity, higher melanoma risk", 1
    elif big == 1:
        text, level = "R/+ carrier: red hair unlikely; more freckling and sun sensitivity, higher melanoma risk", 0
    elif small:
        text, level = "weak r alleles only: small pigmentation effect", 0
    else:
        text, level = "no tested MC1R variant", 0
    if carried:
        text += " [" + ", ".join(carried) + "]"
    return [Finding("appearance", "MC1R red hair", text + _untested(miss), "B", src, level,
                    {"R": big, "r": small, "carried": carried})]


def f_tas2r38(cs):
    ids = ("rs713598", "rs1726866", "rs10246939")
    k = [cs.n(r) for r in ids]
    if None in k:
        return [_na("senses", "TAS2R38 bitter taste", "B", "Kim 2003 Science 299:1221", cs.missing(ids))]
    if k[0] == k[1] == k[2]:
        text = {2: "PAV/PAV: taster, strongest PTC/PROP bitterness", 1: "PAV/AVI: taster",
                0: "AVI/AVI: non-taster of PTC/PROP"}[k[0]]
    else:
        text = "includes a rare haplotype (AAV, AVV, PVI...); taster status uncertain"
    return [Finding("senses", "TAS2R38 bitter taste", text, "B", "Kim 2003 Science 299:1221", 0, {"pav_counts": k})]


def f_hfe(cs):
    src = "Feder 1996 Nat Genet 13:399; Allen 2008 N Engl J Med 358:221"
    c, h = cs.n("rs1800562"), cs.n("rs1799945")
    if c is None:
        return [_na("health", "HFE haemochromatosis", "A", src, ["rs1800562"])]
    h0 = h or 0
    if c == 2:
        text, level = "C282Y/C282Y: haemochromatosis genotype; iron-overload disease in 28.4% of male and 1.2% " \
                      "of female homozygotes; ferritin and transferrin saturation are the follow-up tests", 2
    elif c == 1 and h0 == 1:
        text, level = "C282Y/H63D compound heterozygote: mild iron loading possible, overload disease uncommon", 1
    elif c == 1:
        text, level = "C282Y carrier: no iron overload expected", 0
    elif h0 == 2:
        text, level = "H63D/H63D: clinical iron overload uncommon", 0
    elif h0 == 1:
        text, level = "H63D carrier", 0
    else:
        text, level = "no C282Y or H63D", 0
    if c + h0 > 2:
        text += "; unusual combination, verify genotyping"
    return [Finding("health", "HFE haemochromatosis", text + ("" if h is not None else " (H63D not genotyped)"),
                    "A", src, level, {"C282Y": c, "H63D": h})]


def f_mthfr(cs):
    src = "Frosst 1995 Nat Genet 10:111; van der Put 1998 Am J Hum Genet 62:1044; Hickey 2013 Genet Med 15:153"
    a, b = cs.n("rs1801133"), cs.n("rs1801131")
    if a is None and b is None:
        return [_na("diet", "MTHFR", "B", src, ["rs1801133", "rs1801131"])]
    a0, b0 = a or 0, b or 0
    table = {(2, 0): "677TT: thermolabile enzyme with reduced activity; homocysteine rises mainly when folate is low",
             (1, 1): "677CT/1298AC compound heterozygote: mildly reduced activity",
             (1, 0): "677CT: mildly reduced activity", (0, 2): "1298CC: mildly reduced activity",
             (0, 1): "1298AC: minimal effect", (0, 0): "677CC/1298AA: reference genotype"}
    text = table.get((a0, b0), "unusual combination (677T and 1298C are rarely in cis); verify genotyping")
    text += "; ACMG advises against MTHFR testing for thrombophilia"
    miss = [n for n, k in (("C677T", a), ("A1298C", b)) if k is None]
    return [Finding("diet", "MTHFR", text + _untested(miss), "B", src, 0, {"677T": a, "1298C": b})]


def f_celiac(cs):
    src = "Monsuur 2008 PLoS One 3:e2270"
    dq25, dq8 = cs.n("rs2187668"), cs.n("rs7454108")
    if dq25 is None and dq8 is None:
        return [_na("health", "HLA-DQ celiac", "B", src, ["rs2187668", "rs7454108"])]
    tags = [n for n, k in (("DQ2.5", dq25), ("DQ8", dq8)) if k]
    if tags:
        text = " and ".join(tags) + " tag present: permissive HLA for celiac disease; most carriers never develop it"
    else:
        text = "no DQ2.5/DQ8 tag: celiac disease unlikely (DQ2.2 not assessed)"
    miss = [n for n, k in (("DQ2.5", dq25), ("DQ8", dq8)) if k is None]
    return [Finding("health", "HLA-DQ celiac", text + _untested(miss), "B", src, 0, {"DQ2.5": dq25, "DQ8": dq8})]


def _diplotype(alleles):
    alleles = list(alleles)
    while len(alleles) < 2:
        alleles.insert(0, "*1")
    return "/".join(alleles)


def f_cyp2c19(cs):
    src = "Lee 2022 Clin Pharmacol Ther 112:959 (CPIC)"
    n2, n3, n17 = cs.n("rs4244285"), cs.n("rs4986893"), cs.n("rs12248560")
    if n2 is None or n17 is None:
        return [_na("drugs", "CYP2C19", "A", src, cs.missing(("rs4244285", "rs12248560")))]
    n3 = n3 or 0
    lof, gof = n2 + n3, n17
    if lof + gof > 2:
        return [Finding("drugs", "CYP2C19", "more than two variant alleles; cannot assign a diplotype", "A", src, 1)]
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
    text = f"{dip}: {phen}; relevant to clopidogrel, PPIs, citalopram/escitalopram, voriconazole"
    text += _untested(["*3"] if cs.n("rs4986893") is None else [])
    return [Finding("drugs", "CYP2C19", text, "A", src, level, {"diplotype": dip, "phenotype": phen})]


def f_cyp2c9(cs):
    src = "Johnson 2017 Clin Pharmacol Ther 102:397 (CPIC); Rieder 2005 N Engl J Med 352:2285"
    n2, n3, vk = cs.n("rs1799853"), cs.n("rs1057910"), cs.n("rs9923231")
    parts, level, data = [], 0, {}
    if n2 is None or n3 is None:
        parts.append("CYP2C9 not assessable")
    elif n2 + n3 > 2:
        parts.append("CYP2C9: more than two variant alleles")
        level = 1
    else:
        score = 2 - 0.5 * n2 - 1.0 * n3
        phen = "normal" if score == 2 else "intermediate" if score >= 1 else "poor"
        dip = _diplotype(["*2"] * n2 + ["*3"] * n3)
        parts.append(f"CYP2C9 {dip} ({phen} metabolizer, activity score {score:g})")
        level = max(level, {"normal": 0, "intermediate": 1, "poor": 2}[phen])
        data.update(diplotype=dip, activity_score=score)
    if vk is None:
        parts.append("VKORC1 -1639 not genotyped")
    else:
        parts.append({0: "VKORC1 -1639GG: usual warfarin sensitivity",
                      1: "VKORC1 -1639GA: increased warfarin sensitivity",
                      2: "VKORC1 -1639AA: high warfarin sensitivity"}[vk])
        level = max(level, 1 if vk else 0)
        data["vkorc1_A"] = vk
    if all(p.endswith(("assessable", "genotyped")) for p in parts):
        return [_na("drugs", "CYP2C9 + VKORC1", "A", src, cs.missing(("rs1799853", "rs1057910", "rs9923231")))]
    text = "; ".join(parts) + "; relevant to warfarin, phenytoin, several NSAIDs"
    return [Finding("drugs", "CYP2C9 + VKORC1", text, "A", src, level, data)]


DPYD = (("rs3918290", "*2A", 0.0), ("rs55886062", "*13", 0.0), ("rs67376798", "c.2846A>T", 0.5),
        ("rs56038477", "HapB3", 0.5))


def f_dpyd(cs):
    src = "Amstutz 2018 Clin Pharmacol Ther 103:210 (CPIC)"
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
    if len(miss) == len(DPYD):
        return [_na("drugs", "DPYD", "A", src, [r for r, _, _ in DPYD])]
    if total > 2:
        return [Finding("drugs", "DPYD", "more than two variant alleles; cannot score", "A", src, 2)]
    if score >= 2:
        text, level = "no tested decreased/no-function variant: normal activity assumed", 0
    elif score >= 1:
        text, level = f"activity score {score:g}: intermediate metabolizer; CPIC advises a reduced " \
                      f"fluoropyrimidine (5-FU, capecitabine) starting dose", 2
    else:
        text, level = f"activity score {score:g}: poor metabolizer; CPIC advises avoiding fluoropyrimidines", 2
    if found:
        text += " [" + ", ".join(found) + "; rare variant, confirm clinically]"
    return [Finding("drugs", "DPYD", text + _untested(miss), "A", src, level,
                    {"activity_score": score, "variants": found})]


def f_thiopurine(cs):
    src = "Relling 2019 Clin Pharmacol Ther 105:1095 (CPIC)"
    t2, t3b, t3c, nu = cs.n("rs1800462"), cs.n("rs1800460"), cs.n("rs1142345"), cs.n("rs116855232")
    parts, level, data = [], 0, {}
    if t3b is None or t3c is None:
        parts.append("TPMT not assessable")
    else:
        a3 = min(t3b, t3c)
        nb, nc, n2 = t3b - a3, t3c - a3, t2 or 0
        total = n2 + a3 + nb + nc
        dip = _diplotype(["*2"] * n2 + ["*3A"] * a3 + ["*3B"] * nb + ["*3C"] * nc)
        if total > 2:
            parts.append("TPMT: more than two variant alleles")
            level = 2
        elif total == 0:
            parts.append("TPMT *1/*1: normal metabolizer")
        else:
            phen = "intermediate" if total == 1 else "poor"
            note = "; *3B/*3C in trans (poor) is rare but not excluded" if (a3 == 1 and total == 1) else ""
            parts.append(f"TPMT {dip}: {phen} metabolizer{note}")
            level = 2
        data["tpmt"] = dip
        if t2 is None:
            parts[-1] += " (*2 not genotyped)"
    if nu is None:
        parts.append("NUDT15 not genotyped")
    else:
        parts.append({0: "NUDT15 normal", 1: "NUDT15 *3 carrier: intermediate metabolizer",
                      2: "NUDT15 *3/*3: poor metabolizer"}[nu])
        level = max(level, 2 if nu else 0)
        data["nudt15_T"] = nu
    if all(p.endswith(("assessable", "genotyped")) for p in parts):
        return [_na("drugs", "TPMT + NUDT15", "A", src, cs.missing(("rs1800460", "rs1142345", "rs116855232")))]
    text = "; ".join(parts) + "; relevant to azathioprine, mercaptopurine, thioguanine"
    return [Finding("drugs", "TPMT + NUDT15", text, "A", src, level, data)]


FINDERS = (f_irisplex, f_mc1r, f_tas2r38, f_mthfr, f_cyp2c19, f_cyp2c9, f_dpyd, f_thiopurine,
           f_apoe, f_hfe, f_celiac)


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
DISCLAIMER = ("Research use only. Chip calls of rare variants are often false positives (Weedon 2021 BMJ "
              "372:n214); confirm anything medically relevant with a clinical test before acting on it.")


class Style:
    def __init__(self, on: bool):
        self.on = on

    def __call__(self, codes: str, s: str) -> str:
        return f"\033[{codes}m{s}\033[0m" if self.on and codes else s


def _ev_ok(ev: str, minimum: str) -> bool:
    return EV_ORDER.index(ev) <= EV_ORDER.index(minimum)


def call_text(c: Call):
    v = c.v
    if c.status == "ok":
        d = c.dosage
        text, level = v.text[d], int(v.levels[d])
        if len(c.alleles) == 1:
            text += " (hemizygous)"
        if c.flipped:
            text += " [strand-flipped]"
        if c.unverified:
            text += " [A/T or C/G SNP in a mixed-strand file: orientation unverifiable]"
            level = max(level, 1)
        if c.via == "position":
            text += " [matched by GRCh37 position]"
        elif c.via not in ("", "rsid"):
            text += f" [reported as {c.via}]"
        if v.rare and d:
            text += " [rare variant: confirm clinically]"
            level = max(level, 1)
        return text, level
    if c.status == "conflict":
        return f"genotype {c.raw!r} fits neither strand of {v.ref}/{v.alt}; check file build/format", 1
    return ("no call", -1) if c.status == "nocall" else ("not on this chip", -1)


def _wrap(text: str, width: int, indent: int):
    return textwrap.wrap(text, max(width - indent, 30)) or [""]


def render_text(genome, cs, findings, pgs_results, args, out):
    st = Style(args.color)
    width = shutil.get_terminal_size((110, 24)).columns if args.color else 110
    calls = list(cs.values())
    n = {s: sum(c.status == s for c in calls) for s in ("ok", "nocall", "absent", "conflict")}
    flips = sum(c.flipped for c in calls)
    build = {"36": "NCBI36", "37": "GRCh37", "38": "GRCh38"}.get(genome.build, "build unknown")
    w = out.write
    w(st("1", f"genetic-trait-detector {__version__}") + "\n")
    w(f"file      {genome.path}\n")
    w(f"format    {genome.fmt} | {build}" + (f" ({genome.build_note})" if genome.build_note else "") +
      f" | sex {genome.sex}\n")
    rate = genome.nocalls / genome.markers if genome.markers else 0.0
    w(f"markers   {genome.markers:,} ({genome.nocalls:,} no-calls, {rate:.1%})\n")
    w(f"database  {len(calls)} variants: {n['ok']} called, {n['nocall']} no-call, {n['absent']} not on chip, "
      f"{flips} strand-flipped, {n['conflict']} conflicts\n")
    if genome.strand == "minus":
        w(st("33", "warning   file reports the minus strand; all calls, including A/T and C/G SNPs, were "
                   "complemented") + "\n")
    elif genome.strand == "mixed":
        w(st("33", "warning   file mixes strands; A/T and C/G SNPs cannot be oriented and are flagged") + "\n")
    if genome.build != "37":
        w(f"note      {build}: rsID matching only, no position fallback\n")
    indent = 34
    for cat, (label, color) in CATEGORIES.items():
        if args.category and cat not in args.category:
            continue
        rows = [c for c in calls if c.v.cat == cat and _ev_ok(c.v.ev, args.min_evidence) and
                (c.status == "conflict" or (args.all or (not c.v.hidden and c.status == "ok")))]
        fnds = [f for f in findings if f.cat == cat and _ev_ok(f.ev, args.min_evidence)]
        if not rows and not fnds:
            continue
        w("\n" + st(color + ";1", label) + "\n")
        entries = [(f"  {c.v.rsid:<11} {c.v.gene:<9} {c.genotype:<5}", c.v.ev, *call_text(c), c.v.src) for c in rows]
        entries += [("  " + st("1", f"{f.name:<27}"), f.ev, f.text, f.level, f.src) for f in fnds]
        for head, ev, text, level, src in entries:
            body = _wrap(text, width, indent)
            w(f"{head} [{ev}] " + st(LEVEL_STYLE[level], body[0]) + "\n")
            for extra in body[1:]:
                w(" " * indent + st(LEVEL_STYLE[level], extra) + "\n")
            if not args.no_sources:
                for extra in _wrap(src, width, indent):
                    w(" " * indent + st("2", extra) + "\n")
    for r in pgs_results:
        w("\n" + st("1", f"Polygenic score {r['id']}") + (f" ({r['name']}; {r['trait']})" if r["trait"] else "") + "\n")
        w(f"  {r['variants']:,} variants: {r['matched']:,} matched, {r['missing']:,} missing, {r['flipped']:,} "
          f"strand-flipped, {r['conflicts']:,} conflicts, {r['skipped']:,} skipped (non-SNV or unusable weight)\n")
        w(f"  raw score {r['score']:.6g}\n")
        if "z" in r:
            w(f"  mean-imputed score {r['imputed_score']:.6g} vs expected {r['expected']:.6g} (SD {r['sd']:.4g}): "
              f"z {r['z']:+.2f}, percentile ~{r['percentile']:.0f} (HWE, independent SNPs)\n")
        else:
            w("  no percentile: scoring file lacks effect-allele frequencies for every variant\n")
    w("\nevidence  " + EVIDENCE_SHORT + "\n")
    w(st("2", DISCLAIMER) + "\n")


def _variant_json(c: Call):
    text, level = call_text(c)
    v = c.v
    return {"rsid": v.rsid, "gene": v.gene, "category": v.cat, "ref": v.ref, "effect_allele": v.alt,
            "status": c.status, "genotype": c.genotype if c.status == "ok" else None, "raw": c.raw,
            "dosage": c.dosage, "flipped": c.flipped, "matched_by": c.via or None, "evidence": v.ev,
            "level": level, "interpretation": text, "component": v.hidden, "source": v.src,
            "snpedia": f"https://www.snpedia.com/index.php/{v.rsid}"}


def render_json(genome, cs, findings, pgs_results, args, out):
    doc = {"tool": "genetic-trait-detector", "version": __version__,
           "file": {"path": genome.path, "format": genome.fmt, "build": genome.build or None,
                    "build_note": genome.build_note, "sex": genome.sex, "markers": genome.markers,
                    "nocalls": genome.nocalls},
           "variants": [_variant_json(c) for c in cs.values()],
           "findings": [{"name": f.name, "category": f.cat, "result": f.text, "evidence": f.ev, "level": f.level,
                         "source": f.src, "data": f.data} for f in findings],
           "pgs": pgs_results, "evidence_scale": EVIDENCE, "disclaimer": DISCLAIMER}
    doc["file"]["strand"] = genome.strand
    json.dump(doc, out, indent=2)
    out.write("\n")


def render_tsv(genome, cs, findings, pgs_results, args, out):
    cols = ("kind", "category", "id", "gene", "genotype", "effect_allele", "dosage", "status", "flipped",
            "evidence", "interpretation", "source")
    out.write("\t".join(cols) + "\n")
    for c in cs.values():
        j = _variant_json(c)
        row = ("variant", j["category"], j["rsid"], j["gene"], j["genotype"] or "", j["effect_allele"],
               "" if j["dosage"] is None else j["dosage"], j["status"], int(j["flipped"]), j["evidence"],
               j["interpretation"], j["source"])
        out.write("\t".join(map(str, row)) + "\n")
    for f in findings:
        row = ("finding", f.cat, f.name, "", "", "", "", "ok" if f.level >= 0 else "na", 0, f.ev, f.text, f.src)
        out.write("\t".join(map(str, row)) + "\n")
    for r in pgs_results:
        row = ("pgs", "", r["id"], "", "", "", "", "ok", 0, "", f"score {r['score']:.6g}; "
               f"matched {r['matched']}/{r['variants']}" + (f"; z {r['z']:+.3f}" if "z" in r else ""), r["trait"])
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
            out.write(f"  {v.rsid:<11} {v.gene:<9} {v.ref}>{v.alt}  [{v.ev}] {loc:<13} "
                      f"{'component ' if v.hidden else ''}{st('2', v.src)}\n")
    out.write("\ncompound calls: " + ", ".join(f.__name__[2:] for f in FINDERS) + "\n")
    out.write("alleles: GRCh37 plus strand, ref>effect\n")


def build_parser():
    p = argparse.ArgumentParser(prog="genetic-trait-detector",
                                description="Offline trait and variant report for consumer raw DNA files.",
                                epilog=DISCLAIMER)
    p.add_argument("file", nargs="?", help="raw genotype file (.txt/.csv/.tsv/.vcf, optionally .gz/.bz2/.zip)")
    p.add_argument("-f", "--format", choices=("text", "json", "tsv"), default="text")
    p.add_argument("-c", "--category", action="append", choices=list(CATEGORIES),
                   help="limit output to a category (repeatable)")
    p.add_argument("-e", "--min-evidence", choices=tuple(EV_ORDER), default="D",
                   help="hide entries weaker than this evidence tier (A strongest)")
    p.add_argument("-a", "--all", action="store_true",
                   help="also list component SNPs, no-calls and variants absent from the chip")
    p.add_argument("--pgs", action="append", default=[], metavar="FILE",
                   help="apply a PGS Catalog scoring file (repeatable)")
    p.add_argument("--build", choices=("36", "37", "38"), help="override the detected genome build")
    p.add_argument("--no-color", action="store_true")
    p.add_argument("--no-sources", action="store_true", help="omit literature sources in text output")
    p.add_argument("--list", action="store_true", help="print the variant database and exit")
    p.add_argument("--selftest", action="store_true", help="run built-in tests and exit")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
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
    genome = finalize(load_genome(path, ids, pos), DB, args.build or "")
    cs, findings = analyse(genome)
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
    tp = f_thiopurine(fake({"rs1800462": "CC", "rs1800460": "TT", "rs1142345": "CC", "rs116855232": "CC"}))[0]
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
                ("IrisPlex", tuple(round(f["IrisPlex eye colour"].data["p"][k], 3)
                                   for k in ("blue", "intermediate", "brown")) == (0.119, 0.213, 0.668)),
                ("MC1R R=2", f["MC1R red hair"].data.get("R") == 2),
                ("TAS2R38 PAV/AVI", f["TAS2R38 bitter taste"].text.startswith("PAV/AVI")),
                ("CYP2C19 *2/*17", f["CYP2C19"].data.get("diplotype") == "*2/*17"),
                ("TPMT *1/*3A", f["TPMT + NUDT15"].data.get("tpmt") == "*1/*3A"),
                ("DPYD AS 1", f["DPYD"].data.get("activity_score") == 1.0),
                ("HFE C282Y/C282Y", f["HFE haemochromatosis"].level == 2),
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
