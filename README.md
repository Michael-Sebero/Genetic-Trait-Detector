## How to Install

```
git clone https://github.com/Michael-Sebero/Genetic-Trait-Detector && cd Genetic-Trait-Detector && python3 genetic-trait-detector.py --selftest
```

Requires Python 3.8 or newer and nothing outside the standard library. `colorama` is optional and only helps colour output on older Windows consoles. The self-test runs 36 pass/fail checks and should end with `selftest: 36/36 passed`.

## Description

Genetic-Trait-Detector reads the raw DNA file you download from a consumer genetic testing company and reports which traits your genotypes are associated with. It runs entirely offline: nothing is uploaded and no account or network connection is needed.

The intention of this software is to provide a FOSS alternative for detecting genetic traits. Data breaches occur often and companies shouldn't be trusted with processing or storing genetic data long term. Genetic-Trait-Detector keeps your data on your own machine.

Every variant in the database is stored on the GRCh37 plus strand with its trait allele, a per-genotype interpretation, an evidence tier and a literature citation. Genotypes reported on the opposite strand are corrected automatically and a genotype that fits neither strand is flagged as a conflict rather than guessed.

## Supported files

23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, LivingDNA, Illumina final reports, generic rsid/chromosome/position/genotype tables and single-sample VCF. Files can be plain text or compressed (`.gz`, `.bz2`, `.zip`), so the download can be passed in as-is without unzipping or reformatting.

## This script detects traits associated with

```
Appearance      eye colour (IrisPlex model), red hair (MC1R), blond hair, hair form, skin pigmentation, freckling, male-pattern baldness
Senses          bitter taste (TAS2R38), cilantro soapiness, asparagus smell, photic sneeze, earwax type, ACTN3 muscle fibre type
Diet            lactase persistence, caffeine metabolism, FUT2 secretor status, FTO, MC4R, TCF7L2, vitamin D, triglycerides, MTHFR
Substances      alcohol flush (ALDH2), ADH1B, nicotine dependence (CHRNA5)
Drug response   CYP2C19, CYP2C9, VKORC1, SLCO1B1, CYP3A5, CYP2D6*4, DPYD, TPMT, NUDT15, HLA-B*57:01
Health          APOE, factor V Leiden, prothrombin G20210A, HFE haemochromatosis, CCR5-delta32, APP A673T, celiac HLA-DQ, SH2B3
Neuro           COMT, BDNF, CACNA1C, ANK3
Longevity       FOXO3, APOE
```

`--list` prints the full database of 75 variants with alleles, positions, evidence tiers and sources.

Earlier versions listed autism, schizophrenia, OCD, depression, anxiety and intelligence SNPs. Those were candidate-gene findings that haven't replicated in large studies and no single common SNP predicts these conditions, so they were removed. For polygenic traits, apply a scoring file from the [PGS Catalog](https://www.pgscatalog.org) with `--pgs`.

## How to use Genetic-Trait-Detector

1. Download your raw DNA data from your testing company.
2. Run the script on the downloaded file:

```
python3 genetic-trait-detector.py /path/to/genome.zip
```

Running it without a path prompts for one. The report lists the traits present in your file under `FOUND IN YOUR FILE, TRAIT PRESENT`. The header counts the rest and `-a` also lists traits not present and database SNPs that aren't in your file or were no-calls.

Each result shows the trait, what it means, your genotype as read from the file with the number of copies of the trait allele, the evidence tier and the source:

```
  Blue/light eyes (HERC2)
    Blue/light eyes expected.
    rs12913832 G/G in your file; 2 copies of trait allele G, 2 needed; evidence B (strong)
    Sturm 2008 Am J Hum Genet 82:424; Eiberg 2008 Hum Genet 123:177
```

To check exactly what the script read for any rsID, including ones outside the database:

```
python3 genetic-trait-detector.py /path/to/genome.zip --snp rs1006737
```

Every result title is printed in a colour no other title in the report uses.

```
-f {text,json,tsv}      output format; json and tsv include every assessed variant
-c CATEGORY             limit output to appearance, senses, diet, substances, drugs, health, neuro or longevity (repeatable)
-e {A,B,C,D}            hide results weaker than this evidence tier
-a, --all               also list traits not present, SNPs not in your file and component SNPs
--snp RSID              show exactly what the file holds for an rsID and how it is interpreted (repeatable)
--pgs FILE              apply a PGS Catalog scoring file (repeatable)
--build {36,37,38}      override the detected genome build
--no-color              disable colour (NO_COLOR is also honoured)
--no-sources            omit literature sources
--list                  print the variant database and exit
--selftest              run the built-in tests and exit
```

Evidence tiers:

```
A   clinical grade: guideline-backed or established pathogenic/protective
B   strong association, large effect
C   replicated association, small effect
D   functional variant, trait associations inconsistent
```

### Why isn't a SNP in my file listed?

Having an rsID in your file is not the same as carrying the trait allele. A trait is present only when your genotype contains the allele that defines it, in the number of copies the trait requires. For example, rs1006737 (CACNA1C) is associated with bipolar disorder through its A allele (Ferreira 2008, Nat Genet 40:1056). A file line of `rs1006737  12  2345295  G  G` is read as G/G, which has zero copies of A, so the trait is not present and it isn't listed by default. `-a` shows it under `FOUND IN YOUR FILE, TRAIT NOT PRESENT` and `--snp rs1006737` shows exactly what was read:

```
  Bipolar-associated allele (CACNA1C)
    No CACNA1C risk allele.
    rs1006737 G/G in your file; 0 copies of trait allele A; evidence C (small effect)
```

## DISCLAIMER

This software does not diagnose any disorder or disease and is intended for research and educational use. Each database entry cites the primary literature it is based on. Consumer genotyping chips have a very high false-positive rate for rare variants (Weedon 2021, BMJ 372:n214). Confirm anything medically relevant, including all drug-response and health results, with a clinical-grade test before acting on it.

<p align="center">
  <img src="https://i.postimg.cc/zv6QDz2Q/The-Ooze-MD-09.gif"/>
</p>

## Donations and Contacts

PayPal: https://www.paypal.com/donate/?cmd=_donations&business=YYGU9JWJEE2AG

Email: michaelsebero@disroot.org
