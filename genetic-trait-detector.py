import sys
from colorama import Fore, Style, init

# Initialize colorama for cross-platform color support
init()

class Rsid:
    """Genetic rsid detector class for identifying genetic variants"""

    def __init__(self, rsid: str, alleles: tuple | str, association: str, description: str = "") -> None:
        """Initialize RSID detector with additional description parameter"""
        self.rsid = rsid
        self.alleles = alleles if isinstance(alleles, tuple) else (alleles,)
        self.association = association
        self.description = description or f"Variant associated with {association}"

    def detect(self, line: str) -> tuple[bool, str]:
        """Detect this rsid and return detection status with matched genotype"""
        if self.rsid not in line:
            return False, ""
        
        for allele_format in self.allele_match_set:
            if allele_format in line:
                return True, allele_format
        return False, ""

    def __str__(self) -> str:
        return f"<Rsid[{self.rsid}]({self.association}) - {self.description}>"

    @property
    def allele_match_set(self) -> tuple:
        """Generate SNPedia-compatible allele format variations"""
        def _format_allele(allele: str) -> tuple:
            if ";" in allele:
                a0, a1 = allele.upper().split(";")
                return (f"({a0};{a1})", f"{a0}{a1}", f"{a0}/{a1}", f"{a0}	{a1}")
            return (f"({allele.upper()};{allele.upper()})",)
        return tuple(format for allele in self.alleles for format in _format_allele(allele))

    @property
    def url(self) -> str:
        """Return SNPedia URL for this RSID"""
        return f"https://www.snpedia.com/index.php/{self.rsid}"

    def info(self, matched_allele: str = "") -> str:
        """Generate detailed information string with uniform color"""
        genotype = matched_allele if matched_allele else self.alleles[0]
        color = trait_colors.get(self.association, Fore.WHITE)
        
        # Clean up the genotype display format
        cleaned_genotype = genotype
        if ";" in genotype:
            # Handle formats like "(A;G)" or "A;G"
            if genotype.startswith("(") and genotype.endswith(")"):
                cleaned_genotype = genotype[1:-1]  # Remove parentheses
            parts = cleaned_genotype.split(";")
            cleaned_genotype = f"{parts[0]} {parts[1]}"
        elif "/" in genotype:
            # Handle formats like "A/G"
            parts = genotype.split("/")
            cleaned_genotype = f"{parts[0]} {parts[1]}"
        elif "\t" in genotype:
            # Handle tab-separated format
            parts = genotype.split("\t")
            cleaned_genotype = f"{parts[0]} {parts[1]}"
        
        return (
            f"{color}"
            f"{self.association}: {self.rsid} "
            f"genotype {cleaned_genotype} | "
            f"{self.description} | "
            f"{self.url}"
            f"{Style.RESET_ALL}"
        )

class Trait:
    """Enum-like class for genetic traits"""
    alzheimers = "Alzheimer's Disease"
    autism = "Autism"
    bipolar = "Bipolar Disorder"
    immunity = "Immunity"
    intelligence = "Intelligence"
    longevity = "Longevity"
    metabolism = "Metabolism"
    muscle = "Muscular Performance"
    ocd = "OCD"
    schizophrenia = "Schizophrenia"
    eyes = "Eye Characteristics"
    hair = "Hair Characteristics"
    anxiety = "Anxiety"
    depression = "Depression"
    addiction = "Addiction Risk"

# Unique colors for each trait (visible on black background)
trait_colors = {
    Trait.alzheimers: Fore.YELLOW,          # Yellow
    Trait.autism: Fore.GREEN,              # Green
    Trait.bipolar: Fore.MAGENTA,           # Magenta
    Trait.immunity: Fore.RED,              # Red
    Trait.intelligence: Fore.CYAN,         # Cyan
    Trait.longevity: Fore.BLUE,            # Blue
    Trait.metabolism: Fore.WHITE,          # White
    Trait.muscle: Fore.LIGHTBLACK_EX,      # Gray (visible, replaces black)
    Trait.ocd: Fore.LIGHTMAGENTA_EX,       # Light magenta
    Trait.schizophrenia: Fore.LIGHTYELLOW_EX,  # Light yellow
    Trait.eyes: Fore.LIGHTBLUE_EX,         # Light blue
    Trait.hair: Fore.LIGHTRED_EX,          # Light red
    Trait.anxiety: Fore.LIGHTCYAN_EX,      # Light cyan
    Trait.depression: Fore.LIGHTGREEN_EX,  # Light green
    Trait.addiction: Fore.LIGHTWHITE_EX,   # Light white (replaces dark black)
}

# Updated RSIDs list with corrections and new additions
_RSIDS = [
    # Bipolar Disorder
    ("rs1006737", "A;A", Trait.bipolar, "Calcium channel gene variant linked to bipolar risk"),
    ("rs4027132", "A;A", Trait.bipolar, "Associated with bipolar disorder risk"),
    ("rs7570682", "A;A", Trait.bipolar),
    ("rs1375144", "C;C", Trait.bipolar),
    ("rs683395", "C;T", Trait.bipolar),
    ("rs2609653", "C;T", Trait.bipolar),
    ("rs10982256", "C;C", Trait.bipolar),
    ("rs11622475", "C;C", Trait.bipolar),
    ("rs1344484", "T;T", Trait.bipolar),
    ("rs2953145", "C;G", Trait.bipolar),
    ("rs420259", "T;T", Trait.bipolar),
    ("rs4276227", "C;C", Trait.bipolar),
    ("rs4027132", "A;G", Trait.bipolar, "Heterozygous variant with intermediate bipolar risk"),
    ("rs2609653", "C;C", Trait.bipolar),
    ("rs2953145", "G;G", Trait.bipolar),
    ("rs6448030", "G;G", Trait.bipolar, "Increases bipolar disorder risk"),   
    ("rs10033237", "A;A", Trait.bipolar, "Associated with mood disorders"),   

    # Alzheimer's Disease
    ("rs429358", "C;C", Trait.alzheimers, "APOE4 homozygous variant strongly linked to Alzheimer's risk"),
    ("rs145999145", "A;A", Trait.alzheimers),
    ("rs908832", "A;A", Trait.alzheimers),
    ("rs63750847", "A;A", Trait.alzheimers),
    ("rs429358", "C;T", Trait.alzheimers, "APOE4 heterozygous variant increasing Alzheimer's risk"),
    ("rs145999145", "A;G", Trait.alzheimers),
    ("rs63750847", "A;G", Trait.alzheimers),
    ("rs7412", "C;C", Trait.alzheimers, "APOE ε2/ε2 protective variant"),   
    ("rs4420638", "G;G", Trait.alzheimers, "Near APOC1, increases Alzheimer's risk"),   

    # Autism
    ("rs1858830", "C;C", Trait.autism, "MET gene variant associated with autism risk"),
    ("rs2710102", "C;C", Trait.autism),
    ("rs7794745", "A;T", Trait.autism),
    ("rs1322784", "C;C", Trait.autism),
    ("rs1322784", "C;T", Trait.autism),
    ("rs1322784", "T;T", Trait.autism),
    ("rs265981", "A;G", Trait.autism),
    ("rs4532", "C;T", Trait.autism),
    ("rs686", "A;G", Trait.autism),
    ("rs1143674", "A;A", Trait.autism),
    ("rs6807362", "C;C", Trait.autism),
    ("rs757972971", "A;A", Trait.autism),
    ("rs2217262", "A;A", Trait.autism),
    ("rs6766410", "A;A", Trait.autism),
    ("rs6766410", "A;C", Trait.autism),
    ("rs6766410", "C;C", Trait.autism),
    ("rs1445442", "A;A", Trait.autism),
    ("rs1445442", "A;G", Trait.autism),
    ("rs1445442", "G;G", Trait.autism),
    ("rs2421826", "C;T", Trait.autism),
    ("rs2421826", "T;T", Trait.autism),
    ("rs2421826", "C;C", Trait.autism),
    ("rs1358054", "G;G", Trait.autism),
    ("rs1358054", "G;T", Trait.autism),
    ("rs1358054", "T;T", Trait.autism),
    ("rs536861", "A;A", Trait.autism),
    ("rs536861", "A;C", Trait.autism),
    ("rs536861", "C;C", Trait.autism),
    ("rs722628", "A;A", Trait.autism),
    ("rs722628", "A;G", Trait.autism),
    ("rs722628", "G;G", Trait.autism),
    ("rs1858830", "C;G", Trait.autism),
    ("rs2710102", "G;G", Trait.autism),
    ("rs7794745", "T;T", Trait.autism),
    ("rs265981", "G;G", Trait.autism),
    ("rs4532", "T;T", Trait.autism),
    ("rs686", "A;A", Trait.autism),
    ("rs1143674", "A;A", Trait.autism),
    ("rs757972971", "A;A", Trait.autism),
    ("rs10513025", "T;T", Trait.autism, "Linked to autism spectrum disorder risk"),   
    ("rs1718101", "C;C", Trait.autism, "CNTNAP2 gene, associated with autism susceptibility"),   

    # Schizophrenia
    ("rs27388", "A;A", Trait.schizophrenia),
    ("rs2270641", "G;G", Trait.schizophrenia),
    ("rs4129148", "C;C", Trait.schizophrenia),
    ("rs28694718", "A;A", Trait.schizophrenia),
    ("rs6422441", "C;C", Trait.schizophrenia),
    ("rs28414810", "C;C", Trait.schizophrenia),
    ("rs6603272", "G;G", Trait.schizophrenia),
    ("rs17883192", "C;C", Trait.schizophrenia),
    ("rs165599", "G;G", Trait.schizophrenia),
    ("rs27388", "A;G", Trait.schizophrenia),
    ("rs4129148", "C;G", Trait.schizophrenia),
    ("rs28694718", "A;G", Trait.schizophrenia),
    ("rs6422441", "C;T", Trait.schizophrenia),
    ("rs28414810", "C;G", Trait.schizophrenia),
    ("rs6603272", "G;T", Trait.schizophrenia),
    ("rs17883192", "C;G", Trait.schizophrenia),
    ("rs4633", "C;C", Trait.schizophrenia, "COMT gene, linked to schizophrenia risk"),   
    ("rs13192841", "A;A", Trait.schizophrenia, "Increases schizophrenia susceptibility"),   

    # Longevity
    ("rs3758391", "C;T", Trait.longevity),
    ("rs5882", "A;A", Trait.longevity),
    ("rs1042522", "C;C", Trait.longevity),
    ("rs3803304", "C;C", Trait.longevity),
    ("rs3803304", "C;G", Trait.longevity),
    ("rs3803304", "G;G", Trait.longevity),
    ("rs6873545", "C;C", Trait.longevity),
    ("rs4590183", "C;C", Trait.longevity),
    ("rs1556516", "C;C", Trait.longevity),
    ("rs1556516", "C;G", Trait.longevity),
    ("rs1556516", "G;G", Trait.longevity),
    ("rs7137828", "C;C", Trait.longevity),
    ("rs7137828", "C;T", Trait.longevity),
    ("rs7137828", "T;T", Trait.longevity),
    ("rs1627804", "C;C", Trait.longevity),
    ("rs1627804", "A;A", Trait.longevity),
    ("rs1627804", "A;C", Trait.longevity),
    ("rs7844965", "A;G", Trait.longevity),
    ("rs7844965", "A;A", Trait.longevity),
    ("rs7844965", "G;G", Trait.longevity),
    ("rs61978928", "C;C", Trait.longevity),
    ("rs61978928", "C;T", Trait.longevity),
    ("rs61978928", "T;T", Trait.longevity),
    ("rs28926173", "C;C", Trait.longevity),
    ("rs28926173", "C;T", Trait.longevity),
    ("rs28926173", "T;T", Trait.longevity),
    ("rs146254978", "C;C", Trait.longevity),
    ("rs146254978", "C;T", Trait.longevity),
    ("rs146254978", "T;T", Trait.longevity),
    ("rs139137459", "A;G", Trait.longevity),
    ("rs139137459", "A;A", Trait.longevity),
    ("rs139137459", "G;G", Trait.longevity),
    ("rs3758391", "T;T", Trait.longevity),
    ("rs5882", "A;G", Trait.longevity),
    ("rs1042522", "C;G", Trait.longevity),
    ("rs2542052", "C;C", Trait.longevity, "APOC3 variant, linked to longer lifespan"),   
    ("rs1042714", "G;G", Trait.longevity, "ADRB2, associated with longevity"),   
    ("rs8034191", "C;C", Trait.longevity, "Associated with lung cancer risk and longevity"),
    ("rs2802292", "G;G", Trait.longevity, "FOXO3 variant linked to extended lifespan"),

    # Immunity
    ("rs333", "46373456", Trait.immunity, "CCR5 delta32 variant affecting HIV resistance"),
    ("rs601338", "A;A", Trait.immunity, "FUT2 variant affecting gut microbiome and immunity"),
    ("rs1800896", "T;T", Trait.immunity, "IL10 variant, affects immune response"),   
    ("rs5743810", "C;C", Trait.immunity, "TLR6, linked to infection resistance"),   
    ("rs5743810", "C;T", Trait.immunity, "TLR6, linked to infection resistance"),   
    ("rs5743810", "T;T", Trait.immunity, "TLR6, linked to infection resistance"),   

    # Intelligence
    ("rs28379706", "T;T", Trait.intelligence),
    ("rs28379706", "C;T", Trait.intelligence),
    ("rs28379706", "C;C", Trait.intelligence),
    ("rs363039", "A;G", Trait.intelligence),
    ("rs4680", "A;A", Trait.intelligence, "COMT Val/Val variant, higher COMT activity"),
    ("rs4680", "G;G", Trait.intelligence, "COMT Met/Met variant, better memory function"),
    ("rs363039", "C;C", Trait.intelligence),
    ("rs53576", "G;G", Trait.intelligence, "OXTR variant influencing social cognition"),
    ("rs2251499", "C;C", Trait.intelligence, "Associated with cognitive performance"),   
    ("rs10119", "A;A", Trait.intelligence, "TOMM40, linked to memory function"),   

    # Muscular Performance
    ("rs1815739", "C;C", Trait.muscle, "ACTN3 variant linked to sprint performance"),
    ("rs1805086", "C;C", Trait.muscle),
    ("rs1815739", "C;T", Trait.muscle),
    ("rs1805086", "C;T", Trait.muscle),
    ("rs1815739", "T;T", Trait.muscle),
    ("rs1799752", "D;D", Trait.muscle, "ACE deletion allele linked to endurance"),   
    ("rs4253778", "C;C", Trait.muscle, "PPARGC1A, enhances muscle performance"),   

    # OCD
    ("rs4570625", "G;G", Trait.ocd),
    ("rs4565946", "C;C", Trait.ocd),
    ("rs6313", "A;A", Trait.ocd, "HTR2A, linked to obsessive-compulsive traits"),   
    ("rs301430", "C;C", Trait.ocd, "SLC1A1, associated with OCD risk"),   

    # Metabolism
    ("rs1801131", "A;C", Trait.metabolism),
    ("rs1801131", "C;C", Trait.metabolism),
    ("rs1801133", "C;T", Trait.metabolism),
    ("rs1801133", "T;T", Trait.metabolism),
    ("rs2282679", "A;C", Trait.metabolism),
    ("rs2282679", "C;C", Trait.metabolism),
    ("rs12785878", "G;T", Trait.metabolism),
    ("rs12785878", "T;T", Trait.metabolism),
    ("rs1799945", "F;G", Trait.metabolism),
    ("rs4988235", "C;C", Trait.metabolism),
    ("rs182549", "C;C", Trait.metabolism),
    ("rs2187668", "A;A", Trait.metabolism),
    ("rs2187668", "A;G", Trait.metabolism),
    ("rs5030858", "T;T", Trait.metabolism),
    ("rs72921001", "C;C", Trait.metabolism),
    ("rs7903146", "T;T", Trait.metabolism),
    ("rs7903146", "C;C", Trait.metabolism),
    ("rs7903146", "C;T", Trait.metabolism),
    ("rs662799", "A;G", Trait.metabolism),
    ("rs662799", "G;G", Trait.metabolism),
    ("rs13119723", "A;A", Trait.metabolism),
    ("rs13119723", "A;G", Trait.metabolism),
    ("rs13119723", "G;G", Trait.metabolism),
    ("rs6822844", "G;G", Trait.metabolism),
    ("rs3184504", "C;T", Trait.metabolism),
    ("rs3184504", "T;T", Trait.metabolism),
    ("rs1042713", "G;G", Trait.metabolism, "ADRB2 variant affecting energy expenditure"),
    ("rs9939609", "A;A", Trait.metabolism, "FTO variant linked to obesity risk"),
    ("rs4988235", "T;T", Trait.metabolism, "LCT variant linked to lactose tolerance"),
    ("rs17782313", "C;C", Trait.metabolism, "MC4R, increases obesity risk"),   
    ("rs1421085", "C;C", Trait.metabolism, "FTO, linked to higher BMI"),   

    # Eye Characteristics
    ("rs12913832", "A;A", Trait.eyes, "HERC2 variant strongly linked to blue eyes"),
    ("rs12913832", "A;G", Trait.eyes, "HERC2 variant linked to intermediate eye color"),
    ("rs12913832", "G;G", Trait.eyes, "HERC2 variant linked to brown eyes"),
    ("rs28938473", "T;T", Trait.eyes),
    ("rs61753033", "T;T", Trait.eyes),
    ("rs61753034", "T;T", Trait.eyes),
    ("rs4778241", "A;A", Trait.eyes),
    ("rs4778241", "A;C", Trait.eyes),
    ("rs4778241", "C;C", Trait.eyes),
    ("rs7495174", "A;A", Trait.eyes),
    ("rs1129038", "A;A", Trait.eyes),
    ("rs1129038", "A;G", Trait.eyes),
    ("rs1129038", "G;G", Trait.eyes),
    ("rs916977", "A;A", Trait.eyes),
    ("rs916977", "A;G", Trait.eyes),
    ("rs916977", "G;G", Trait.eyes),
    ("rs1667394", "A;A", Trait.eyes),
    ("rs12203592", "T;T", Trait.eyes, "IRF4 variant affecting eye color variation"),
    ("rs16891982", "C;C", Trait.eyes, "SLC45A2, influences light eye color"),   
    ("rs1393350", "A;A", Trait.eyes, "TYR, linked to lighter eye shades"),   

    # Hair Characteristics
    ("rs6152", "A;A", Trait.hair),
    ("rs6152", "A;G", Trait.hair),
    ("rs6152", "A;", Trait.hair),
    ("rs6152", "G;G", Trait.hair),
    ("rs1805009", "C;C", Trait.hair),
    ("rs1805009", "C;G", Trait.hair),
    ("rs1805007", "C;T", Trait.hair, "MC1R variant linked to red hair"),
    ("rs1805007", "T;T", Trait.hair),
    ("rs1805008", "C;T", Trait.hair),
    ("rs1805008", "T;T", Trait.hair),
    ("rs1805006", "A;A", Trait.hair),
    ("rs1805006", "A;C", Trait.hair),
    ("rs11547464", "A;A", Trait.hair),
    ("rs11547464", "A;G", Trait.hair),
    ("rs35264875", "T;T", Trait.hair),
    ("rs7349332", "T;T", Trait.hair),
    ("rs11803731", "T;T", Trait.hair),
    ("rs17646946", "A;A", Trait.hair),
    ("rs1667394", "A;A", Trait.hair),
    ("rs16891982", "G;G", Trait.hair, "SLC45A2 variant influencing hair pigmentation"),
    ("rs17822931", "C;C", Trait.hair, "ABCC11 variant linked to earwax type and hair texture"),
    ("rs17822931", "C;T", Trait.hair, "Dry earwax, no body odor, likely Asian ancestry"),
    ("rs12821256", "C;C", Trait.hair, "KITLG, associated with blond hair"),   
    ("rs2378249", "G;G", Trait.hair, "PADI3, linked to hair texture"),   

    # Depression
    ("rs6265", "C;C", Trait.depression, "BDNF Val66Met variant linked to depression risk"),
    ("rs25531", "T;T", Trait.depression, "SERT gene variant affecting serotonin transport"),
    ("rs1049353", "A;A", Trait.depression, "CNR1, associated with depression risk"),   
    ("rs909525", "T;T", Trait.depression, "TPH2, linked to serotonin regulation"),   

    # Anxiety
    ("rs1360780", "T;T", Trait.anxiety, "FKBP5 variant associated with stress response"),
    ("rs3810366", "G;G", Trait.anxiety, "CRHR1 variant linked to stress response"),
    ("rs6295", "C;C", Trait.anxiety, "HTR1A, increases anxiety risk"),   
    ("rs6311", "A;A", Trait.anxiety, "HTR2A, linked to stress response"),   

    # Addiction Risk
    ("rs1800497", "A;A", Trait.addiction, "DRD2 Taq1A A1/A1 variant, higher addiction risk"),
    ("rs1800497", "A;G", Trait.addiction, "DRD2 Taq1A, intermediate addiction risk"),  # Added intermediate
    ("rs1799971", "A;G", Trait.addiction, "OPRM1 variant linked to opioid receptor function"),
    ("rs662", "A;A", Trait.addiction, "PON1 variant associated with alcohol metabolism"),
    ("rs1229984", "C;C", Trait.addiction, "ADH1B variant affecting alcohol processing"),
    ("rs16969968", "A;A", Trait.addiction, "CHRNA5 variant associated with nicotine dependence"),
    ("rs2832407", "A;A", Trait.addiction, "GRIK1, increases alcohol dependence risk"),   
    ("rs6277", "G;G", Trait.addiction, "DRD2, linked to reward sensitivity"),   
]

RSIDS = [Rsid(*x) for x in _RSIDS]

def parse_file(path: str) -> dict[str, str]:
    """Parse genetic data file into lookup dictionary"""
    try:
        with open(path, 'r') as genes:
            return {x.split("\t")[0]: x.strip() 
                   for x in genes.readlines() 
                   if x.strip() and x.startswith("rs")}
    except FileNotFoundError:
        print(f"{Fore.RED}Error: File not found at {path}{Style.RESET_ALL}")
        sys.exit(1)
    except Exception as e:
        print(f"{Fore.RED}Error parsing file: {str(e)}{Style.RESET_ALL}")
        sys.exit(1)

def scan_genes(rsid_dict: dict[str, str]) -> list[tuple[Rsid, str]]:
    """Scan genetic data for matching RSIDs and return detected variants with genotypes"""
    detected = []
    for rsid in RSIDS:
        if rsid.rsid in rsid_dict:
            detected_flag, matched_allele = rsid.detect(rsid_dict[rsid.rsid])
            if detected_flag:
                detected.append((rsid, matched_allele))
    return detected

if __name__ == "__main__":
    filename = sys.argv[1] if len(sys.argv) > 1 else input("Enter DNA data file location: ")
    rsid_dict = parse_file(filename)
    detected_traits = scan_genes(rsid_dict)
    
    if not detected_traits:
        print(f"{Fore.YELLOW}No matching genetic variants found in the provided data.{Style.RESET_ALL}")
    else:
        print(f"\n{Fore.WHITE}Found {len(detected_traits)} genetic variants:{Style.RESET_ALL}")
        for rsid, allele in detected_traits:
            print(rsid.info(allele))
            print()
