"""
genetic trait detection library
"""
import sys
from colored import fg, attr

class Rsid:
    """genetic rsid detector class"""

    def __init__(
        self,
        rsid: str,
        alleles: tuple | str,
        association: str,
    ) -> None:
        """constructor"""
        self.rsid = rsid
        self.alleles = alleles if isinstance(alleles, tuple) else (alleles,)
        self.association = association

    def detect(self, line: str) -> bool:
        """detect this rsid in the given line of dna data"""
        return self.rsid in line and any(y in line for y in self.allele_match_set)

    def __str__(self) -> str:
        """string representation"""
        return f"<Rsid[{self.rsid}]({self.association})>"

    def __repr__(self) -> str:
        """repl representation"""
        return self.__str__()

    @property
    def single_allele(self) -> bool:
        """if this is a single allele match"""
        return len(self.alleles) == 1

    @property
    def allele_match_set(self) -> tuple:
        """form a set of different formats of alleles"""

        def _set(allele: str) -> tuple:
            """form a set of allele formats"""
            if ";" in allele:
                a0, a1 = allele.upper().split(";")
                return (f"{a0}	{a1}", f"{a0}|{a1}", f"{a0}{a1}")
            return (allele,)

        tuples = ()
        for allele in self.alleles:
            tuples += _set(allele)
        return tuples

    @property
    def url(self) -> str:
        """return the url to snpedia"""
        return f"https://www.snpedia.com/index.php/{self.rsid}"

    @property
    def info(self) -> str:
        """return an informational string about this rsid when detected"""
        allele = (
            f"with the genotype of ({self.alleles[0]}) " if self.single_allele else ""
        )
        return (
            f"You have a trait which is associated with {trait_colors[self.association]}{self.association}{attr(0)} | "
            f"{self.rsid} {allele}{self.url}"
        )
class Trait:
    """traits enum class"""

    alzheimers = "Alzheimer's Disease"
    autism = "Autism"
    bipolar_disorder = "Bipolar Disorder"
    immunity = "Immunity"
    intelligence = "Intelligence"
    longevity = "Longevity"
    metabolism = "Metabolism"
    muscular_performance = "Muscular Performance"
    ocd = "OCD"
    schizophrenia = "Schizophrenia"
    eyes = "Eyes"
    hair = "Hair"
    
trait_colors = {
    Trait.alzheimers: fg('light_yellow'),
    Trait.autism: fg('green'),
    Trait.bipolar_disorder: fg('yellow'),
    Trait.immunity: fg('red'),
    Trait.intelligence: fg('cyan'),
    Trait.longevity: fg('blue'),
    Trait.metabolism: fg('dark_turquoise'),
    Trait.muscular_performance: fg('light_green'),
    Trait.ocd: fg('light_magenta'),
    Trait.schizophrenia: fg('magenta'),
    Trait.eyes: fg('light_blue'),
    Trait.hair: fg('light_red'),
}

_RSIDS = [
    ("rs1006737", "A;A", Trait.bipolar_disorder),
    ("rs4027132", "A;A", Trait.bipolar_disorder),
    ("rs7570682", "A;A", Trait.bipolar_disorder),
    ("rs1375144", "C;C", Trait.bipolar_disorder),
    ("rs683395", "C;T", Trait.bipolar_disorder),
    ("rs2609653", "C;T", Trait.bipolar_disorder),
    ("rs10982256", "C;C", Trait.bipolar_disorder),
    ("rs11622475", "C;C", Trait.bipolar_disorder),
    ("rs1344484", "T;T", Trait.bipolar_disorder),
    ("rs2953145", "C;G", Trait.bipolar_disorder),
    ("rs420259", "T;T", Trait.bipolar_disorder),
    ("rs4276227", "C;C", Trait.bipolar_disorder),
    ("rs4027132", "A;G", Trait.bipolar_disorder),
    ("rs2609653", "C;C", Trait.bipolar_disorder),
    ("rs2953145", "G;G", Trait.bipolar_disorder),
    
    ("rs429358", "C;C", Trait.alzheimers),
    ("rs145999145", "A;A", Trait.alzheimers),
    ("rs908832", "A;A", Trait.alzheimers),
    ("rs63750847", "A;A", Trait.alzheimers),
    ("rs429358", "C;T", Trait.alzheimers),
    ("rs145999145", "A;G", Trait.alzheimers),
    ("rs63750847", "A;G", Trait.alzheimers),
    
    ("rs1858830", "C;C", Trait.autism),
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
    
    ("rs333", ("46373456",), Trait.immunity),
    
    ("rs28379706", "T;T", Trait.intelligence),
    ("rs28379706", "C;T", Trait.intelligence),
    ("rs28379706", "C;C", Trait.intelligence),
    ("rs363039", "A;G", Trait.intelligence),
    ("rs4680", "A;A", Trait.intelligence),
    ("rs363039", "C;C", Trait.intelligence),
    
    ("rs1815739", "C;C", Trait.muscular_performance),
    ("rs1805086", "C;C", Trait.muscular_performance),
    ("rs1815739", "C;T", Trait.muscular_performance),
    ("rs1805086", "C;T", Trait.muscular_performance),
    ("rs1815739", "T;T", Trait.muscular_performance),
    
    ("rs4570625", "G;G", Trait.ocd),
    ("rs4565946", "C;C", Trait.ocd),
    
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
    
    ("rs12913832", "A;A", Trait.eyes),
    ("rs12913832", "A;G", Trait.eyes),
    ("rs12913832", "G;G", Trait.eyes),
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
    
    ("rs6152", "A;A", Trait.hair),
    ("rs6152", "A;G", Trait.hair),
    ("rs6152", "A;", Trait.hair),
    ("rs6152", "G;G", Trait.hair),
    ("rs1805009", "C;C", Trait.hair),
    ("rs1805009", "C;G", Trait.hair),
    ("rs1805007", "C;T", Trait.hair),
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
]
RSIDS = [Rsid(x[0], x[1], x[2]) for x in _RSIDS]


def parse_file(path: str) -> dict[str, str]:
    """parse a file into a dict that we can perform lookups with"""
    with open(path) as genes:
        return {x.split("\t")[0]: x for x in genes.readlines() if x.startswith("rs")}


def scan_genes(rsid_dict: dict[str, str]) -> list[Rsid]:
    """scan the given gene dictionary, returning a list of RSIDS that were detected"""
    return [x for x in RSIDS if x.rsid in rsid_dict and x.detect(rsid_dict[x.rsid])]


if __name__ == "__main__":
    filename = (
        sys.argv[1]
        if len(sys.argv) > 1
        else input("Enter the location of the DNA data file: ")
    )
    rsid_dict = parse_file(filename)
    rsids = scan_genes(rsid_dict)
    for rsid in rsids:
        print(rsid.info + "\n")
