"""Utilities for mlst."""
from pathlib import Path
from typing import Union
from pprint import pformat

class SpeciesOptions:
    _species_options = {
        "Achromobacter": "achromobacter",
        "Acinetobacter baumannii": "abaumannii, abaumannii_2",
        "Aeromonas": "aeromonas",
        "Aggregatibacter actinomycetemcomitans": "aactinomycetemcomitans",
        "Anaplasma phagocytophilum": "aphagocytophilum",
        "Arcobacter": "arcobacter",
        "Aspergillus fumigatus": "afumigatus",
        "Bacillus cereus": "bcereus",
        "Bacillus licheniformis": "blicheniformis",
        "Bacillus subtilis": "bsubtilis",
        "Bacteroides fragilis": "bfragilis",
        "Bartonella bacilliformis": "bbacilliformis",
        "Bartonella henselae": "bhenselae",
        "Bartonella washoensis": "bwashoensis",
        "Bordetella": "bordetella",
        "Borrelia": "borrelia",
        "Brachyspira": "brachyspira",
        "Brachyspira hampsonii": "bhampsonii",
        "Brachyspira hyodysenteriae": "bhyodysenteriae",
        "Brachyspira pilosicoli": "bpilosicoli",
        "Brucella intermedia": "bintermedia",
        "Brucella": "brucella",
        "Burkholderia cepacia": "bcepacia",
        "Burkholderia pseudomallei": "bpseudomallei",
        "Campylobacter concisus": "cconcisus",
        "Campylobacter fetus": "cfetus",
        "Campylobacter helveticus": "chelveticus",
        "Campylobacter hyointestinalis": "chyointestinalis",
        "Campylobacter insulaenigrae": "cinsulaenigrae",
        "Campylobacter jejuni": "cjejuni",
        "Campylobacter lanienae": "clanienae",
        "Campylobacter lari": "clari",
        "Campylobacter sputorum": "csputorum",
        "Campylobacter upsaliensis": "cupsaliensis",
        "Candida albicans": "calbicans",
        "Candida glabrata": "cglabrata",
        "Candida krusei": "ckrusei",
        "Candida tropicalis": "ctropicalis",
        "Candidatus liberibacter": "cliberibacter",
        "Carnobacterium maltaromaticum": "cmaltaromaticum",
        "Citrobacter freundii": "cfreundii",
        "Chlamydiales": "chlamydiales",
        "Clostridium botulinum": "cbotulinum",
        "Clostridium difficile": "cdifficile",
        "Clostridium perfringens": "cperfringens",
        "Clostridium septicum": "csepticum",
        "Clonorchis sinensis": "csinensis",
        "Corynebacterium diphtheriae": "cdiphtheriae",
        "Cronobacter": "cronobacter",
        "Cutibacterium acnes": "cacnes",
        "Dichelobacter nodosus": "dnodosus",
        "Edwardsiella": "edwardsiella",
        "Enterobacter cloacae": "ecloacae",
        "Enterococcus faecalis": "efaecalis",
        "Enterococcus faecium": "efaecium",
        "Escherichia coli": "ecoli, ecoli_2",
        "Flavobacterium psychrophilum": "fpsychrophilum",
        "Gallibacterium anatis": "ganatis",
        "Glaesserella parasuis": "gparasuis",
        "Geotrichum": "geotrichum",
        "Helicobacter cinaedi": "hcinaedi",
        "Helicobacter pylori,": "hpylori",
        "Helicobacter suis": "hsuis",
        "Haemophilus influenzae": "hinfluenzae",
        "Kingella kingae": "kkingae",
        "Klebsiella aerogenes": "kaerogenes",
        "Klebsiella oxytoca": "koxytoca", 
        "Klebsiella pneumoniae": "kpneumoniae", 
        "Kudoa septempunctata": "kseptempunctata", 
        "Lactococcus lactis": "llactis",
        "Leptospira": "leptospira, leptospira_2, leptospira_3",
        "Ligilactobacillus salivarius": "lsalivarius",
        "Listeria monocytogenes": "lmonocytogenes",
        "Macrococcus caseolyticus": "mcaseolyticus",
        "Mammaliicoccus sciuri": "msciuri",
        "Mannheimia haemolytica": "mhaemolytica",
        "Melissococcus plutonius": "mplutonius",
        "Microsporum canis": "mcanis",
        "Moraxella catarrhalis": "mcatarrhalis",
        "Mycobacteria": "mycobacteria",
        "Mycobacterium bovis": "mbovis",
        "Mycobacterium massiliense": "mmassiliense",
        "Mycobacteroides abscessus": "mabscessus",
        "Mycoplasma agalactiae": "magalactiae",
        "Mycoplasma anserisalpingitidis": "manserisalpingitidis",
        "Mycoplasma flocculare": "mflocculare",
        "Mycoplasma gallisepticum": "mgallisepticum, mgallisepticum_2",
        "Mycoplasma hominis": "mhominis",
        "Mycoplasma hyopneumoniae": "mhyopneumoniae",
        "Mycoplasma hyorhinis": "mhyorhinis",
        "Mycoplasma iowae": "miowae",
        "Mycoplasma pneumoniae": "mpneumoniae",
        "Mycoplasma synoviae": "msynoviae",
        "Neisseria": "neisseria",
        "Orientia tsutsugamushi": "otsutsugamushi",
        "Ornithobacterium rhinotracheale": "orhinotracheale",
        "Paenibacillus larvae": "plarvae",
        "Pasteurella multocida": "pmultocida, pmultocida_2",
        "Pediococcus pentosaceus": "ppentosaceus",
        "Photobacterium damselae": "pdamselae",
        "Piscirickettsia salmonis": "psalmonis",
        "Porphyromonas gingivalis": "pgingivalis",
        "Propionibacterium acnes": "pacnes",
        "Pseudomonas aeruginosa": "paeruginosa",
        "Pseudomonas fluorescens": "pfluorescens",
        "Pseudomonas putida": "pputida",
        "Rhodococcus": "rhodococcus",
        "Riemerella anatipestifer": "ranatipestifer",
        "Salmonella enterica": "senterica",
        "Shewanella": "shewanella",
        "Sinorhizobium": "sinorhizobium",
        "Staphylococcus aureus": "saureus",
        "Staphylococcus epidermidis": "sepidermidis",
        "Staphylococcus haemolyticus": "shaemolyticus",
        "Staphylococcus hominis": "shominis",
        "Staphylococcus lugdunensis": "slugdunensis",
        "Staphylococcus pseudintermedius": "spseudintermedius",
        "Stenotrophomonas maltophilia": "smaltophilia",
        "Streptococcus agalactiae": "sagalactiae",
        "Streptococcus bovis": "sbovis",
        "Streptococcus canis": "scanis",
        "Staphylococcus chromogenes": "schromogenes",
        "Streptococcus dysgalactiae": "sdysgalactiae",
        "Streptococcus gallolyticus": "sgallolyticus",
        "Streptococcus oralis": "soralis",
        "Streptococcus pneumoniae": "spneumoniae",
        "Streptococcus pyogenes": "spyogenes",
        "Streptococcus suis": "ssuis",
        "Syspastospora parasitica":"sparasitica",
        "Streptococcus thermophilus": "sthermophilus, sthermophilus_2",
        "Streptococcus uberis": "suberis",
        "Streptococcus zooepidemicus": "szooepidemicus",
        "Streptomyces": "streptomyces",
        "Taylorella": "taylorella",
        "Tenacibaculum": "tenacibaculum",
        "Treponema pallidum": "tpallidum",
        "Trichomonas vaginalis": "tvaginalis",
        "Ureaplasma": "ureaplasma",
        "Vibrio": "vibrio",
        "Vibrio cholerae": "vcholerae, vcholerae_2",
        "Vibrio parahaemolyticus": "vparahaemolyticus",
        "Vibrio tapetis": "vtapetis",
        "Vibrio vulnificus": "vvulnificus",
        "Wolbachia": "wolbachia",
        "Xylella fastidiosa": "xfastidiosa",
        "Yersinia pseudotuberculosis": "ypseudotuberculosis",
        "Yersinia ruckeri": "yruckeri"
    }

    def __init__(self, species_options=_species_options):
        self.species_options = species_options

    def is_species_valid(self, species: str) -> bool:
        """Check if provided argument is valid species code."""
        for value in self.species_options.values():
            if species in value:
                return True
        return False

    def print_species_options(self) -> None:
        options = pformat(self.species_options, indent=0)
        chars_to_remove = ["{", "}", "'", ","]
        for char in chars_to_remove:
            options = options.replace(char, "")
        print(options)


class InputMlstyper:
    """Class to store input to run mlst."""
    def __init__(
            self,
            infile: Path,
            species: str,
            database: Path,
            tmp_dir: Path,
            method_path: Union[Path, None],
            outdir_mlstyper: Path = Path.cwd(),
            outdir_mlst_runner: Path = Path.cwd(),
            extented_output: bool = True,
            quiet: bool = True,
            kma_matrix: bool = False,
            save_tmp: bool = False,
            depth: float = 5.0,
            csv_fieldnames: list[str] = ['id', 'sequence_type']
    ):
        self.infile = infile
        self.species = species
        self.database = database
        # The files in tmp_dir are deleted.
        self.tmp_dir = tmp_dir
        # Path to blastn.
        self.method_path = method_path
        # Directory to store the files created by mlstyper. For this script, we
        # need only the results.json file and can be deleted after extracting
        # the sequence type. Therefore, if you don't need the output of
        # mlstyper, provide the path to tmp_dir
        self.outdir_mlstyper = outdir_mlstyper
        self.outdir_mlst_runner = outdir_mlst_runner
        self.extented_output = extented_output
        self.quiet = quiet
        self.kma_matrix = kma_matrix
        self.save_tmp = save_tmp
        self.depth = depth
        # Fieldnames for the `results.csv` file that save the extracted
        # sequence type from the `results.json` created by mlstyper.
        self.csv_fieldnames = csv_fieldnames
