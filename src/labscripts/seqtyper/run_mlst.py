import argparse
from importlib import resources
from pathlib import Path
import subprocess
import json
import csv
import sys

from Bio import SeqIO

from labscripts import mlst
from labscripts import seqtyper
from labscripts.mlst.mlst import InputMlstyper, mlstyper

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
    "Vibrio vulnificus ": "vvulnificus",
    "Wolbachia": "wolbachia",
    "Xylella fastidiosa": "xfastidiosa",
    "Yersinia pseudotuberculosis": "ypseudotuberculosis",
    "Yersinia ruckeri": "yruckeri"
}

def parse_command_line():
    parser = argparse.ArgumentParser(
        add_help=False,
        prog="seqtyper",
        formatter_class=argparse.RawTextHelpFormatter,
        description="Run mlst with one or more FASTA sequences.",
        # epilog=f"List of species options:\n{', '.join(_species_options)}"
    )
    helper = parser.add_argument_group("Help")
    helper.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit."
    )
    # -- SUBPARSERS -----------------------------------------------------------
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(
        dest='command', title='Positional arguments', help='Commands',
    )
    # Create subparser for the 'run' seqtyper command
    run = subparsers.add_parser(
        'run', help='Run seqtyper', add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
        )
    subparsers.add_parser('list_sp', help='List species options')

    # -- SUBPARSER run --------------------------------------------------------
    # Make arguments groups.
    run_helper = run.add_argument_group("Help")
    run_required = run.add_argument_group("Required")
    run_optional = run.add_argument_group("Optional")

    # -- ARGUMENTS ------------------------------------------------------------
    # Help argument.
    run_helper.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit."
    )
    # Required arguments.
    run_required.add_argument(
        "-i", "--input", nargs="+", 
        required=True,
        help=("Path to file with FASTA sequence(s).")
    )
    run_required.add_argument(
        "-s", "--species",
        required=True,
        help=("Species database used from MLST prediction.")
    )
    # Optional arguments.
    run_optional.add_argument(
        "-o", "--outdir", default=Path.cwd(),
        help=(
            "Output directory.\n" +
            "Default: current working directory."
        )
    )
    run_optional.add_argument(
        "-mp", "--method_path",
        help=(
            "Path to blastn.\n" +
            "If you don't have blastn in the environmental path, use this\n" +
            "flag to enter the path of blastn."
        )
    )

    args = parser.parse_args()
    check_command_line_arguments(args)
    return args

def is_species_valid(species: str) -> bool:
    for value in _species_options.values():
        if species in value:
            return True
    return False

def check_command_line_arguments(args) -> None:
    if args.command == 'list_sp':
        for key, value in _species_options.items():
            print(f"{key}: {value}.", end=" ")
        print()
        sys.exit(0)
    if not Path(args.input[0]).exists():
        sys.exit(f'Error: {args.input[0]} does not exists')
    if not Path(args.input[0]).is_file():
        sys.exit(f'Error: {args.input[0]} is not a file')
    if not Path(args.outdir).exists():
        sys.exit(f'Error: {args.outdir} does not exists')
    if not Path(args.outdir).is_dir():
        sys.exit(f'Error: {args.outdir} is not a directory')
    if args.method_path and not Path(args.method_path).exists():
        sys.exit(f'Error: {args.method_path} does not exists')
    if args.method_path and not Path(args.method_path).is_file():
        sys.exit(f'Error: {args.method_path} is not a file')
    if not is_species_valid(args.species):
        sys.exit(f'Error: {args.species} is not a valid species option.')

def is_single_fasta_file(infile: Path) -> bool:
    """Check if fasta file has only one fasta sequence."""
    with open(infile, 'r') as f:
        counter = 0
        for line in f:
            if '>' in line:
                counter += 1
            if counter == 2:
                return False
    return True

def extract_sequence_type_from_json(infile: Path) -> str:
    """Get the sequence type from the json file generated by mlst.py"""
    with open(infile, 'r') as f:
        data = json.load(f)
    st = data["mlst"]["results"]["sequence_type"]
    return st

def run_mlstyper_single_fasta(
        input_mlstyper: InputMlstyper, output_path: Path
    ) -> None:
    """Run mlst with single fasta sequence."""
    # Get record id.
    record = SeqIO.read(input_mlstyper.infile[0], 'fasta')
    record_id = record.id
    # Run mlst.
    mlstyper(input_mlstyper)
    # Get st from data.json.
    st = extract_sequence_type_from_json(input_mlstyper.tmp_dir / "data.json")
    # Headers for results.csv
    fieldnames = ['id', 'sequence_type']
    # Make csv file with results.
    with open(output_path / 'results.csv', 'w') as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({'id': record_id, 'sequence_type': st})
    # Delete everything in the tmp folder.
    folder_path = input_mlstyper.tmp_dir / '*'
    subprocess.run(f'rm -rf {folder_path}', shell=True)

def run_mlstyper_multiple_fasta(
        input_mlstyper: InputMlstyper, output_path: Path
    ) -> None:
    """Run mlst with a file with multiple fasta sequences."""
    # Save path to infile.
    path_infile = input_mlstyper.infile[0]
    # Headers for results.csv
    fieldnames = ['id', 'sequence_type']
    # open results.csv to save results.
    output = open(output_path / 'results.csv', 'a')
    # Make a DictWriter object to facilitate saving results.
    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    # Iterate over fasta sequences.
    for sequence in SeqIO.parse(path_infile, 'fasta'):
        # Get record id.
        record_id = sequence.id
        # Path to save the sequence to analyze
        path_sequence = input_mlstyper.tmp_dir / 'sequence.fasta'
        # Write sequence.fasta to analyze
        SeqIO.write(sequence, path_sequence, 'fasta')
        # Change the path to infile in the input_mlstyper class. The path is
        # provided as a list because this is how the mlst script from cge works
        input_mlstyper.infile = [str(path_sequence)]
        # Run mlst.
        mlstyper(input_mlstyper)
        # Get st from data.json
        st = extract_sequence_type_from_json(
            input_mlstyper.tmp_dir / 'data.json'
        )
        writer.writerow({'id': record_id, 'sequence_type': st})
        # Emtpy the tmp folder for the next analysis.
        folder_path = input_mlstyper.tmp_dir / '*'
        subprocess.run(f'rm {folder_path}', shell=True)
    # Close results.csv
    output.close()


def get_sequence_types() -> None:
    # Get path to mlst database
    mlst_package = resources.files(mlst)
    mlst_db = mlst_package / 'mlst_db'
    # Get path to tmp directory
    run_mlst_path = resources.files(seqtyper)
    path_tmp_dir = run_mlst_path / 'tmp'
    # Get user input
    args = parse_command_line()
    # Path to fasta file
    infile = args.input[0]
    # Initialize InputMlstyper
    input_mlstyper = InputMlstyper(
        infile=args.input,
        species=args.species,
        database=mlst_db,
        tmp_dir=path_tmp_dir,
        method_path=args.method_path,
        outdir=path_tmp_dir,
        extented_output=False,
        quiet=True
    )

    if is_single_fasta_file(infile):
        run_mlstyper_single_fasta(input_mlstyper, Path(args.outdir))
    else:
        run_mlstyper_multiple_fasta(input_mlstyper, Path(args.outdir))
    print(f'Done!\nYour results are in: {args.outdir}')

if __name__ == "__main__":
    get_sequence_types()
