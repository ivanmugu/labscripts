import argparse
from importlib import resources
from pathlib import Path
import subprocess
import json
import csv
import sys

from Bio import SeqIO

from labscripts import mlst
from labscripts import run_mlst
from labscripts.mlst.mlst import InputMlstyper, mlstyper

_species_options = [
    "aactinomycetemcomitans", "abaumannii", "abaumannii_2", "achromobacter",
    "aeromonas", "afumigatus", "aphagocytophilum", "arcobacter",
    "bbacilliformis", "bcepacia", "bcereus", "bfragilis", "bhampsonii",
    "bhenselae", "bhyodysenteriae", "bintermedia", "blicheniformis",
    "bordetella", "borrelia", "bpilosicoli", "bpseudomallei", "brachyspira",
    "brucella", "bsubtilis", "bwashoensis", "cacnes", "calbicans",
    "cbotulinum", "cconcisus", "cdifficile", "cdiphtheriae", "cfetus",
    "cfreundii", "cglabrata", "chelveticus", "chlamydiales",
    "chyointestinalis", "cinsulaenigrae", "cjejuni", "ckrusei", "clanienae",
    "clari", "cliberibacter", "cmaltaromaticum", "cperfringens", "cronobacter",
    "csepticum", "csinensis", "csputorum", "ctropicalis", "cupsaliensis",
    "dnodosus", "ecloacae", "ecoli", "ecoli_2", "edwardsiella", "efaecalis",
    "efaecium", "fpsychrophilum", "ganatis", "geotrichum", "gparasuis",
    "hcinaedi", "hinfluenzae", "hpylori", "hsuis", "kaerogenes", "kkingae",
    "koxytoca", "kpneumoniae", "kseptempunctata", "leptospira", "leptospira_2",
    "leptospira_3", "llactis", "lmonocytogenes", "lsalivarius", "mabscessus",
    "magalactiae", "manserisalpingitidis", "mbovis", "mcanis", "mcaseolyticus",
    "mcatarrhalis", "mflocculare", "mgallisepticum", "mgallisepticum_2",
    "mhaemolytica", "mhominis", "mhyopneumoniae", "mhyorhinis", "miowae",
    "mmassiliense", "mplutonius", "mpneumoniae", "msciuri", "msynoviae",
    "mycobacteria", "neisseria", "orhinotracheale", "otsutsugamushi", "pacnes",
    "paeruginosa", "pdamselae", "pfluorescens", "pgingivalis", "plarvae",
    "pmultocida", "pmultocida_2", "ppentosaceus", "pputida", "psalmonis",
    "ranatipestifer", "rhodococcus", "sagalactiae", "saureus", "sbovis",
    "scanis", "schromogenes", "sdysgalactiae", "senterica", "sepidermidis",
    "sgallolyticus", "shaemolyticus", "shewanella", "shominis",
    "sinorhizobium", "slugdunensis", "smaltophilia", "soralis", "sparasitica",
    "spneumoniae", "spseudintermedius", "spyogenes", "ssuis", "sthermophilus",
    "sthermophilus_2", "streptomyces", "suberis", "szooepidemicus",
    "taylorella", "tenacibaculum", "tpallidum", "tvaginalis", "ureaplasma",
    "vcholerae", "vcholerae_2", "vibrio", "vparahaemolyticus", "vtapetis",
    "vvulnificus", "wolbachia", "xfastidiosa", "ypseudotuberculosis",
    "yruckeri"
]

def parse_command_line():
    parser = argparse.ArgumentParser(
        add_help=False,
        prog="run_mlst",
        formatter_class=argparse.RawTextHelpFormatter,
        description="Run mlst with one or more FASTA sequences.",
        epilog=f"List of species options:\n{', '.join(_species_options)}"
    )
    # Make arguments groups.
    helper = parser.add_argument_group("Help")
    required = parser.add_argument_group("Required")
    optional = parser.add_argument_group("Optional")

    # -- ARGUMENTS ------------------------------------------------------------
    # Help argument.
    helper.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit."
    )
    # Required arguments.
    required.add_argument(
        "-i", "--input", nargs="+", required=True,
        help=("Path to file with FASTA sequence(s).")
    )
    required.add_argument(
        "-s", "--species", required=True,
        help=("Species database used from MLST prediction.")
    )
    # Optional arguments.
    optional.add_argument(
        "-o", "--outdir", default=Path.cwd(),
        help=(
            "Output directory.\n" +
            "Default: current working directory."
        )
    )
    optional.add_argument(
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


def check_command_line_arguments(args) -> None:
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
    if args.species not in _species_options:
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
    subprocess.run(f'rm {folder_path}', shell=True)

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

def main():
    # Get path to mlst database
    mlst_package = resources.files(mlst)
    mlst_db = mlst_package / 'mlst_db'
    # Get path to tmp directory
    run_mlst_path = resources.files(run_mlst)
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

if __name__ == "__main__":
    main()