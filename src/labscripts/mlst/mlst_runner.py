"""Module to run mlst from cge.

It can be run using a file with one or more FASTA sequences or a folder
carrying one or more FASTA sequences.
"""
import argparse
from importlib import resources
from pathlib import Path
import subprocess
import json
import csv
import sys
import os

from Bio import SeqIO

from labscripts import mlst
from labscripts.mlst.mlst_utils import SpeciesOptions, InputMlstyper
from labscripts.mlst.mlst_cge import mlstyper


def parse_command_line():
    parser = argparse.ArgumentParser(
        add_help=False,
        prog="mlst",
        formatter_class=argparse.RawTextHelpFormatter,
        description="Run mlst with one or more FASTA sequences.",
        epilog="Note:\nCurrently, the script works only with FASTA sequences."
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
        description="mlst runner.",
        formatter_class=argparse.RawTextHelpFormatter,
        )
    # Create subparser to list species options.
    subparsers.add_parser(
        'list_sp',
        description=(
            "List species options with valid input for `--species` flag."
        ),
        help='List species options'
    )

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
        # Note: if planning to use reads, include `nargas=+` and do the
        # corresponding changes in the run_mlstyper functions.
        "-i", "--input", required=True,
        help=("Path to file or folder with FASTA sequence(s).")
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


def has_fasta_extension(file_name: str) -> bool:
    """Check if file has FASTA extension."""
    extension = file_name.split('.')[-1]
    if (
        extension == 'fasta' or extension == 'fna' or extension == 'ffn' or
        extension == 'faa' or extension == 'frn' or extension == 'fa'
    ):
        return True
    return False

def mk_list_fasta_files(directory: Path) -> list[Path]:
    """Make a list of FASTA files from directory."""
    docs = os.listdir(directory)
    fastas = [directory / doc for doc in docs if has_fasta_extension(doc)]
    return fastas

def list_has_fasta_files(files_list: list[Path]) -> bool:
    """Check if a list of file paths has a FASTA file."""
    for doc in files_list:
        if has_fasta_extension(str(doc)):
            return True
    return False

def check_command_line_arguments(args) -> None:
    species_options = SpeciesOptions()
    if args.command == 'list_sp':
        species_options.print_species_options()
        sys.exit(0)
    if not Path(args.input).exists():
        sys.exit(f'Error: {args.input} does not exist.')
    if Path(args.input).is_file() and not has_fasta_extension(args.input):
        sys.exit(f'Error: {args.input} is not FASTA file.')
    if Path(args.input).is_dir():
        # check if folder has fasta files.
        ls_fasta = mk_list_fasta_files(Path(args.input))
        if not list_has_fasta_files(ls_fasta):
            sys.exit(f'Error: {args.input} does not have any FASTA file.')
    if not Path(args.outdir).exists():
        sys.exit(f'Error: {args.outdir} does not exist.')
    if not Path(args.outdir).is_dir():
        sys.exit(f'Error: {args.outdir} is not a directory.')
    if args.method_path and not Path(args.method_path).exists():
        sys.exit(f'Error: {args.method_path} does not exist.')
    if args.method_path and not Path(args.method_path).is_file():
        sys.exit(f'Error: {args.method_path} is not a file.')
    if not species_options.is_species_valid(args.species):
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
    st = data["mlst_cge"]["results"]["sequence_type"]
    return st

def run_mlstyper_single_fasta(input_mlstyper: InputMlstyper) -> None:
    """Run mlst with single fasta sequence."""
    # Get record id.
    record = SeqIO.read(input_mlstyper.infile, 'fasta')
    record_id = record.id
    # Change the path to infile in the input_mlstyper class. The path is
    # provided as a list because this is how the mlst script from cge works
    input_mlstyper.infile = [str(input_mlstyper.infile)]
    # Run mlst.
    mlstyper(input_mlstyper)
    # Get st from data.json.
    st = extract_sequence_type_from_json(input_mlstyper.tmp_dir / "data.json")
    # Headers for results.csv
    fieldnames = ['id', 'sequence_type']
    # Open and close results file to remove any existing file with same name.
    open(input_mlstyper.outdir_mlst_runner / 'results.csv', 'w').close()
    # Make csv file with results.
    with open(input_mlstyper.outdir_mlst_runner / 'results.csv', 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({'id': record_id, 'sequence_type': st})
    # Delete everything in the tmp folder.
    folder_path = input_mlstyper.tmp_dir / '*'
    subprocess.run(f'rm -rf {folder_path}', shell=True)

def run_mlstyper_multiple_fasta(input_mlstyper: InputMlstyper) -> None:
    """Run mlst with a file with multiple fasta sequences."""
    # Save path to infile.
    path_infile = input_mlstyper.infile
    # Open and close results file to remove any existing file with same name.
    open(input_mlstyper.outdir_mlst_runner / 'results.csv', 'w').close()
    # open results.csv to save results.
    output = open(input_mlstyper.outdir_mlst_runner / 'results.csv', 'a')
    # Make a DictWriter object to facilitate saving results.
    writer = csv.DictWriter(output, fieldnames=input_mlstyper.csv_fieldnames)
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

def run_mlstyper_list_fasta(input_mlstyper: InputMlstyper) -> None:
    """Run mlst with a list of FASTA files.

    FASTA files with more than one sequence are ignored.
    """
    # Save path to directory with FASTA files.
    path_dir = input_mlstyper.infile
    # Open and close results file to remove any existing file with same name.
    open(input_mlstyper.outdir_mlst_runner / 'results.csv', 'w').close()
    # Open file to save results.
    output = open(input_mlstyper.outdir_mlst_runner / 'results.csv', 'a')
    # Make a DictWriter object to facilitate saving results.
    writer = csv.DictWriter(output, fieldnames=input_mlstyper.csv_fieldnames)
    writer.writeheader()
    # Iterate over fasta file paths to perform mlst.
    for fasta in mk_list_fasta_files(path_dir):
        if not is_single_fasta_file(fasta):
            continue
        # Get record id.
        record_id = SeqIO.read(fasta, 'fasta').id
        # Change the path to infile in the input_mlstyper class. The path is
        # provided as a list because this is how the mlst script from cge works
        input_mlstyper.infile = [str(fasta)]
        # Run mlst.
        mlstyper(input_mlstyper)
        # Get st from data.json
        st = extract_sequence_type_from_json(
            input_mlstyper.tmp_dir / 'data.json'
        )
        writer.writerow({'id': record_id, 'sequence_type': st})
        # Empty the tmp folder for the next analysis.
        folder_path = input_mlstyper.tmp_dir / '*'
        subprocess.run(f'rm {folder_path}', shell=True)
    # Close output file.
    output.close()


def get_sequence_types() -> None:
    # Get path to mlst database
    mlst_package = resources.files(mlst)
    mlst_db = mlst_package / 'mlst_db'
    # Get path to tmp directory
    path_tmp_dir = mlst_package / 'tmp'
    # Get user input
    args = parse_command_line()
    # Path to fasta file
    infile = Path(args.input)
    # Initialize InputMlstyper
    input_mlstyper = InputMlstyper(
        infile=Path(args.input),
        species=args.species,
        database=mlst_db,
        tmp_dir=path_tmp_dir,
        method_path=args.method_path,
        outdir_mlstyper=path_tmp_dir,
        outdir_mlst_runner=Path(args.outdir),
        extented_output=False,
        quiet=True
    )

    if input_mlstyper.infile.is_file() and is_single_fasta_file(infile):
        run_mlstyper_single_fasta(input_mlstyper)
    elif input_mlstyper.infile.is_file() and not is_single_fasta_file(infile):
        run_mlstyper_multiple_fasta(input_mlstyper)
    elif input_mlstyper.infile.is_dir():
        run_mlstyper_list_fasta(input_mlstyper)
    print(f'Done!\nYour results are in: {args.outdir}')

if __name__ == "__main__":
    get_sequence_types()
