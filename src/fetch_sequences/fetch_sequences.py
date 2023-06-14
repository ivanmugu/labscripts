"""Download DNA sequences from the nuccore database."""
import sys
import argparse
from argparse import Namespace
from typing import Union
import re
from pathlib import Path

from Bio import Entrez, SeqIO

from utils.utils import (
    check_argparse_mandatory_arguments, check_infile, check_output_folder
)
from extract_sequences.extract_sequences import (
    extract_fasta_files, extract_gb_files
)


class UserInput:
    """Class to store user input."""
    def __init__(
            self,
            infile: Union[Path, None] = None,
            sequence_type: Union[str, None] = None,
            email: Union[Path, None] = None,
            output_folder: Path = Path.cwd(),
            output_name: Union[str, None] = None,
            path_output_file: Union[Path, None] = None,
            split_sequences: Union[bool, None] = None,
    ):
        self.infile = infile
        self.sequence_type = sequence_type
        self.email = email
        self.output_folder = output_folder
        self.output_name = output_name
        self.path_output_file = path_output_file
        self.split_sequences = split_sequences


def parse_command_line_input() -> UserInput:
    """Parse command line arguments."""
    # Initate parser
    parser = argparse.ArgumentParser(
        add_help=False,
        prog='fetch_sequences',
        formatter_class=argparse.RawTextHelpFormatter,
        description='Download DNA sequences from the nuccore database.'
    )
    # Make arguments groups.
    helper = parser.add_argument_group('Help')
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    # Make help argument.
    helper.add_argument(
        '-h', '--help', action='help',
        help='Show this help message and exit.'
    )
    # Make required arguments.
    required.add_argument(
        '-i', '--input', required=True,
        help='Path to txt file with list of accession numbers.'
    )
    required.add_argument(
        '-t', '--type', required=True, help=(
            'Sequence type to retrieve.\nOptions: `gb` or `fasta`.'
        )
    )
    required.add_argument(
        '-e', '--email', required=True, help='Provide email to NCBI.'
    )
    # Make optional argument.
    optional.add_argument(
        '-o', '--output_folder',
        help='Path to output folder.\nDefault: current working directory.'
    )
    optional.add_argument(
        '-n', '--output_name',
        help=(
            'Name of the output file holding all the retrieved sequences.\n' +
            'Default: depending the sequence type the default name will\n' +
            '`sequences.fasta` or `sequences.gb`.'
        )
    )
    optional.add_argument(
        '-s', '--split_sequences', action="store_true",
        help=(
            'Save sequences as independent files.\n' +
            'If not provided, the sequences will be concatenated in one file.'
        )
    )
    # Parse the command line arguments
    args = parser.parse_args()
    # Make sure user provided all required arguments.
    check_argparse_mandatory_arguments(args, mandatory_group_name=required)
    # Parse command line arguments and check their correctness.
    user_input = parse_command_line_arguments(args)

    return user_input


def parse_command_line_arguments(args: Namespace) -> UserInput:
    """Get command line input and save it in an UserInput class."""
    # Initiate UserInput.
    user_input = UserInput()
    # Check if input is valid.
    user_input.infile = check_infile(args.input)
    # Check if type is valid.
    user_input.sequence_type = check_type(args.type)
    # Get user email.
    user_input.email = check_email(args.email)
    # Check if user provided output folder and if it is valid.
    user_input.output_folder = check_output_folder(
        args.output_folder
    )
    # Check if user provided output_name.
    user_input.output_name = check_output_name(args.output_name)
    # Make path to output file.
    user_input.path_output_file = (
        user_input.output_folder / (
            f'{user_input.output_name}.{user_input.sequence_type}'
        )
    )
    # Check if user wants to create independent files per retrived sequence.
    user_input.split_sequences = args.split_sequences

    return user_input


def check_type(sequence_type: Union[str, None]) -> Union[str, None]:
    """Check if sequence type to retrieve is valid."""
    if sequence_type == 'gb':
        return 'gb'
    elif sequence_type == 'fasta':
        return 'fasta'
    else:
        sys.exit(f'Error: `{sequence_type}` is not a valid type')


def check_email(email: Union[str, None]) -> Union[str, None]:
    """Get email and check if it's valid."""
    # Check if email is valid.
    if not email:
        return None
    if is_valid_email(email):
        return email
    else:
        sys.exit(f'Error: invalid email `{email}`')


def is_valid_email(email: str) -> bool:
    regex = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return re.match(regex, email) is not None


def check_output_name(output_name: str) -> str:
    """Check output name provided by user."""
    # If not output name, output file will be named fetched_sequences.
    if not output_name:
        return 'fetched_sequences'
    else:
        return output_name


def make_acc_number_batches(input_file: Path, batch_size: int = 100) -> list:
    """Make batches of accession numbers."""
    # Open input file to make batches of accession numbers.
    with open(input_file, 'r') as f:
        # Make list of accession numbers.
        acc_list = f.readlines()
        # Make acccession numbers batches.
        batches = []
        for i, item in enumerate(acc_list):
            # Start batch.
            if i % batch_size == 0:
                batch = item.replace("\n", "")
            # End batch and append to batches.
            elif (i % batch_size) == (batch_size - 1):
                batch = batch + "," + item.replace("\n", "")
                batches.append(batch)
            # Add items to batch.
            else:
                batch = batch + "," + item.replace("\n", "")
        # If last batch is smaller than batch_size append it to batches.
        if (batch.count(',') + 1) != batch_size:
            batches.append(batch)

    return batches


def fetcher(batches: list, rettype: str, path_output_file: Path) -> None:
    """Fetch DNA sequences from batches of accession numbers."""
    # Open file to save sequences in append mode.
    output = open(path_output_file, 'a')
    # Variable to keep track of downloaded accession numbers.
    end = 0 
    # Fetch sequences by batches.
    for batch in batches:
        # Start and end of batch.
        start = end + 1
        end = end + batch.count(',') + 1
        print(f"Downloading sequence {start} to {end}")
        # Access nuccore database to retrieve sequences.
        with Entrez.efetch(
            db='nuccore', id=batch, rettype=rettype, retmode='text',
            style='master'
        ) as handle:
            # Parse and save the sequences.
            record = SeqIO.parse(handle, rettype)
            SeqIO.write(record, output, rettype)
    # Close output file.
    output.close()


def extractor(
        sequence_file: Path, sequence_type: str, output_folder: Path
    ) -> None:
    """Extract GenBank or FASTA sequences from a file into independet files."""
    if sequence_type == 'gb':
        extract_gb_files(
            inputfile=sequence_file, output_folder=output_folder
        )
    if sequence_type == 'fasta':
        extract_fasta_files(
            inputfile=sequence_file, output_folder=output_folder
        )


def main():
    # Parse command line arguments.
    user_input = parse_command_line_input()
    # Provide email to NCBI
    Entrez.email = user_input.email
    # Make batches of accession numbers.
    batches = make_acc_number_batches(
        input_file=user_input.infile, batch_size=10
    )
    print('\nConecting to nuccore database to donwload sequences.\n')
    # fetch sequences.
    fetcher(
        batches=batches,
        rettype=user_input.sequence_type,
        path_output_file=user_input.path_output_file
    )
    # If requested, split sequences into individual files
    if user_input.split_sequences:
        extractor(
            sequence_file=user_input.path_output_file,
            sequence_type=user_input.sequence_type,
            output_folder=user_input.output_folder
        )
    # If everything finished correctly, print output messaage
    print('Done!')
    print(
        f'You will find your `{user_input.sequence_type}` file(s) in:\n'+
        f'{user_input.output_folder}'
    )

if __name__ == '__main__':
    main()
