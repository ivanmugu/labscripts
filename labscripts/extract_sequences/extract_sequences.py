"""Make independent sequence files from a file that has multiple sequences.

The target file can have either GenBank or FASTA sequences.
"""
import sys
import argparse
from argparse import Namespace
from pathlib import Path
from typing import Union


class UserInput:
    """Class to store user input."""
    def __init__(
            self,
            infile: Union[Path, None] = None,
            sequence_type: Union[str, None] = None,
            output_folder: Union[Path, None] = None,
    ):
        self.infile = infile
        self.sequence_type = sequence_type
        self.output_folder = output_folder


def parse_command_line_input() -> UserInput:
    """Parse command line input."""
    # Initiate parser
    parser = argparse.ArgumentParser(
        add_help=False,
        prog="extract_sequences",
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Make independent sequence files from a file that has multiple " +
            "sequences."
        )
    )
    # Make arguments groups.
    helper = parser.add_argument_group('Help')
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    # Make help argument.
    helper.add_argument(
        '-h', '--help', action='help', help='Show this help message and exit.'
    )
    # Make required arguments.
    required.add_argument(
        '-i', '--input', required=True,
        help=(
            'Path to input file.\n' +
            'Make sure to provide a file with a correct extension.\n' +
            'GenBank valid extensions: gb, and gbk.\n' +
            'FASTA valid extensions: fasta, fna, ffn, faa, frn, and fa.'
        )
    )
    # Make optional arguments.
    optional.add_argument(
        '-o', '--output',
        help='Path to output folder.\nDefault: current working directory.'
    )
    # Parse command line arguments and check their correctness.
    user_input = parse_command_line_arguments(parser.parse_args())
    # Make sure user provided all required arguments.
    check_if_user_provided_required_arguments(user_input)

    return user_input


def parse_command_line_arguments(command_line_input: Namespace) -> UserInput:
    """Parse command line arguments and check their correctness."""
    # Initate UserInput class.
    user_input = UserInput()
    # Check if user provided infile and if it is valid.
    user_input.infile = check_infile(command_line_input.input)
    # Get input file extension
    user_input.sequence_type = get_input_file_extension(user_input.infile)
    # Check if output folder is valid.
    user_input.output_folder = check_output_folder(command_line_input.output)
    return user_input


def check_infile(infile: Union[str, None]) -> Union[Path, None]:
    """Get infile and check if it's valid."""
    # if user provided input file check if valid.
    if infile:
        infile = Path(infile)
        if not infile.exists():
            sys.exit(f'Error: `{infile}` does not exist')
        if not infile.is_file():
            sys.exit(f'Error: `{infile}` is not a file')
        return infile
    else:
        sys.exit('Error: the `--input` argument is mandatory.')


def get_input_file_extension(infile: Path) -> str:
    """Get input file extension and check if it's valid."""
    extension = infile.name.split('.')[-1:][0]
    if extension == 'gb' or extension == 'bgk':
        return 'gb'
    elif (
        extension == 'fasta' or extension == 'fna' or extension == 'ffn'
        or extension == 'faa' or extension == 'frn' or extension == 'fa'
    ):
        return 'fasta'
    else:
        sys.exit(
            f"Error: submitted input file `{infile.name}` doesn't have a " + 
            "valid extension."
        )


def check_output_folder(output_folder: Union[str, None]) -> Path:
    """Get output folder and check if it's valid."""
    # If user didn't provide output folder, provide current working directory.
    if not output_folder:
        return Path.cwd()
    else:
        output_folder = Path(output_folder)
    # Check if output folder is valid.
    if not output_folder.exists():
        sys.exit(f"Error: {output_folder} does not exist")
    elif not output_folder.is_dir():
        sys.exit(f"Error: {output_folder} does not exist")
    else:
        return output_folder


def extract_gb_files(inputfile: Path, output_folder: Path) -> None:
    """Extract sequences from a GenBank file."""
    with open(inputfile, 'r') as f:
        # Read file by lines until line is empty (EOF).
        while line := f.readline():
            # If line is header open a new file to save the sequence.
            if "LOCUS" in line:
                # Get accession number for naming file.
                name = line.split()[1]
                # Make output path.
                output_file = output_folder / f"{name}.gb"
                # Open new file for writting.
                writter = open(output_file, 'w')
                # Write gb header.
                writter.write(line)
            # If `//` is at the begining of the line close file.
            elif "/" in line[0] and "/" in line[1]:
                writter.close()
            # If line is empty, do nothing.
            elif line == '\n':
                continue
            # Else, concatenate lines into current working file.
            else:
                writter.write(line)


def extract_fasta_files(inputfile: Path, output_folder: Path) -> None:
    """Extract sequences from a FASTA file."""
    with open(inputfile, 'r') as f:
        # counter to keep track of sequences.
        counter = 0
        # Read file by lines.
        while line := f.readline():
            # If line is header open a new file to save the sequence.
            if line[0] == '>':
                # Get accession number for naming file.
                name = get_acc_number_from_fasta_header(line)
                # Make output path.
                output_file = output_folder / f"{name}.fa"
                # If sequence is not the first one, close file before opening a
                # new one.
                if counter > 0:
                    writter.close()
                # Open new file for writting.
                writter = open(output_file, 'w')
                # Write header.
                writter.write(line)
                counter += 1
            # Else, concatenate lines into current working file.
            else:
                writter.write(line)
        # Close last file.
        writter.write(line)


def get_acc_number_from_fasta_header(header: str) -> str:
    """Extract accession number from a fasta sequence header."""
    name = header.split()[0]
    name = name.split('.')[0]
    return name[1:]


def extractor() -> UserInput:
    """Extract GenBank or FASTA sequences from a file."""
    user_input = parse_command_line_input()
    if user_input.sequence_type == 'gb':
        extract_gb_files(
            inputfile=user_input.infile,
            output_folder=user_input.output_folder
        )
    if user_input.sequence_type == 'fasta':
        extract_fasta_files(
            inputfile=user_input.infile,
            output_folder=user_input.output_folder
        )
    return user_input


if __name__ == "__main__":
    user_input = extractor()
    print('Done!')
    print(f'Your new files are in: {user_input.output_folder}')
