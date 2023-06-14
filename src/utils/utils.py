"""Helper functions."""
import sys
from pathlib import Path
from typing import Union
from argparse import Namespace, _ArgumentGroup


def check_argparse_mandatory_arguments(
        arguments: Namespace, mandatory_group_name: _ArgumentGroup
    ) -> None:
    """Check if user provided the argparse mandatory arguments."""
    # Retrieve a list of the mandatory arguments.
    mandatory_arguments = [
        action.dest for action in mandatory_group_name._group_actions
    ]
    # Iterate over list of aguments to check the value of the attribute.
    for arg in mandatory_arguments:
        if not getattr(arguments, arg):
            sys.exit(f'Error: the `{arg}` argument is mandatory')


def check_infile(infile: Union[str, None]) -> Union[Path, None]:
    """Check if infile is valid and return it as a Path."""
    # If user provided input file, check if valid.
    infile = Path(infile)
    if not infile.exists():
        sys.exit(f'Error: `{infile}` does not exist')
    if not infile.is_file():
        sys.exit(f'Error: `{infile}` is not a file')
    return infile


def check_output_folder(output_folder: Union[str, None]) -> Union[Path, None]:
    """Check path to output folder provided by user and return it as Path.

    If user doesn't provide an output folder, this function returns the output
    foder as the current working directory.
    """
    # If user didn't provide output folder, provide current working directory.
    if not output_folder:
        return Path.cwd()
    else:
        output_folder = Path(output_folder)
    # Check if output folder is valid.
    if not output_folder.exists():
        sys.exit(f"Error: `{output_folder}` does not exist")
    elif not output_folder.is_dir():
        sys.exit(f"Error: `{output_folder}` is not a directory")
    else:
        return output_folder
