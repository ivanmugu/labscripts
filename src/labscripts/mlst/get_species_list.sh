#!/usr/bin/env bash
# =============================================================================
# 
#      Author: Ivan Munoz Gutierrez
#     Program: get_species_list.sh
#     Created: 2023/06/24
#     Updated: 2023/06/24
# Description: Script to get the optional species names to run mlst from the
#              mlst_db directory.
#
# =============================================================================

# Get path to directory that contains the folder with the alleles for typing
# all the different species.
directory_name=$1

# Iterate over the subdirectories in alleles to get the name of subdirectories
for file in $directory_name/*; do
    # Check if the file is a directory.
    if [ -d "$file" ]; then
        # echo species option.
        echo "${file##*/}"
    fi
done
