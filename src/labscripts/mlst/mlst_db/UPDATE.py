#!/usr/bin/env python3
# This scripts updates an already installed MLST database
# The database most have been cloned - git clone https://bitbucket.org/genomicepidemiology/mlst-2.0_db.git
# and it must have been installed - ./INSTALL.py
#
# KMA should be installed before running this script
# The scripts assumes that it is placed in the git repository together with the MLST species directories and the config file
# Usage: ./UPDATE.py [/path/to/kma_index]
# The script can easily be set up to be used for automatic update. For this, use the following syntax:
# ./UPDATE.py [/path/to/kma_index] not_interactive

import shutil, os, sys, subprocess

if len(sys.argv) >= 2 and sys.argv[1] != "not_interactive":
    kma_index = sys.argv[1]
else:
    kma_index = "kma_index"
if sys.argv[-1] != "not_interactive":
    # Check if executable kma_index program is installed, if not promt the user for path
    while shutil.which(kma_index) is None:
        ans = input("Please input path to executable kma_index program or enter 'q'/'quit' to exit:")
        if ans == "q" or ans == "quit":
            print("Exiting!\n\n \
                   Please install executable KMA programs in order to install this database.\n\n \
                   KMA can be obtained from bitbucked:\n\n\
                   git clone https://bitbucket.org/genomicepidemiology/kma.git\n\n\
                   KMA programs must afterwards be compiled:\n\n\
                   gcc -O3 -o kma KMA.c -lm -lpthread\n\
                   gcc -O3 -o kma_index KMA_index.c -lm")
        break
    kma_index = ans
# Check that kma_index program is executable
elif shutil.which(kma_index) == None:
    sys.exit("KMA indexing program, {}, is not executable".format(kma_index))

# Get newest changes from remote repository (git pull)
git_pull_cmd = "git pull"
process = subprocess.Popen(git_pull_cmd, shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
result, err = process.communicate()
result = result.decode('UTF-8')
error = err.decode('UTF-8')

if error != "":
    sys.exit(error)

# Index databases again
# Use config_file to go through database dirs
if "Already up-to-date." not in result:
    config_file = open("config", "r")
    for line in config_file:
        if line.startswith("#"):
            continue
        else:
            line = line.rstrip().split("\t")
            species_dir = line[0] 
            # for each dir index the fasta files
            kma_index_cmd = "{0} -i {1}/{1}.fsa -o {1}/{1}".format(kma_index, species_dir)
            os.system(kma_index_cmd)
    config_file.close() 
else: 
    sys.stdout.write("Already up-to-date.")
