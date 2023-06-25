#!/usr/bin/env python
''' Validate Database for MLST-2.+'''
import sys, os, re

if len(sys.argv) == 2:
   db_path = sys.argv[1]
else:
   db_path = os.path.dirname(os.path.realpath(__file__))

# Check existence of config file
db_config_file = '{}/config'.format(db_path)
if not os.path.exists(db_config_file):
   sys.exit("Error: The database config file could not be found!")

dbs = []
with open(db_config_file) as config_file:
   for line in config_file:
      line = line.rstrip()
      if line == '': continue
      if line.startswith('#'):
         # Check if important files are given in config file
         if 'important files are:' in line.lower():
            important_files = ["{}/{}".format(db_path, db_file.strip())
                               for db_file in line.split('are:')[-1].split(',')]
            # Check all files exist
            for path in important_files:
               if not os.path.exists(path):
                  sys.exit('Error: {} not found!'.format(path))
         continue
      line_data = line.split('\t')
      if len(line_data) != 3:
         sys.exit(("Error: Invalid line in the database config file!\n"
                   "A proper entry requires 3 tab separated columns!\n{}").format(line))
      species_dir = line_data[0].strip()
      species_name = line_data[1].split('#')[0]
      loci_list = line_data[2].split(",")
      # Check if species fasta file is present
      species_path = "{}/{}/{}.fsa".format(db_path, species_dir, species_dir)
      if not os.path.exists(species_path):
         sys.exit(("Error: The database file ({}) could not be found!").format(
                   species_path))
      species_file = open(species_path, "r")
      for line in species_file:
          if line.startswith(">"):
              locus = re.search(">(\w+)[-|_]\d+", line).group(1)
              if locus not in loci_list:
                  sys.exit(("Error: The database file ({}) does not contain the correct locus header\n{} not found in loci list:{}!").format(
                            species_path, locus, ", ".join(loci_list)))                  
      species_file.close()
      # Check if profile file is present
      profile_path = "{}/{}/{}.tsv".format(db_path, species_dir, species_dir)
      if not os.path.exists(profile_path):
         sys.exit(("Error: The database file ({}) could not be found!").format(
                   profile_path))
      dbs.append((species_name, locus))

if len(dbs) == 0:
   sys.exit("Error: No databases were found in the database config file!")
else:
   print("Validation passed. Database is valid.")
