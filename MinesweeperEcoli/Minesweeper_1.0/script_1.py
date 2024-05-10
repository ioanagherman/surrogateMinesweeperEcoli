#!/usr/bin/env python

# Author: Joshua Rees, joshua.rees@bristol.ac.uk, Life Sciences, University of Bristol
# Modified: Ioana Gherman, ig13470@bristol.ac.uk, Engineering Mathematics, University of Bristol


"""
Stage 1: single gene knockouts are conducted to identify low / no essentiality genes (whose knockout does not prevent cell division) for Stage 2.
Creates single gene knock out simulation bash scripts / exp list / ko list from list of genes.
Expects: Gene file in INPUT_script_1 FOLDER in format gene_name, one per line
Output: bash scripts and / exp list / ko list and gene list text into OUTPUT_script_1
"""

import os
import pandas as pd

# File paths
GENE_INFO_CSV = "gene_name_id_touse.csv"
GENES_TXT = "INPUT_script_1/genes.txt"
USER_INPUT_TXT = "INPUT_script_1/user_input.txt"
KO_LIST = "OUTPUT_script_1/ko{}_ko.list"
EXPERIMENT_SCRIPT = "OUTPUT_script_1/mineinputko{}.sh"
DELETION_LOG = "OUTPUT_final/deletionlog.txt"
#This table contains the name and ID of each gene modelled in the current WCM version
df = pd.read_csv("gene_name_id_touse.csv", header=0, sep=";")
# Functions

def splashscreen():
    """ Fancy splash screen for Minesweeper scripts """
    pass

def user_input():
    """ Ask user for variables to modify template script """

    linecount = 0

    if os.path.isfile(USER_INPUT_TXT):
        linecount = len(open(USER_INPUT_TXT).readlines())

    if linecount < 10:
        input_values = []
        input_values.append(("LAUNCHPAD", "What is the path to your my_lauchpad.yaml file?"))
        input_values.append(("PROJECT_PATH", "What is the path to your wcEcoli folder?"))
        input_values.append(("RUNSCRIPTSFOLDER", "What is the path to your runscripts folder?"))
        input_values.append(("SEEDSNUMBER", "How many seeds per KO would you like to run? (This will apply to all stages)"))

        with open(USER_INPUT_TXT, 'w', newline='\n') as target:
            for key, question in input_values:
                value = input(question)
                target.write(f"{key} {value}\n")

        print("\nYour input values have been saved in INPUT_script_1/user_input.txt - you can update them there if needed.")
        with open(DELETION_LOG, "a+", newline='\n') as log:
            log.write("\nYour input values have been saved in INPUT_script_1/user_input.txt - you can update them there if needed.\n")
    else:
        print("\nYou have pre-existing values in INPUT_script_1/user_input.txt - you can update them there if needed.")

def create_gene_list():
    """ Convert given txt file into correct gene format, create a usable list
        Used by create_scripts() """

    input_genes_txt = GENES_TXT

    with open(input_genes_txt) as gene_file:
        gene_list = [line.strip() for line in gene_file.readlines()]

    return gene_list

def create_scripts(gene_list):
    """ Convert template script to bash script using user input, create exp and ko txt files """

    user_input_values = open(USER_INPUT_TXT, "r+").readlines()[:10]
    value_list = [line.split()[1] for line in user_input_values]

    launchpad_path, project_path, runscripts_path, seeds_number = value_list

    experimentscript = EXPERIMENT_SCRIPT.format("")
    exp_script = open(experimentscript, 'w', newline='\n')


    for i, row in df.iterrows():
        gene = row['Gene']
        rna = row['RNA']
        ko_index = int(row['KO index'])

        if gene in gene_list:
            exp_script.write(f'DESC="KOName {gene} KDRNA {rna}"  \\\n')
            exp_script.write(f'VARIANT="geneKnockout" FIRST_VARIANT_INDEX={ko_index} LAST_VARIANT_INDEX={ko_index} \\\n')
            exp_script.write(f'CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \\\n')
            exp_script.write(f'SINGLE_DAUGHTERS=0 N_GENS=1 N_INIT_SIMS={seeds_number} \\\n')
            exp_script.write(f'MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 \\\n')
            exp_script.write(f'LAUNCHPAD_FILE="{launchpad_path}/my_launchpad.yaml" \\\n')
            exp_script.write(f'python {runscripts_path}/fireworks/fw_queue.py\n\n')

    exp_script.write('nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 40 --maxjobs_queue 100 \n')

    print(f'\nCreated bash script (*.sh), in OUTPUT_script_1.')
    print(f'\n> Now you should copy the bash script to your wcEcoli/runscripts folder.')

    with open(DELETION_LOG, "a+", newline='\n') as log:
        log.write("\nCreated bash script/s (*.sh), *_exp.list/s, and *_ko.list/s in OUTPUT_script_1.\n")
        log.write("\n> Now you should copy the bash script to your wcEcoli/runscripts folder.\n")

def main():
    splashscreen()
    user_input()
    gene_list = create_gene_list()
    create_scripts(gene_list)

if __name__ == "__main__":
    main()
