#!/usr/bin/env python

# Author: Joshua Rees, joshua.rees@bristol.ac.uk, Life Sciences, University of Bristol
# Modified: Ioana Gherman, ig13470@bristol.ac.uk, Engineering Mathematics, University of Bristol

"""
Stage 2: n deletion segments, ranging in size from 100% to 5% of the low / no essentiality genes, are generated.
Deletion segments that do not prevent division go to Stage 3.
Create deletion segments from inputko*_endtimes.txt, output segments bash scripts and text files.
Expects: Place N inputko*.txt files in INPUT_script_2 folder.
Output: deletion segments, deletion segments list, single bash script and exp / ko lists, in OUTPUT_script_2 folder
"""

# Imports
import os
import fnmatch
import re
import pandas as pd
import pickle

USER_INPUT_TXT = "INPUT_script_1/user_input.txt"
EXPERIMENT_SCRIPT = "OUTPUT_script_2/mineinputko{}.sh"
DELETION_LOG = "OUTPUT_final/deletionlog.txt"


def splashscreen():
    """ Fancy splash screen for Lego scripts """
    pass

def process_endtimes_file(n, endtimes, unsortedsimresults):
    iteration = str(n)
    nthendtimes = endtimes.format(iteration)

    with open(nthendtimes) as endtimes_file:
        endtimes_content = endtimes_file.read()
        endtimes_content = re.sub(r"(\d)\n^(\w)", r"\1\t\2", endtimes_content, flags=re.MULTILINE)

    nthunsortedsimresults = unsortedsimresults.format(iteration)
    with open(nthunsortedsimresults, "w", newline='\n') as unsorted_file:
        unsorted_file.write(endtimes_content)

    result_list = endtimes_content.split('\n')
    result_list_copy = ['No_Result'] * len(result_list)
    for line in result_list:
        results = line.split("\t")
        sim = int(results[0])
        sim -= 1
        time = results[1]
        outcome = "\t".join(results[2:])
        del result_list_copy[sim]
        result_list_copy.insert(sim, f"{time}\t{outcome}")

    return result_list_copy

def interpretResults():
    """ Match the simulation results of Stage 1 to the knocked out gene, create unsortedsimresults and matchedgeneresults.txt. """

    genelist = "INPUT_script_1/genes.txt"
    endtimes = "INPUT_script_2/inputko_endtimes.txt"
    unsortedsimresults = "OUTPUT_script_2/unsortedsimresults{}.txt"
    matchedresults = "OUTPUT_script_2/matchedgeneresults.txt"
    full_results_list = []

    numberofendtimes = len(fnmatch.filter(os.listdir('INPUT_script_2/'), '*endtimes.txt'))
    for n in range(1, numberofendtimes + 1):
        full_results_list.extend(process_endtimes_file(n, endtimes, unsortedsimresults))

    with open(genelist) as genelist_file:
        gene_list = [line.rstrip() for line in genelist_file.readlines()]

    matched_results_list = zip(gene_list, full_results_list)
    with open(matchedresults, "w", newline='\n') as matchedresults_file:
        for line in matched_results_list:
            matchedresults_file.write('\t'.join(str(s) for s in line) + '\n')

    print('\nHave matched the simulation results with gene codes and saved the results in OUTPUT_script_2/matchedgeneresults.txt.\n')

    deletionlog = "OUTPUT_final/deletionlog.txt"
    with open(deletionlog, "a", newline='\n') as log:
        log.write("\nHave matched the simulation results with gene codes and saved the results in OUTPUT_script_2/matchedgeneresults.txt.\n")

def createNEList():
    """ Filter the Stage 1 matched results, finding those that divided, while avoiding excluded genes.
        Creates nonessential.txt """

    matchedresults = "OUTPUT_script_2/matchedgeneresults.txt"
    nonessential = "OUTPUT_script_2/nonessential.txt"
    exclusionlist = "INPUT_script_2/exclusionlist.txt"

    with open(matchedresults) as generesults_file:
        generesults_list = [line.rstrip() for line in generesults_file.readlines()]

    non_essential_genes = [line.split('\t')[0] for line in generesults_list if line.split('\t')[2]=='Divided']

    with open(exclusionlist) as exclusion_file:
        exclusiongenes_list = [line.rstrip() for line in exclusion_file.readlines()]

    non_essential_genes = [gene for gene in non_essential_genes if gene not in exclusiongenes_list]

    with open(nonessential, "w", newline='\n') as neresults_file:
        neresults_file.write("\n".join(non_essential_genes) + '\n')

    print('\nHave filtered the results (deletions that still produce dividing cells) and saved in OUTPUT_script_2/nonessential.txt.\n')
    deletionlog = "OUTPUT_final/deletionlog.txt"
    with open(deletionlog, "a", newline='\n') as log:
        log.write("\nHave filtered the results (deletions that still produce dividing cells) and saved in OUTPUT_script_2/nonessential.txt.\n")

def outputToOwnFile(segment_n, segment_s, divisionsegment_list):
    divisionsegment_name = f'divisionsegment{str(segment_n).replace(".", "_")}{segment_s}.txt'
    output_divisionsegment = open(os.path.join('OUTPUT_script_2', divisionsegment_name), "w", newline='\n')
    output_divisionsegment.writelines(divisionsegment_list)
    output_divisionsegment.close()

    deletionlog = "OUTPUT_final/deletionlog.txt"
    with open(deletionlog, "a", newline='\n') as log:
        log.write(f"\nHave created the {divisionsegment_name} deletion segment in OUTPUT_script_2/{divisionsegment_name}\n")

def outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, trigger):
    alldivisionsegments_txt = "OUTPUT_script_2/alldivisionsegments.txt"
    alldivisionsegments_codes_txt = "OUTPUT_script_2/alldivisionsegments_codes.txt"

    if segment_n == 100:
        alldivisionsegments_codes.append(f"{segment_n}\n")
        alldivisionsegments_list.append(f"{segment_n}\n")
    else:
        alldivisionsegments_codes.append(f"{segment_n}{segment_s}\n")
        alldivisionsegments_list.append(f"{segment_n}{segment_s}\n")

    alldivisionsegments_list.extend(temp_list)
    alldivisionsegments_list.append('\n')

    if trigger == 'Yes':
        with open(alldivisionsegments_txt, "w", newline='\n') as output_alltxt:
            output_alltxt.writelines(alldivisionsegments_list)

        with open(alldivisionsegments_codes_txt, "w", newline='\n') as output_alltxt2:
            output_alltxt2.writelines(alldivisionsegments_codes)

        deletionlog = "OUTPUT_final/deletionlog.txt"
        with open(deletionlog, "a", newline='\n') as log:
            log.write("\nHave created the 56 deletion segments, see OUTPUT_script_2/alldivisionsegments.txt\n")
    else:
        pass

    return alldivisionsegments_codes, alldivisionsegments_list

def createDivisionSegments():
    nonessential = "OUTPUT_script_2/nonessential.txt"
    alldivisionsegments_list = []
    alldivisionsegments_codes = []

    neresults_gene_list = []
    with open(nonessential, "r") as neresults:
        for line in neresults:
            neresults_gene_list.append(line.strip() + ' ')

    percentages = [100, 90, 80, 70, 60, 50, 33, 25, 12.5, 10, 5]
    dict_alphabet = {i: chr(97 + i) for i in range(26)}
    dict_out = {}

    for p in percentages:
        if p == 100:
            dict_out[f"{p}a"] = neresults_gene_list
        elif p >= 50:
            dict_out[f"{p}a"] = neresults_gene_list[:int(p / 100 * len(neresults_gene_list))]
            dict_out[f"{p}b"] = neresults_gene_list[int((100 - p) / 100 * len(neresults_gene_list)):]
        else:
            len_segments = int(p * len(neresults_gene_list) / 100)
            for x in range(len(neresults_gene_list) // len_segments):
                dict_out[f"{p}{dict_alphabet[x]}"] = neresults_gene_list[x * len_segments:(x + 1) * len_segments]

    for k, temp_list in dict_out.items():
        segment_n = float(k[:-1])
        segment_s = k[-1]

        outputToOwnFile(segment_n, segment_s, temp_list)

        if segment_n == percentages[-1] and segment_s == dict_alphabet[int(100 / percentages[-1]) - 1]:
            alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'Yes')
        else:
            alldivisionsegments_codes, alldivisionsegments_list = outputToLists(segment_n, segment_s, alldivisionsegments_codes, alldivisionsegments_list, temp_list, 'No')


def createKOlist():
    '''We use this method to create a pickle file of multiple genes to be
    deleted at once. This pickled file will be read in
    .../wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py'''

    with open('OUTPUT_script_2/alldivisionsegments.txt') as f:
       lines = f.readlines()

    df_names = pd.read_csv("gene_name_id_touse.csv", header=0, sep=';')

    to_sim_name = []
    to_sim_ids = []
    c=0

    for l1 in lines:
        list_kos = l1.split(' ')
        if not list_kos[0][0].isdigit():
            l = [x.replace('\n', '').replace('\t', '').replace(' ', '') for x in list_kos if x!='\n' and not x[0].isdigit()]
            to_sim_id = tuple()
            to_sim_name.append(tuple(l))
            for x in l:
                try:
                    df_names[df_names["Gene"] == x]["KO index"].values[0]
                except:
                    print('Gene not found ', x)

            to_sim_id = tuple([df_names[df_names["Gene"]==x]["KO index"].values[0] for x in l])
            to_sim_ids.append(to_sim_id)
    with open("stage2_sims", "wb") as fp:   #Pickling
        pickle.dump(to_sim_ids, fp)

def create_scripts():
    """ Convert template script to bash script using user input, create exp and ko txt files """

    with open('OUTPUT_script_2/alldivisionsegments.txt') as f:
       lines = f.readlines()
    sims_to_run = int(len(lines)/2)

    user_input_values = open(USER_INPUT_TXT, "r+").readlines()[:10]
    value_list = [line.split()[1] for line in user_input_values]
    launchpad_path, project_path, runscripts_path, seeds_number = value_list

    experimentscript = EXPERIMENT_SCRIPT.format("")
    exp_script = open(experimentscript, 'w', newline='\n')


    for i in range (0, sims_to_run, 5):
        first_index=i
        last_index=i+4
        if last_index>=sims_to_run:
            last_index=sims_to_run-1

        exp_script.write(f'VARIANT="double_knockout" FIRST_VARIANT_INDEX={first_index} LAST_VARIANT_INDEX={last_index} \\\n')
        exp_script.write(f'CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \\\n')
        exp_script.write(f'SINGLE_DAUGHTERS=0 N_GENS=6 N_INIT_SIMS={seeds_number} \\\n')
        exp_script.write(f'MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 LOG_TO_DISK_EVERY=30 \\\n')
        exp_script.write(f'LAUNCHPAD_FILE="{launchpad_path}/my_launchpad.yaml" \\\n')
        exp_script.write(f'python {runscripts_path}/fireworks/fw_queue.py\n\n')

    exp_script.write('nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 40 --maxjobs_queue 100 \n')

    print(f'\nCreated bash script (*.sh), in OUTPUT_script_2.')
    print(f'\n> Now you should copy the bash script to your wcEcoli/runscripts folder.')
    print(f'\n> Remember to read stage2_sims in /wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py.')

    with open(DELETION_LOG, "a+", newline='\n') as log:
        log.write("\n> Created bash script/s in OUTPUT_script_2.\n")
        log.write("\n> Now you should copy the bash script to your wcEcoli/runscripts folder.\n")
        log.write("\n> Remember to read stage2_sims in /wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py.")

def main():
    splashscreen()
    interpretResults()
    createNEList()
    createDivisionSegments()
    #createKOlist()
    create_scripts()

if __name__ == "__main__":
    main()
