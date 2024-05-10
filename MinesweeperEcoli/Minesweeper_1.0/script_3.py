#!/usr/bin/env python

# Author: Joshua Rees, joshua.rees@bristol.ac.uk, Life Sciences, University of Bristol
# Modified: Ioana Gherman, ig13470@bristol.ac.uk, Engineering Mathematics, University of Bristol


from sys import exit
from itertools import chain, combinations
import re
import os
import time
import pandas as pd
import pickle

# Global constants
NAME = 'mine'
JOB = f'{NAME}conquer'
POWERSETCOMBOGLOBAL = []
YELLOW_RAN = []
BLUE_RAN = []
SIMNAMESTART = len(NAME)
SIMNAMEEND = (len(JOB))
SIMNAME = JOB[SIMNAMESTART:SIMNAMEEND]
OUTPUT_FOLDER_2 = 'OUTPUT_script_2'
INPUT_FOLDER_3 = 'INPUT_script_3'
OUTPUT_FOLDER_3 = 'OUTPUT_script_3'
INPUT_FOLDER_4X = 'INPUT_script_4X'
OUTPUT_FINAL_FOLDER = 'OUTPUT_final'
USER_INPUT_TXT = "INPUT_script_1/user_input.txt"
EXPERIMENT_SCRIPT = "OUTPUT_script_3/mineinputko{}.sh"
DELETION_LOG = "OUTPUT_final/deletionlog.txt"

def splash_screen():
    """ Fancy splash screen for Lego scripts """
    pass


def already_run_check():
    ''' Checks to see if input files are present and output files have been generated, checks with the user'''

    inputfile = f'{INPUT_FOLDER_3}/divideko_endtimes.txt'
    outputfile = f'{OUTPUT_FOLDER_3}/{NAME}conquer_blue_0.sh'

    if os.path.isfile(inputfile) and os.path.isfile(outputfile):
        print("You have expected files in your INPUT and OUTPUT folders for this script.")
        print(
            "\nThis suggests you have already run this script on this input. If you run it again the order of the powersets will randomize.")
        print("\nIf you are mid-simulation, this may make the interpretation of your results incorrect.")
        response = input("\nDo you want to run this script? (yes/no)\n> ")
        if response == 'no':
            exit(0)


def interpret_results():
    """ Match the simulation results of Stage 2 to the % board and create unsorted sim results and matched results text files. """
    boardlist = f'{OUTPUT_FOLDER_2}/alldivisionsegments_codes.txt'
    endtimes = f'{INPUT_FOLDER_3}/divideko_endtimes.txt'
    matched_results = f'{OUTPUT_FOLDER_3}/matchedresults.txt'

    # Copy endtimes results to unsortedsimresults.txt
    with open(endtimes) as endtimes_file:
        result_list = [line.rstrip() for line in endtimes_file.readlines()]

    sorted_result = sorted(result_list, key=lambda x: int(x.split('\t')[0]))

    # Zip full results list and gene_list.txt
    with open(boardlist, "r+") as third_input:
        board_list = [line.rstrip() for line in third_input.readlines()]

    full_results_list = zip(board_list, sorted_result)

    # Save matched results in output directory
    with open(matched_results, "w+", newline='\n') as second_output:
        second_output.truncate()

        # Convert zip-created list of tuples, line by line into string
        for line in full_results_list:
            second_output.write('\t'.join(str(s) for s in line) + '\n')

    print('\nMatched the simulation results with % boards and saved the results in "matchedresults.txt".')


def successCheckAndVariants():
    '''Check boards for success (i.e division) using matchedresults.txt (ordered largest to smallest when generated) and save top 3 largest deletions.
        Top 3 are assigned to colours (red, yellow, blue) and used as three variants / avenues of deletion going forward.
        Creates successfulboards.txt + variants.txt.'''

    matchedresults = "OUTPUT_script_3/matchedresults.txt"
    variants = "OUTPUT_script_3/variants.txt"
    successfulboards = "OUTPUT_script_3/successfulboards.txt"
    divisionsegment_path = "OUTPUT_script_2/divisionsegment{}.txt"
    deletionlog = "OUTPUT_final/deletionlog.txt"

    # add additional names to names list if you want to increase number of variants
    # also add additional names in createScripts(), line ~367
    names = ['red_0', 'yellow_0', 'blue_0']

    boards = open(matchedresults, "r+")
    board_list = [line.rstrip() for line in boards.readlines()]
    variants_list = []
    success_list = []
    counter = 0

    for line in board_list:
        line = line.split("\t")
        if line[4]=='Divided':
            line = line[0]
            success_list.append(line)
            if counter < 3:
                variants_list.append(line)
                counter = counter + 1

    counter = 0
    for line in success_list:
        line = line.replace('.', '_')
        success_division = divisionsegment_path.format(line)
        success_division_contents = open(success_division).read()
        success_division_contents = re.sub(r"\n", r"", success_division_contents)
        line = line.replace('_', '.')
        del success_list[counter]
        new_line = line + '\t' + success_division_contents + '\n'
        success_list.insert(counter, new_line)
        counter = counter + 1

    # save succesful boards in output dir
    successfulboards_txt = open(successfulboards, "w+", newline='\n')
    successfulboards_txt.truncate()
    # convert zip created list of tuples, line by line into string
    for line in success_list:
        successfulboards_txt.write(line)
    successfulboards_txt.close()

    named_variants = zip(names, variants_list)
    # save variants in output dir
    variants_txt = open(variants, "w+", newline='\n')
    log = open(deletionlog, "a+", newline='\n')
    log.write('\nLargest initial deletion segments:\n')
    variants_txt.truncate()
    # convert zip created list of tuples, line by line into string
    for line in named_variants:
        variants_txt.write('\t'.join(str(s) for s in line) + '\n')
        log.write('\t'.join(str(s) for s in line) + '\n')
    variants_txt.close()
    log.close()

    print(
        '\nHave filtered the % boards results for successful division and largest 3 division segments, and saved in OUTPUT_script_3/successfulboards.txt + variants.txt.')
    deletionlog = "OUTPUT_final/deletionlog.txt"
    log = open(deletionlog, "a+", newline='\n')
    log.write("\nHave filtered the % boards results for successful division and largest 3 division segments.\n")
    log.close()

    return success_list

def powerset(iterable):
    ''' Daughter function of outputToLists()
		Generates a powerset, all possible unique combinations of a set, in this case our top 3 largest deletions (colour variants) and matching deletion segments
		e.g. [1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
		From this StackOverflow thread: https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements
		Which pulls from this python documentation: https://docs.python.org/3/library/itertools.html#itertools-recipes
		'''
    s = list(iterable)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def outputToLists(components, boardnames, genecodes, red_yellow_blue):
    ''' Daughter function of variantCombinations()
        Uses powerset() to generate combinations
        Records the names of % boards in the combinations, and matches with the combined deleted genes
        Outputs COLOUR_combosandgenes.txt'''

    combos_genes = "OUTPUT_script_3/{}_combosandgenes.txt"
    successfulcomponents = []
    powersetcombos_key = []
    genesdeleted = []
    powersetcombos = []

    ### This checks the list of segments that have been logic matched to the largest deletion, against a
    ###	a list of segments that produce a successfully dividing cell. Only those on both lists get kept
    for component in components:
        for board in boardnames:
            if component == board:
                successfulcomponents.append(component)
            else:
                continue

    print('Components to create powerset ', components)
    print('All successful boardnames ', boardnames)
    print('Successful components ', successfulcomponents)

    ### see http://book.pythontips.com/en/latest/enumerate.html for use of enumerate
    ### use it here to loop over the powerset of the successful components, skipping the first entry in the powerset (always = 0 components combination)
    for i, combo in enumerate(powerset(successfulcomponents), 1):
        # for each combination of matched segments, record the names of the segments, to give a conjoined name for the combination
        # this is to avoid running the same combos multiple times. It might happen that some expected combos
        # in blue or yellow are missing because they were already ran in red. However, this will be saved in YELLOW_RAN and BLUE_RAN
        if sorted(combo) in POWERSETCOMBOGLOBAL and combo:
            if red_yellow_blue.startswith('yellow'):
                YELLOW_RAN.append(combo)
            elif red_yellow_blue.startswith('blue'):
                BLUE_RAN.append(combo)

        if sorted(combo) not in POWERSETCOMBOGLOBAL and combo:
            POWERSETCOMBOGLOBAL.append(sorted(combo))
            combo = ' '.join(combo)
            powersetcombos_key.append(combo)

            # then using the seperate names of the segments
            combosplit = combo.split(" ")
            genesdeleted = []

            for item in combosplit:
                # empty strings are considered FALSE normally.
                # Flip It - Is the string empty ( item.strip() )? IF NOT makes Yes = TRUE.
                # == if the item is just an empty string, append a space, and go to next iteration of loop
                if not item.strip():
                    genesdeleted.append(" ")
                    continue

                # use the name of the segment to find the index location (using the boardnames list)
                # use the index to recall the gene codes of that segment / board
                else:
                    x = boardnames.index(item)
                    genesdeleted.append(genecodes[x])

            # record the genes
            genesdeleted_txt = ''.join(genesdeleted)
            powersetcombos.append(genesdeleted_txt)

    # combine the recorded names and recorded genes
    combosandgenes = zip(powersetcombos_key, powersetcombos)
    combosandgenes_name = combos_genes.format(red_yellow_blue)
    combosandgenes_txt = open(combosandgenes_name, "w+", newline='\n')
    combosandgenes_txt.truncate()
    # convert zip created list of tuples, line by line into string
    for i, comb in enumerate(combosandgenes):
        combosandgenes_txt.write(powersetcombos_key[i] + '\n')
        combosandgenes_txt.write(powersetcombos[i] + '\n')
    combosandgenes_txt.close()

    #Write to text file the combos that were meant to be deleted in yellow and blue but were already ran in red
    if red_yellow_blue.startswith('yellow'):
        with open('OUTPUT_script_3/yellow_already_ran.txt', 'w') as f:
            for item in YELLOW_RAN:
                f.write(str(item) + '\n')
    elif red_yellow_blue.startswith('blue'):
        with open('OUTPUT_script_3/blue_already_ran.txt', 'w') as f:
            for item in BLUE_RAN:
                f.write(str(item) + '\n')

    ##Uncomment the line bellow if it is commented. It saves the combos to a pickle file to be used by the WCM.
    ##It takes a very long time to pickle the combos because there are tens of thousands of them. Run this script
    ##overnight if you want to pickle the combos.
    #createKOlist(powersetcombos, combosandgenes_name)
    print(f'\nCreated and saved {combosandgenes_name}.')
    deletionlog = "OUTPUT_final/deletionlog.txt"
    log = open(deletionlog, "a+", newline='\n')
    log.write(f"\nCreated and saved {combosandgenes_name}.\n")
    log.close()


def variantCombinations(successlist):
    '''The Variants (three largest deletion segments) are matched with all other dividing, non-overlapping segments
        Using logic outlined below. Removing 100% of genes ends Minesweeper. Removing 90% of the genes ends this script, moving onto the next'''
    variants = "OUTPUT_script_3/variants.txt"
    top3 = open(variants, "r+")
    top3_list = [line.rstrip() for line in top3.readlines()]
    boardnames = []
    genecodes = []

    for line in successlist:
        boardname, genecode = line.split("\t")
        boardnames.append(boardname)
        genecode = re.sub(r"\n", r"", genecode)
        genecodes.append(genecode)

    for line in top3_list:
        line = line.split("\t")
        red_yellow_blue = line[0]
        boardcode = line[1]

        if boardcode == '100.0':
            print("Final results have been written in OUTPUT_final/finalresults.txt")
            finalresults = "OUTPUT_final/finalresults.txt"
            finalresults_txt = open(finalresults, "w+", newline='\n')
            finalresults_txt.truncate()
            for line in successlist:
                finalresults_txt.write(line)
                break
            finalresults_txt.close()

            deletionlog = "OUTPUT_final/deletionlog.txt"
            log = open(deletionlog, "a+", newline='\n')
            log.write(
                "\n\nMinesweeper is finished :) Final results have been written in OUTPUT_final/finalresults.txt. ")
            log.close()

            exit()

        elif boardcode == '90.0a':
            ## If already removed 90% of the genes, going to create resulting files without running the code / simulations,
            ## to spoof the next script into working.
            explist_path = "OUTPUT_script_3/mineconquer_red_0_exp.list"
            combolist = "OUTPUT_script_3/red_0_combosandgenes.txt"
            endtimes = "INPUT_script_4X/conquerko_red_0_endtimes.txt"

            explist_path_txt = open(explist_path, "w+", newline='\n')
            explist_path_txt.truncate()
            explist_path_txt.write('conquer\n')
            explist_path_txt.write('wildtype\n')
            explist_path_txt.write('mutant\n')
            explist_path_txt.close()

            combolist_txt = open(combolist, "w+", newline='\n')
            combolist_txt.truncate()
            combolist_txt.write('\n')
            # As it's the top result, we can just read the first line from successfulboards.txt
            topline = open('OUTPUT_script_3/successfulboards.txt').readline()
            combolist_txt.write(topline)
            combolist_txt.close()

            endtimes_txt = open(endtimes, "w+", newline='\n')
            endtimes_txt.truncate()
            endtimes_txt.write('1\t11.22\nNon Essential Divided\n')
            endtimes_txt.write('2\t9.183\nNon Essential Divided\n')
            endtimes_txt.write('1\t13.89\nDNA/RNA/Protein/Metabolic NoDivision\n')
            endtimes_path_txt.close()

            print(
                "This stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
            deletionlog = "OUTPUT_final/deletionlog.txt"
            log = open(deletionlog, "a+", newline='\n')
            log.write(
                "\n\nThis stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
            log.close()

            exit()

        elif boardcode == '90.0b':
            ## If already removed 90% of the genes, going to create resulting files without running the code / simulations,
            ## to spoof the next script into working.
            explist_path = "OUTPUT_script_3/mineconquer_red_0_exp.list"
            combolist = "OUTPUT_script_3/red_0_combosandgenes.txt"
            endtimes = "INPUT_script_4X/conquerko_red_0_endtimes.txt"

            explist_path_txt = open(explist_path, "w+", newline='\n')
            explist_path_txt.truncate()
            explist_path_txt.write('conquer\n')
            explist_path_txt.write('wildtype\n')
            explist_path_txt.write('mutant\n')
            explist_path_txt.close()

            combolist_txt = open(combolist, "w+", newline='\n')
            combolist_txt.truncate()
            combolist_txt.write('\n')
            # As it's the top result, we can just read the first line from successfulboards.txt
            topline = open('OUTPUT_script_3/successfulboards.txt').readline()
            combolist_txt.write(topline)
            combolist_txt.close()

            endtimes_txt = open(endtimes, "w+", newline='\n')
            endtimes_txt.truncate()
            endtimes_txt.write('1\t11.22\nNon Essential Divided\n')
            endtimes_txt.write('2\t9.183\nNon Essential Divided\n')
            endtimes_txt.write('1\t13.89\nDNA/RNA/Protein/Metabolic NoDivision\n')
            endtimes_path_txt.close()

            print(
                "This stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
            deletionlog = "OUTPUT_final/deletionlog.txt"
            log = open(deletionlog, "a+", newline='\n')
            log.write(
                "\n\nThis stage has finished. You should run the next stage.\nYou are carrying forward one variant (red), your results have been placed in INPUT_script_4X/comboresults.txt\n")
            log.close()

            exit()

        ### Logic of matching segments starts here
        elif boardcode == '80.0a':
            components = ['80.0a', '12.5g', '12.5h', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '80.0b':
            components = ['80.0b', '12.5a', '12.5b', '5.0a', '5.0b', '5.0c', '5.0d']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '70.0a':
            components = ['70.0a', '12.5f', '12.5g', '12.5h', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '70.0b':
            components = ['70.0b', '12.5a', '12.5b', '12.5c', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '60.0a':
            components = ['60.0a', '33.0c', '12.5e', '12.5f', '12.5g', '12.5h', '5.0m', '5.0n', '5.0o', '5.0p', '5.0q',
                          '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '60.0b':
            components = ['60.0b', '33.0a', '12.5a', '12.5b', '12.5c', '12.5d', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e',
                          '5.0f', '5.0g', '5.0h']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '50.0a':
            components = ['50.0a', '33.0c', '12.5e', '12.5f', '12.5g', '12.5h', '5.0k', '5.0l', '5.0m', '5.0n', '5.0o',
                          '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '50.0b':
            components = ['50.0b', '33.0a', '12.5a', '12.5b', '12.5c', '12.5d', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e',
                          '5.0f', '5.0g', '5.0h', '5.0i', '5.0j']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '33.0a':
            components = ['33.0a', '33.0b', '33.0c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '5.0g', '5.0h',
                          '5.0i', '5.0j', '5.0k', '5.0l', '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s',
                          '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '33.0b':
            components = ['33.0b', '33.0a', '33.0c', '12.5a', '12.5b', '12.5f', '12.5g', '12.5h', '10.0d', '10.0g',
                          '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0n', '5.0o', '5.0p', '5.0q',
                          '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '33.0c':
            components = ['33.0c', '33.0a', '33.0b', '12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '5.0a', '5.0b',
                          '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l', '5.0m',
                          '5.0n']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '25.0a':
            components = ['25.0a', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '10.0c', '5.0f', '5.0g',
                          '5.0h', '5.0i', '5.0j', '5.0k', '5.0l', '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r',
                          '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '25.0b':
            components = ['25.0b', '12.5a', '12.5b', '12.5e', '12.5f', '12.5g', '12.5h', '10.0c', '5.0a', '5.0b',
                          '5.0c', '5.0d', '5.0e', '5.0k', '5.0l', '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r',
                          '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '25.0c':
            components = ['25.0c', '12.5a', '12.5b', '12.5c', '12.5d', '12.5g', '12.5h', '10.0h', '5.0a', '5.0b',
                          '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0p', '5.0q', '5.0r',
                          '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '25.0d':
            components = ['25.0d', '12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '10.0h', '5.0a', '5.0b',
                          '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l', '5.0m',
                          '5.0n', '5.0o']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5a':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '5.0c', '5.0d',
                          '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l', '5.0m', '5.0n', '5.0o',
                          '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5b':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '10.0b', '10.0c',
                          '5.0a', '5.0b', '5.0c', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l', '5.0m',
                          '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5c':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '10.0c', '10.0d',
                          '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l', '5.0m',
                          '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5d':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '5.0a', '5.0b',
                          '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0k', '5.0l', '5.0m', '5.0n', '5.0o',
                          '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5e':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '5.0a', '5.0b',
                          '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0m', '5.0n', '5.0o',
                          '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5f':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '10.0g', '10.0h',
                          '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k',
                          '5.0l', '5.0m', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5g':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5i', '12.5h', '10.0h', '10.0i',
                          '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k',
                          '5.0l', '5.0m', '5.0n', '5.0o', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '12.5h':
            components = ['12.5a', '12.5b', '12.5c', '12.5d', '12.5e', '12.5f', '12.5g', '12.5h', '5.0a', '5.0b',
                          '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l', '5.0m',
                          '5.0n', '5.0o', '5.0p', '5.0q', '5.0r']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0a':
            components = ['10.0a', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l',
                          '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0b':
            components = ['10.0b', '5.0a', '5.0b', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l',
                          '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0c':
            components = ['10.0c', '5.0a', '5.0b', '5.0c', '5.0d', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k', '5.0l',
                          '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0d':
            components = ['10.0d', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0i', '5.0j', '5.0k', '5.0l',
                          '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0e':
            components = ['10.0e', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0k', '5.0l',
                          '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0f':
            components = ['10.0f', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j',
                          '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0g':
            components = ['10.0g', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j',
                          '5.0k', '5.0l', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0h':
            components = ['10.0h', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j',
                          '5.0k', '5.0l', '5.0m', '5.0n', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0i':
            components = ['10.0i', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j',
                          '5.0k', '5.0l', '5.0m', '5.0n', '5.0o', '5.0p', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode == '10.0j':
            components = ['10.0j', '5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j',
                          '5.0k', '5.0l', '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

        elif boardcode in ['5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k',
                           '5.0l', '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']:
            components = ['5.0a', '5.0b', '5.0c', '5.0d', '5.0e', '5.0f', '5.0g', '5.0h', '5.0i', '5.0j', '5.0k',
                          '5.0l', '5.0m', '5.0n', '5.0o', '5.0p', '5.0q', '5.0r', '5.0s', '5.0t']
            outputToLists(components, boardnames, genecodes, red_yellow_blue)

def createKOlist(powersetcombos, combosandgenes_name):
    '''We use this method to create a pickle file of multiple genes to be
    deleted at once. This pickled file will be read in
    .../wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py'''

    df_names = pd.read_csv("gene_name_id_touse.csv", header=0, sep=';')
    color = combosandgenes_name.split('/')[1].split('_')[0]
    to_sim_name = []
    to_sim_ids = []
    c=0
    st = time.time()
    for l1 in powersetcombos:
        c+=1
        list_kos = l1.split(' ')
        if not list_kos[0][0].isdigit():
            l = [x.replace('\n', '').replace('\t', '').replace(' ', '') for x in list_kos if x!='\n' and x!='']
            to_sim_id = tuple()
            to_sim_name.append(tuple(l))
            for x in l:
                try:
                    df_names[df_names["Gene"] == x]["KO index"].values[0]
                except:
                    print('Gene not found ', x)

            to_sim_id = tuple([df_names[df_names["Gene"]==x]["KO index"].values[0] for x in l])
            if len(l) != len(to_sim_id):
                print('length missmatch')
            to_sim_ids.append(to_sim_id)


    create_scripts(c, color)
    with open("stage3_sims_"++color, "wb") as fp:   #Pickling
        pickle.dump(to_sim_ids, fp)
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

def create_scripts(c, combosandgenes_name):
    """ Convert template script to bash script using user input, create exp and ko txt files """

    user_input_values = open(USER_INPUT_TXT, "r+").readlines()[:10]
    value_list = [line.split()[1] for line in user_input_values]
    launchpad_path, project_path, runscripts_path, seeds_number = value_list

    experimentscript = EXPERIMENT_SCRIPT.format('_last_'+combosandgenes_name)
    exp_script = open(experimentscript, 'w', newline='\n')


    for i in range (0, c, 12):
        first_index=i
        last_index=i+11
        if last_index>=c:
            last_index=c-1

        exp_script.write(f'VARIANT="double_knockout" FIRST_VARIANT_INDEX={first_index} LAST_VARIANT_INDEX={last_index} \\\n')
        exp_script.write(f'CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \\\n')
        exp_script.write(f'SINGLE_DAUGHTERS=0 N_GENS=1 N_INIT_SIMS={seeds_number} \\\n')
        exp_script.write(f'MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 LOG_TO_DISK_EVERY=60 \\\n')
        exp_script.write(f'LAUNCHPAD_FILE="{launchpad_path}/my_launchpad.yaml" \\\n')
        exp_script.write(f'python {runscripts_path}/fireworks/fw_queue.py\n\n')

    exp_script.write('nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 40 --maxjobs_queue 100 \n')

    print(f'\nCreated bash script (*.sh), in OUTPUT_script_3.')
    print(f'\n> Now you should copy the bash script to your wcEcoli/runscripts folder.')
    print(f'\n> Remember to read stage3_sims in /wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py.')

    with open(DELETION_LOG, "a+", newline='\n') as log:
        log.write("\n> Created bash script/s in OUTPUT_script_3.\n")
        log.write("\n> Now you should copy the bash script to your wcEcoli/runscripts folder.\n")
        log.write("\n> Remember to read stage3_sims in /wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py.")



def main():
    # Display the splash screen
    splash_screen()

    # Check if the script has already been run on the current input
    already_run_check()

    # Interpret the results from Stage 2
    interpret_results()

    # Define the list of genes for which you want to create powerset combinations
    successlist = successCheckAndVariants()
    variantCombinations(successlist)


if __name__ == "__main__":
    main()