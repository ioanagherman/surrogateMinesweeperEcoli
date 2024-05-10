"""
Step 4 to X: The largest deletion segment determines the remaining low / no essentiality genes that have not been deleted.
According to ranked results from INPUT_script_4X/conquerko_red_0_endtimes.txt, repeat for yellow and blue.
These remaining genes are divided into eight groups, a powerset generated for these eight groups, and each combination from the powerset
individually appended to the current largest deletion combination and simulated.
If none of these simulations produces a dividing cell, the remaining genes are appended as single knockouts to the current largest deletion
combination and simulated. The individual remaining genes that don’t produce a dividing cell are temporarily excluded and a reduced remaining gene list produced.
Expects: Place N conquerko_COLOUR_N_endtimes.txt in INPUT_script_4X folder.
Output: either bash scripts and exp / ko lists in OUTPUT_script_4X folder OR final results in OUTPUT_final
"""

from sys import exit
from itertools import chain, combinations
import re
from datetime import datetime
import math
import os
import glob
from typing import List, Any
import ast
import pandas as pd
import time
import pickle

NAME = 'mine'
JOB = NAME + 'gapko{}'
SIMNAMESTART = len(NAME)
SIMNAMEEND = (len(JOB))-2
SIMNAME = JOB[SIMNAMESTART:SIMNAMEEND]
USER_INPUT_TXT = "INPUT_script_1/user_input.txt"
EXPERIMENT_SCRIPT = "OUTPUT_script_4X/mineinputko{}.sh"

def splashscreen():
    """ Fancy splash screen for Lego scripts """
    pass

def alreadyRunCheck():
    rounds = []
    for file in os.listdir("INPUT_script_4X/"):
        if file!='.DS_Store' and file!='desktop.ini':
            rounds.append(int(file.split('_')[-2]))
    roundnumber = max(rounds)
    roundnumberadded = roundnumber + 1

    if roundnumberadded == roundnumber:
        print(f"The latest endtimes.txt is numbered {roundnumber}, and the latest bash file is numbered {roundnumber}." )
        print("\nThis increment of 1 suggests you have already run this script on this input. If you run it again the order of the powersets will randomise." )
        print("\nIf you are mid simulation this make the interpretation of your results incorrect." )
        response = input("\nDo you want to run this script? yes / no\n> " )
        if response == 'no':
            exit(0)
    print('Round number is ', roundnumber)
    return roundnumber


def interpretResults():
    """ Match the simulation results of Stage 3 / Stage 4 to the specific powerset combination, create unsortedsimresults and matchedresults.txt. """

    # What round are we on? Input from stage 3 = 1, input from stage 4 = n + 1, should work as the number of files decreases (i.e. red, yellow finish, leaving blue)
    previousroundnumber = alreadyRunCheck()
    roundnumber = previousroundnumber+1

    # search for round results, should work as number of files decreases
    # How many endtime files we have in the input
    filesearch = 'INPUT_script_4X/*' + str(previousroundnumber) + '_endtimes.txt'
    filesearching = glob.glob(filesearch)

    names = ['_' + name.split('/')[1].split('_')[1] + "_" + str(previousroundnumber) for name in filesearching]

    # input
    if previousroundnumber == 0:
        combolist = "OUTPUT_script_3/{}_combosandgenes.txt"
        endtimes = "INPUT_script_4X/conquerko{}_endtimes.txt"
    elif previousroundnumber > 0:
        combolist = "OUTPUT_script_4X/{}_eightsandgenes.txt"  # transition from combosandgenes to eightsandgenes
        endtimes = "INPUT_script_4X/gapko{}_endtimes.txt"
    previousroundnumber = str(previousroundnumber)

    # output
    matchedresults = "OUTPUT_script_4X/matchedresults{}.txt"
    deletionlog = "OUTPUT_final/deletionlog.txt"

    for name in names:
        nthendtimes = endtimes.format(name)

        # put the endtime results in a list
        firstinput = open(nthendtimes).read()
        firstinput = firstinput.split('\n')[:-1]
        name_variant = name[1:]
        # put the combination of segments in a list
        combolist_path = combolist.format(name_variant)
        secondinput = open(combolist_path).read()
        secondinput = re.sub(r"\w\w\w\w\w ", "", secondinput, flags=re.MULTILINE)
        secondinput = re.sub(r"\w\w\w\w ", "", secondinput, flags=re.MULTILINE)
        secondinput = re.sub(r"\w\w\w ", "", secondinput, flags=re.MULTILINE)
        secondinput = secondinput.split('\n')[:-1]

        if roundnumber==1:
            if len(firstinput) != int(len(secondinput)/2)-1:
                print('ERROR: Your input does not match the previous step simulations')
        else:
            if len(firstinput) != int(len(secondinput))-1:
                print('ERROR: Your input does not match the previous step simulations')

        # match the segments if with the results
        matched = []
        for row in firstinput:
            row = row.split('\t')
            if roundnumber==1:
                row1 = [str(row[0])] + [str(secondinput[int(row[0])*2])] + row[1:]
            else:
                row1 = [str(row[0])] + [str(secondinput[int(row[0])+1].split('\t\t')[0])] + row[1:]
            matched.append(row1)

        nthmatchedresults = matchedresults.format(name)
        secondoutput = open(nthmatchedresults, "w+", newline='\n')
        secondoutput.truncate()
        for line in matched:
            secondoutput.write('\t'.join(str(s) for s in line) + '\n')
        secondoutput.close()

    print( f'\n\nHave matched the simulation results for {names} and saved the results in OUTPUT_script_4X/matchedresults.txt.')
    log = open(deletionlog, "a+", newline='\n')
    log.write(
        f"\n\nHave matched the simulation results for {names} and saved the results in OUTPUT_script_4X/matchedresults.txt.")
    log.close()

    return names, str(roundnumber), previousroundnumber

def divided(line):
    linetocheck = line
    divided = 'Divided'
    if linetocheck.endswith(divided):
        return True
    else:
        return False

def createDividingTxt(names):
    matchedresults = "OUTPUT_script_4X/matchedresults{}.txt"
    dividingresults = "OUTPUT_script_4X/dividing{}.txt"
    deletionlog = "OUTPUT_final/deletionlog.txt"
    dividingcounter_list = []

    for name in names:
        dividingcounter = 0
        nthmatchedresults = matchedresults.format(name)
        results = open(nthmatchedresults,"r+")
        results_list = [line.rstrip() for line in results.readlines()]
        dividing_results = []
        divided = 'Divided'

        for line in results_list:
            if line.endswith(divided):
                dividing_results.append(line)
                dividingcounter = dividingcounter + 1

        # save non_essential_genes in output dir
        nthdividingresults = dividingresults.format(name)
        neresults = open(nthdividingresults,"w+", newline ='\n')
        neresults.truncate()
        for line in dividing_results:
            line = line + '\n'
            neresults.write(line)
        neresults.close()

        print(f'\nHave filtered the dividing results and saved in OUTPUT_script_4X/{nthdividingresults}.')
        log = open(deletionlog, "a+", newline='\n')
        log.write(f"\nHave filtered the dividing results and saved in OUTPUT_script_4X/{nthdividingresults}.")
        log.close()
        dividingcounter_list.append(dividingcounter)

    return dividingcounter_list

def remainingGenes(names, dividingcounter_list, previousroundnumber):
    '''Calculate remaining genes from dividing results. Three Options:
        No division :: Remaining genes = re-record prior round results because the genes will be individually appended to the largest deletion
        Division + Single KO Append Flag :: prior round was an appended single KO round, Remaining genes = reduced set that divided (when singly appended)
        Division + No Flag :: Remaining genes = of those that divided, select smallest number remaining (code starts at line 487)
        Passed to endingDecision() with dividingcounter_list to determine next step.
    '''

    dividingcombos = "OUTPUT_script_4X/dividing{}.txt"
    remaininggenestxt = "OUTPUT_script_4X/remaininggenes{}.txt"
    deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
    deletionlog = "OUTPUT_final/deletionlog.txt"
    negenes = "OUTPUT_script_2/nonessential.txt"
    # format with names (but just red, yellow, blue = name_variant2)
    singleko_trigger = 'OUTPUT_final/singleko{}.txt'
    dividing_already_ran = "OUTPUT_script_4X/ran_in_another_color/dividing{}.txt"


    previousroundnumber = int(previousroundnumber)
    if previousroundnumber == 0:
        matchingcomboandgenes = "OUTPUT_script_3/{}_combosandgenes.txt"
    elif previousroundnumber > 0:
        matchingcomboandgenes = "OUTPUT_script_4X/{}_eightsandgenes.txt"  # transition from combosandgenes to eightsandgenes

    for name in names:
        x = names.index(name)
        counter = dividingcounter_list[x]
        genecodes = []
        combonames = []
        matchedlist = []
        negenes_list = []
        negenes_deleted = []

        if counter == 0:
            #if this is the first round, don't consider this colour
            if previousroundnumber == 0:
                log = open(deletionlog, "a+", newline='\n')
                log.write(
                    f"\n\nNo dividing in-silico cells were produced in this {name} round, so we will not consider this colour further.\n")
                log.close()
                noremaininggenes = open(remaininggenestxt.format(name), "w+", newline='\n')
                noremaininggenes.write("")
                noremaininggenes.close()


                copypastdeletedgenes = deletedgenestxt.format(name)
                copypastdeletedgenestxt = open(copypastdeletedgenes, "w+", newline='\n')
                copypastdeletedgenestxt.write("")
                copypastdeletedgenestxt.close()

            else:
                # we have to re-record the previous round results such that genes are individually appended to the largest deletion in the next step
                name_variant = name[1:]
                pastname_components = name_variant.split('_')
                pastname_name = pastname_components[0]
                pastname_number = pastname_components[1]
                pastname_number = int(pastname_number)
                pastname_number = pastname_number - 1
                pastname_number = str(pastname_number)
                pastname_full = '_' + pastname_name + '_' + pastname_number


                # As no new division, copy remaining / deleted genes from past file to current file
                nthremaininggenes = remaininggenestxt.format(pastname_full)
                remaininggenes_results = open(nthremaininggenes).read()
                copypastremaininggenes = remaininggenestxt.format(name)
                copypastremaininggenestxt = open(copypastremaininggenes, "w+", newline='\n')
                copypastremaininggenestxt.write(remaininggenes_results)
                copypastremaininggenestxt.close()

                # deletedgenes_results = past results, copy past deleted genes = copying past results into current results
                nthdeletedgenes = deletedgenestxt.format(pastname_full)
                deletedgenes_results = open(nthdeletedgenes).read()
                copypastdeletedgenes = deletedgenestxt.format(name)
                copypastdeletedgenestxt = open(copypastdeletedgenes, "w+", newline='\n')
                copypastdeletedgenestxt.write(deletedgenes_results)
                copypastdeletedgenestxt.close()

                log = open(deletionlog, "a+", newline='\n')
                log.write(
                    f"\n\nNo dividing in-silico cells were produced in this {name} round, re-recording {name} last round's results.\n")
                log.close()
                continue
        elif counter > 0:
            # this is the case in which some deletion combination divided. There are 3 options here:
            #1. This is round 0 and we append the 8 segments. This is equivalent to option 3.
            #2. This is not round 0 and a deletion including some of the 8 combinations divided. In this  case we update the largest deletion and the remaining genes based on this largest deletion.
            #3. We already ran single KO appended to the largest deletion, some divided, so the remaining genes list is updated
            name_variant2 = name[1:-2]
            nthsingleko_trigger = singleko_trigger.format(name_variant2)
            name_variant = name[1:]

            if os.path.isfile(nthsingleko_trigger):
                # single gene knockouts appended to the largest deletion, have been run, and produced division
                # this happened because the prior eights simulation round produced no division
                # so create deleted genes from prior deleted genes (which were passed from eights to singleko to here) (i.e. what the single kos were appended to)
                pastname_components = name_variant.split('_')
                pastname_name = pastname_components[0]
                pastname_number = pastname_components[1]
                pastname_number = int(pastname_number)
                pastname_number = pastname_number - 1
                pastname_number = str(pastname_number)
                pastname_full = '_' + pastname_name + '_' + pastname_number
                # deletedgenes_results = past results, copy past deleted genes = copying past results into current results
                nthdeletedgenes = deletedgenestxt.format(pastname_full)
                deletedgenes_results = open(nthdeletedgenes).read()
                copypastdeletedgenes = deletedgenestxt.format(name)
                copypastdeletedgenestxt = open(copypastdeletedgenes, "w+", newline='\n')
                copypastdeletedgenestxt.write(deletedgenes_results)
                copypastdeletedgenestxt.close()

                # and create a reduced set of remaining genes (i.e. combining the appended single kos)
                # temporarily exclude genes that did not divide when appended singly to largest deletion
                # load last rounds (singlekos) eights and combos
                # split to genecodes and code names
                # match dividing sims to genecodes via code names

                # open dividing_COLOUR_N (get combination names that divided)
                dividingtxt = dividingcombos.format(name)
                dividingallcolumns = open(dividingtxt, "r+")
                dividingallcolumns_list = [line.rstrip() for line in dividingallcolumns.readlines()]


                # match only combination names that divided, with gene codes, via lists produced from combosandgenes/eightsandgenes
                reducedremaininggenes = []
                for line in dividingallcolumns_list:
                    results = line.split("\t")
                    genestring = results[1]
                    reducedremaininggenes.append(genestring)

                # get a list of all the genes from matchedresults from previous round
                all_genes_prev_round = open('OUTPUT_script_4X/matchedresults'+name+'.txt', "r+")
                all_genes_prev_round = [line.rstrip() for line in all_genes_prev_round.readlines()]
                genecodes = [x.split('\t')[1] for x in all_genes_prev_round]

                # create a list of genes that did not divide when appended to largest deletion
                nondividingappended = [item for item in genecodes if item not in reducedremaininggenes]
                print('Non dividing genes ', len(nondividingappended))
                remaininggenesn = len(reducedremaininggenes)
                print('Dividing genes ', remaininggenesn)

                log = open(deletionlog, "a+", newline='\n')
                log.write(
                    f"\n In {name_variant} there were {remaininggenesn} genes found that could be appended to largest deletion and still divide, and {nondividingappended} genes that couldn't.\n")
                log.write(f"The reduced remaining genes are: {reducedremaininggenes}\n")
                log.write(
                    f"The temporarily excluded genes (potential conditional essentials) are: {nondividingappended}\n")
                log.close()

                remaininggenes = reducedremaininggenes

                nthremaininggenes = remaininggenestxt.format(name)
                remaininggenes_results = open(nthremaininggenes, "w+", newline='\n')
                remaininggenes_results.truncate()
                for line in remaininggenes:
                    line = line + ' '
                    remaininggenes_results.write(line)
                remaininggenes_results.close()

                continue

            elif not os.path.isfile(nthsingleko_trigger):
                # best deletion is selected, written to new deleted gene list
                # remaining genes is generated by comparing non essential genes with the new deleted gene list
                matchingcombosandgenestxt = matchingcomboandgenes.format(name_variant)
                combonames, genecodes = get_matched_combos_genes(matchingcombosandgenestxt, name_variant, previousroundnumber, combonames, genecodes)

                dividingtxt = dividingcombos.format(name)
                bestresult, longestlength = find_longest_dividing(dividingtxt, combonames, genecodes, matchedlist)
                print('Longest dividing ', bestresult, longestlength)

                # Create files with the combinations from OUTPUT_script_3/color_already_ran.txt (that were already ran in another color)
                # that divided when ran in another color. This is to make sure we don't miss a maximum deletion.
                # This assumes that in stage3 we first create the red combos, then the yellow combos, and then the blue combos.
                if os.path.isfile('OUTPUT_script_3/{}_already_ran.txt'.format(name_variant.split('_')[0])) and previousroundnumber == 0:
                    color_already_ran = open('OUTPUT_script_3/{}_already_ran.txt'.format(name_variant.split('_')[0]), "r+")
                    # color_already_ran has a list of tuples of combinations, where each tuple is on its own row. Read this into a list of lists.
                    color_already_ran_list = [line.rstrip() for line in color_already_ran.readlines()]
                    os.makedirs('OUTPUT_script_4X/ran_in_another_color', exist_ok=True)

                    if name_variant.split('_')[0] == 'yellow':
                        # check in red
                        to_check = open(dividingcombos.format('_red_0'), "r+")
                        dividingallcolumns_list_y = [line.rstrip() for line in to_check.readlines()]
                        combos_red = [line.split('\t')[1] for line in dividingallcolumns_list_y]
                    elif name_variant.split('_')[0] == 'blue':
                        # check in both red and yellow
                        to_check1 = open(dividingcombos.format('_red_0'), "r+")
                        dividingallcolumns_list1 = [line.rstrip() for line in to_check1.readlines()]
                        combos_red = [line.split('\t')[1] for line in dividingallcolumns_list1]
                        to_check2 = open(dividingcombos.format('_yellow_0'), "r+")
                        dividingallcolumns_list2 = [line.rstrip() for line in to_check2.readlines()]
                        combos_yellow = [line.split('\t')[1] for line in dividingallcolumns_list2]

                    #write to appropriate files
                    combonames_from_red, genecodes_from_red, combonames_from_yellow, genecodes_from_yellow  = [], [], [], []
                    for line in color_already_ran_list:
                        line = line.split("\t")
                        combostring = ' '.join(el for el in ast.literal_eval(line[0]))

                        if name_variant.split('_')[0] == 'yellow' and combostring in combos_red:
                            dividingallcolumns = open(dividing_already_ran.format(name+'_fom_red'), "a+", newline='\n')
                            index = combos_red.index(combostring)
                            dividingallcolumns.write(dividingallcolumns_list_y[index].split('/t')[0]+ '\n')

                            combonames_from_red.append(combostring)
                            all_combos_and_genes = [line.rstrip() for line in open('OUTPUT_script_3/red_0_combosandgenes.txt', "r+").readlines()]
                            combos = [x for i, x in enumerate(all_combos_and_genes) if i%2==0]
                            genecodes_red = [x for i, x in enumerate(all_combos_and_genes) if i%2!=0]
                            genecodes_from_red.append(genecodes_red[combos.index(combostring)])
                            folder_name_from_red = dividing_already_ran.format(name + '_fom_red')


                        elif name_variant.split('_')[0] == 'blue':
                            if combostring in combos_red:
                                dividingallcolumns = open(dividing_already_ran.format(name+'_from_red'), "a+", newline='\n')
                                index = combos_red.index(combostring)
                                dividingallcolumns.write(dividingallcolumns_list1[index].split('/t')[0]+ '\n')

                                combonames_from_red.append(combostring)
                                all_combos_and_genes = [line.rstrip() for line in
                                                        open('OUTPUT_script_3/red_0_combosandgenes.txt',
                                                             "r+").readlines()]
                                combos = [x for i, x in enumerate(all_combos_and_genes) if i % 2 == 0]
                                genecodes_red = [x for i, x in enumerate(all_combos_and_genes) if i % 2 != 0]
                                genecodes_from_red.append(genecodes_red[combos.index(combostring)])
                                folder_name_from_red = dividing_already_ran.format(name + '_from_red')

                            elif combostring in combos_yellow:
                                dividingallcolumns = open(dividing_already_ran.format(name+'_from_yellow'), "a+", newline='\n')
                                index = combos_yellow.index(combostring)
                                dividingallcolumns.write(dividingallcolumns_list2[index].split('/t')[0]+ '\n')

                                combonames_from_yellow.append(combostring)
                                all_combos_and_genes = [line.rstrip() for line in
                                                        open('OUTPUT_script_3/yellow_0_combosandgenes.txt',
                                                             "r+").readlines()]
                                combos = [x for i, x in enumerate(all_combos_and_genes) if i % 2 == 0]
                                genecodes_yellow = [x for i, x in enumerate(all_combos_and_genes) if i % 2 != 0]
                                genecodes_from_yellow.append(genecodes_yellow[combos.index(combostring)])
                                folder_name_from_yellow = dividing_already_ran.format(name + '_from_yellow')

                    if name_variant.split('_')[0] == 'yellow':
                        bestresult_other_color, longestlength_other_color = find_longest_dividing(folder_name_from_red, combonames_from_red, genecodes_from_red, [])
                        if len(bestresult_other_color.split('\t')[1].split(' ')) > len(bestresult.split('\t')[1].split(' ')):
                            bestresult = bestresult_other_color
                            longestlength = longestlength_other_color
                            print('There was a longer result that was not run in yellow because it was already ran in red. This was the result: ', bestresult_other_color, longestlength_other_color)
                    elif name_variant.split('_')[0] == 'blue':
                        bestresult_other_color_yellow, longestlength_other_color_yellow = find_longest_dividing(folder_name_from_yellow, combonames_from_yellow, genecodes_from_yellow, [])
                        bestresult_other_color_red, longestlength_other_color_red = find_longest_dividing(folder_name_from_red, combonames_from_red, genecodes_from_red, [])
                        if len(bestresult_other_color_yellow.split('\t')[1].split(' ')) > len(bestresult.split('\t')[1].split(' ')):
                            bestresult = bestresult_other_color_yellow
                            longestlength = longestlength_other_color_yellow
                            print('There was a longer result that was not run in blue because it was already ran in yellow. This was the result: ', bestresult_other_color_yellow, longestlength_other_color_yellow)
                        if len(bestresult_other_color_red.split('\t')[1].split(' ')) > len(bestresult.split('\t')[1].split(' ')):
                            bestresult = bestresult_other_color_red
                            longestlength = longestlength_other_color_red
                            print('There was a longer result that was not run in blue because it was already ran in red. This was the result: ', bestresult_other_color_red, longestlength_other_color_red)



                # compare best result vs OUTPUT_script_1/gene_list.txt
                negenes_txt = open(negenes, "r+")
                full_negenes_list = [line.rstrip() for line in negenes_txt.readlines()]
                for line in full_negenes_list:
                    line = line.split("\t")
                    line = line[0]
                    line = re.sub(r"\n", r"", line)
                    negenes_list.append(line)
                bestgeneresult = bestresult.split("\t")
                bestgeneresult = bestgeneresult[1]
                bestgeneresult = re.sub(r"\n", "", bestgeneresult)
                bestresult_list = bestgeneresult.split(' ')

                for gene in negenes_list:
                    for unorderedgene in bestresult_list:
                        if gene == unorderedgene:
                            negenes_deleted.append(gene)

                remaininggenes = [item for item in negenes_list if item not in negenes_deleted]
                remaininggenesn = len(remaininggenes)
                negenes_deletedn = len(negenes_deleted)

                log = open(deletionlog, "a+", newline='\n')
                log.write(
                    f"\n\nA {name} combination deleted {negenes_deletedn} genes, leaving {remaininggenesn} remaining genes.\n")
                # log.write(f"The largest number of genes deleted was by combination: {combostring}\n")
                # cannot get to work > i.e. Blue_0 prints 12.5h 25b 12.5b (sim 26), but the genes are from 25b 12.5g 12.5f (sim 18)
                log.write(f"The genes deleted: {negenes_deleted}\n")
                log.write(f"The remaining genes are: {remaininggenes}\n")
                log.close()

                nthremaininggenes = remaininggenestxt.format(name)
                remaininggenes_results = open(nthremaininggenes, "w+", newline='\n')
                remaininggenes_results.truncate()
                for line in remaininggenes:
                    line = line + ' '
                    remaininggenes_results.write(line)
                remaininggenes_results.close()

                nthdeletedgenes = deletedgenestxt.format(name)
                deletedgenes_results = open(nthdeletedgenes, "w+", newline='\n')
                deletedgenes_results.truncate()
                for line in negenes_deleted:
                    line = line + ' '
                    deletedgenes_results.write(line)
                deletedgenes_results.close()

def get_matched_combos_genes(matchingcombosandgenestxt, name_variant, previousroundnumber, combonames, genecodes):
    matchingcombosandgenesallcolumns = open(matchingcombosandgenestxt, "r+")
    matchingcombosandgenesallcolumns_list = [line.rstrip() for line in
                                             matchingcombosandgenesallcolumns.readlines()]

    # open combosandgenes / eightsandgenes (get combination name and gene codes)
    for line in matchingcombosandgenesallcolumns_list:
        if not line.strip():
            continue
        else:
            if previousroundnumber == 0:
                if line[0].isdigit():
                    combonames.append(line)
                else:
                    genecodes.append(line)
            else:
                line = line.split('\t')
                combo = line[0]
                combonames.append(combo)
                genecode = line[1] + line[2]
                genecodes.append(genecode)
    return combonames, genecodes

def find_longest_dividing(dividingtxt, combonames, genecodes, matchedlist):
    # open dividing_COLOUR_N (get combination names that divided)
    dividingallcolumns = open(dividingtxt, "r+")
    dividingallcolumns_list = [line.rstrip() for line in dividingallcolumns.readlines()]

    # match only combination names that divided, with gene codes, via lists produced from combosandgenes/eightsandgenes
    for line in dividingallcolumns_list:
        results = line.split("\t")
        combostring = results[1]
        genelocation = combonames.index(combostring)
        genestring = genecodes[genelocation]
        matchedresult = combostring + '\t' + genestring + '\n'
        matchedlist.append(matchedresult)

    # get best (i.e. most genes deleted) out of matched list
    longestlength = 0
    for line in matchedlist:
        line = line.split("\t")
        combostring = line[0]
        combostring = re.sub(r"\n", r"", combostring)
        genestring = line[1]
        genestring = re.sub(r"\n", r"", genestring)
        templist = genestring.split(' ')
        # covert list to set then list = list(set(t))
        # this removes order, but order should come from remaining genes anyway / negenes.txt
        templist = list(set(templist))
        lencount = len(templist)
        genestring = ' '.join(templist)
        if lencount > longestlength:
            bestresult = combostring + '\t' + genestring + '\n'
            longestlength = lencount
    return bestresult, longestlength

def powerset(iterable):
	''' Daughter function of outputToLists()
		Generates a powerset, all possible unique combinations of a set, in this case our top 3 largest deletions (colour variants) and matching deletion segments
		e.g. [1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
		From this StackOverflow thread: https://stackoverflow.com/questions/464864/how-to-get-all-possible-combinations-of-a-list-s-elements
		Which pulls from this python documentation: https://docs.python.org/3/library/itertools.html#itertools-recipes
		'''

	s = list(set(iterable))
	return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def outputToLists(eightsnames_list, deletedgenes, name_variant, roundnumbername):
    ''' Daughter function of eightPanelGroupingsGeneration()
        Uses powerset() to generate combinations
        Records the names of 8 groups in the combinations, and matches with the genes within group
        Outputs COLOUR_eightsandgenes.txt (complete powerset)
    '''

    # Have to code around naming change that occurs here (number change, and _ (underscore) to handle):
    # eightssegments uses name_variant (== previousroundnumber)
    # eightsandgenes uses roundnumbername (== current round number)
    # combosandgenes / eightsandgenes / eightssegments are only files that use name as prefix
    # not insert / suffix, drop the first _ from _COLOUR_N = COLOUR_N using name_variant

    # produce power sets of eightsnames
    # use to generate deletion lines by opening eightssegment (containing genes) with associated eightsnames
    eightssegment = "OUTPUT_script_4X/eightssegments/{}remaining_{}.txt"  # name_variant, eightsname
    eightsandgenes = "OUTPUT_script_4X/{}_eightsandgenes.txt"  # roundnumbername
    deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
    deletionlog = "OUTPUT_final/deletionlog.txt"
    powersetcombos_key = []
    genesdeleted = []
    deletedgenes_list = []
    powersetcombos = []

    ### see http://book.pythontips.com/en/latest/enumerate.html for use of enumerate
    ### use it here to loop over the powerset of the eight groups, skipping the first entry in the powerset (always = 0 components combination)
    for i, combo in enumerate(powerset(eightsnames_list), 1):
        # for each combination of 8 groups, record the names of the groups, to give a conjoined name for the combination
        combo = ' '.join(combo)
        powersetcombos_key.append(combo)

        # then using the seperate names of the groups
        combosplit = combo.split(" ")
        genesdeleted = []
        for item in combosplit:
            # empty strings are considered FALSE normally.
            # Flip It - Is the string empty ( item.strip() )? IF NOT makes Yes = TRUE.
            # == if the item is just an empty string, append a space, and go to next iteration of loop
            if not item.strip():
                genesdeleted.append(" ")
                continue

            # look up genes from eightssegment, using eightsname
            else:
                eightssegment_name = eightssegment.format(name_variant, item)
                eightssegmenttxt = open(eightssegment_name).read()
                genesdeleted.append(eightssegmenttxt)

        genesdeleted_txt = ''.join(genesdeleted)
        powersetcombos.append(genesdeleted_txt)


    name = '_' + name_variant
    nthdeletedgenes = deletedgenestxt.format(name)
    with open(nthdeletedgenes, 'r') as deletedgenes_txt:
        deletedgenes_txtall = deletedgenes_txt.read().replace('\n', '')
    linecount = 1

    for line in powersetcombos:
        if linecount == 1:
            deletedgenes_list.append('')
            linecount = linecount + 1
        elif linecount > 1:
            for line in deletedgenes_txtall:
                deletedgenes_list.append(deletedgenes_txtall)

    # PRODUCE eightsandgenes using roundnumbername RATHER than name_variant
    # CHANGE to eightsandgenes
    eightsandgenes_list = zip(powersetcombos_key, deletedgenes_list, powersetcombos)
    eightsandgenes_name = eightsandgenes.format(roundnumbername)
    eightsandgenes_txt = open(eightsandgenes_name, "w+", newline='\n')
    eightsandgenes_txt.truncate()
    # convert zip created list of tuples, line by line into string
    for line in eightsandgenes_list:
        eightsandgenes_txt.write('\t'.join(str(s) for s in line) + '\n')
    eightsandgenes_txt.close()

    tempscript = open(eightsandgenes_name).read()
    tempscript = re.sub(r",'", ", '", tempscript, flags=re.MULTILINE)
    tempscript = re.sub(r",[\s]*'", ", '", tempscript, flags=re.MULTILINE)
    eightsandgenes_txt = open(eightsandgenes_name, 'w+', newline='\n')
    eightsandgenes_txt.truncate()
    eightsandgenes_txt.write(tempscript)

    #createKOlist(eightsandgenes_name, name, roundnumbername)
    #createKOlist(eightsandgenes_name, name, roundnumbername)
    print(f'\nCreated and saved {eightsandgenes_name}.')
    log = open(deletionlog, "a+", newline='\n')
    log.write(f"\nCreated and saved {eightsandgenes_name}.")
    log.close()


def eightsegmentwrite(name_variant, eights_name, eightsegment_list):
    ''' Daughter function of eightPanelGroupingsGeneration()
        Outputs COLOUR_remaining_N.txt (each of the 8 groups individually)'''

    eightssegment = "OUTPUT_script_4X/eightssegments/{}remaining_{}.txt"  # name_variant, eightsname
    eightssegmenttxt = eightssegment.format(name_variant, eights_name)
    output_eightsegment = open(eightssegmenttxt, "w+", newline='\n')
    output_eightsegment.truncate()
    for line in eightsegment_list:
        output_eightsegment.write(line+' ')
    output_eightsegment.close()

def split_list(lst, m):
    n = len(lst)
    quotient = n // m
    remainder = n % m
    result = []
    start = 0
    for i in range(m):
        end = start + quotient + (i < remainder)
        result.append(lst[start:end])
        start = end
    return result

def eightPanelGroupingsGeneration(remaininggenes, deletedgenes, name, roundnumbername):
    ''' Daughter function of endingDecision()
        Calls eightsegmentwrite() to save each of the eight groups of genes individually
        Calls outputToLists() to save complete powerset.
    '''

    eightssegment = "OUTPUT_script_4X/eightssegments/{}remaining_{}.txt"  # name_variant, eightsname
    deletionlog = "OUTPUT_final/deletionlog.txt"
    clonelist = remaininggenes
    name_variant = name[1:]
    roundnumbername_variant = roundnumbername[1:]
    name_variant2 = name[1:-2]

    whole = len(remaininggenes)  # math.trunc used to round down > prevent overlapping gene knockouts

    nofgroups = 8
    # number of groups and their names in a list for def(outputToLists)
    eightsname_list = []
    list_segments = split_list(clonelist, nofgroups)
    for x in range(1, nofgroups + 1):
        eightsname_list.append(str(x))
        eightsegmentwrite(name_variant, str(x), list_segments[x-1])

    lastcreated = eightssegment.format(name_variant, str(x))
    print(f'\n\nLast created 1/8ths deletion segments of remaining genes = {lastcreated}.')
    log = open(deletionlog, "a+", newline='\n')
    log.write(f"\n\nLast created 1/8ths deletion segments of remaining genes = {lastcreated}.")
    log.close()

    outputToLists(eightsname_list, deletedgenes, name_variant, roundnumbername_variant)


def runappendsingleKOS(remaininggenes_list, deletedgenes_list, name, roundnumbername):
    ''' Daughter function of endingDecision()
        single gene knockouts are appended to the largest deletion with a successful dividing cell, to be run
        this is due to the prior eights simulation round producing no division (see createDividingTxt() and remainingGenes())
        Ouputs a singlekos version of COLOUR_eightsandgenes.txt
    '''

    eightsandgenes = "OUTPUT_script_4X/{}_eightsandgenes.txt"  # roundnumbername
    deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
    deletionlog = "OUTPUT_final/deletionlog.txt"
    deletedgenes = []
    remaininggenes = []
    codename = []

    # produce deleted genes as a text block per line (not single gene per line aka from the list) (with starting line set to blank to lineup with powerset format)
    nthdeletedgenes = deletedgenestxt.format(name)
    with open(nthdeletedgenes, 'r') as deletedgenes_txt:
        deletedgenes_txtall = deletedgenes_txt.read().replace('\n', '')

    remaininggenes_list_len = len(remaininggenes_list) + 1

    linecount = 1
    for line in range(0, remaininggenes_list_len):
        if linecount == 1:
            deletedgenes.append('')
            linecount = linecount + 1
        elif linecount > 1:
            for line in deletedgenes_txtall:
                deletedgenes.append(deletedgenes_txtall)

    # produce remaining genes list (with starting line set to blank to lineup with powerset format)
    linecount = 1
    linenumber = 0
    for line in range(0, remaininggenes_list_len):
        if linecount == 1:
            remaininggenes.append('')
            linecount = linecount + 1
        elif linecount > 1:
            gene = remaininggenes_list[linenumber]
            remaininggenes.append(gene)
            linenumber = linenumber + 1

    # produce names for each line
    linecount = 1
    linenumber = 0
    for line in remaininggenes:
        if linecount == 1:
            remaininggenes.append('')
            linecount = linecount + 1
        elif linecount > 1:
            gene = remaininggenes[linenumber]
            gene = re.sub(r"'", "", gene)
            gene = re.sub(r",", "", gene)
            codename.append(gene)
            linenumber = linenumber + 1

    # PRODUCE eightsandgenes using roundnumbername RATHER than name_variant
    # CHANGE to eightsandgenes
    roundnumbername_variant = roundnumbername[1:]
    eightsandgenes_list = zip(codename, deletedgenes, remaininggenes)
    eightsandgenes_name = eightsandgenes.format(roundnumbername_variant)
    eightsandgenes_txt = open(eightsandgenes_name, "w+", newline='\n')
    eightsandgenes_txt.truncate()
    # convert zip created list of tuples, line by line into string
    for line in eightsandgenes_list:
        eightsandgenes_txt.write('\t'.join(str(s) for s in line) + '\n')
    eightsandgenes_txt.close()

    tempscript = open(eightsandgenes_name).read()
    tempscript = re.sub(r",'", ", '", tempscript, flags=re.MULTILINE)
    tempscript = re.sub(r",[\s]*'", ", '", tempscript, flags=re.MULTILINE)
    eightsandgenes_txt = open(eightsandgenes_name, 'w+', newline='\n')
    eightsandgenes_txt.truncate()
    eightsandgenes_txt.write(tempscript)

    #createKOlist(eightsandgenes_name, name, roundnumbername)
    print(f'\nCreated and saved a (singlekos) version of an eightsandgenes: {eightsandgenes_name}.')
    log = open(deletionlog, "a+", newline='\n')
    log.write(f"\nCreated and saved a (singlekos) version of an eightsandgenes: {eightsandgenes_name}.")
    log.close()


def endingDecision(names, dividingcounter_list, JOB, SIMNAME, previousroundnumber, roundnumber):
    ''' Given the outcome of createDividingTxt() and remainingGenes() determine next step for remaining variants (e.g. Red, Yellow, Blue).
        Four Options per variant:
        If remaining genes =< 16 in last simulation round, triggering final round, record final results (if all three variants complete = end of Minesweeper)
        If conducted singlegene kos for remaining genes, individually appended to largest deletion, in last simulation round
            and none divided, record final results
        If conducted a normal eight group powerset deletion in last simulation round, and none divided, start an appended single kos round
        If conducted a normal eights group powerset deletion in last simulation round, and some divided,
            continue to another round / final round of eights (depending on n of remaining genes)

        Call eightPanelGroupingsGeneration(), runandcombinesingleKOs(), createScripts() as needed to:
        Outputs results in OUTPUT_final and next simulation round files in /bashexpkofiles folder
        '''

    remaininggenestxt = "OUTPUT_script_4X/remaininggenes{}.txt"
    deletedgenestxt = "OUTPUT_script_4X/deletedgenes{}.txt"
    # format with names (but just red, yellow, blue = name_variant2)
    singleko_trigger = 'OUTPUT_final/singleko{}.txt'
    finalround_trigger = 'OUTPUT_final/finalround{}.txt'
    finaloutput = 'OUTPUT_final/{}_finalresult.txt'
    deletionlog = "OUTPUT_final/deletionlog.txt"

    # name == previousroundnumber_name
    roundnumber_names = []

    for name in names:
        colourextract = name.split('_')[1]
        if colourextract == 'red':
            roundnumber_names.append('_red_' + roundnumber)
        elif colourextract == 'yellow':
            roundnumber_names.append('_yellow_' + roundnumber)
        elif colourextract == 'blue':
            roundnumber_names.append('_blue_' + roundnumber)

    for name in names:
        x = names.index(name)
        counter = dividingcounter_list[x]

        name_variant2 = name[1:-2]
        ### use name_variant2 > as colour only for triggers
        nthsingleko_trigger = singleko_trigger.format(name_variant2)
        nthfinalround_trigger = finalround_trigger.format(name_variant2)
        nthfinaloutput = finaloutput.format(name_variant2)

        # get final results from remaininggenes / deleted genes
        remaininggenes_list = []
        nthremaininggenes = remaininggenestxt.format(name)
        remaininggenes_txt = open(nthremaininggenes, "r+")
        for line in remaininggenes_txt:
            linecomponents = line.split(' ')
            for components in linecomponents:
                remaininggenes_list.append(components)
        remaininggenes_list = list(filter(None, remaininggenes_list))

        deletedgenes_list = []
        nthdeletedgenes = deletedgenestxt.format(name)
        deletedgenes_txt = open(nthdeletedgenes, "r+")
        for line in deletedgenes_txt:
            linecomponents = line.split(' ')
            for components in linecomponents:
                deletedgenes_list.append(components)
        deletedgenes_list = list(filter(None, deletedgenes_list))

        remaininggenesn = len(remaininggenes_list)
        print(f"Name: {name} Dividing counter: {counter} remaininggenes: {remaininggenes_list} n of remaininggenes: {remaininggenesn}")

        ## if reached remaining genes =< 16 in last simulation round, triggering final round, record final results
        if os.path.isfile(nthfinalround_trigger):
            finaloutputtxt = open(nthfinaloutput, "w+", newline='\n')
            finaloutputtxt.write(f"{name_variant2} deleted {deletedgenes_list}\n")
            finaloutputtxt.write(f"{name_variant2} did not delete {remaininggenes_list}\n")
            finaloutputtxt.close()

            print(
                f'\n\n{name_variant2} has finished its final round. See OUTPUT_final/{name_variant2}_finalresult.txt for result.')
            log = open(deletionlog, "a+", newline='\n')
            log.write(
                f"\n\n{name_variant2} has finished its final round. See OUTPUT_final/{name_variant2}_finalresult.txt for result.")
            log.close()

            continue

        ## if conducted singlegene kos for remaining genes individually appended to largest deletion in last simulation round, and non divided, record final results
        elif (os.path.isfile(nthsingleko_trigger) and counter == 0) or (counter==0 and roundnumber==0):
            finaloutputtxt = open(nthfinaloutput, "w+", newline='\n')
            finaloutputtxt.write(f"{name_variant2} deleted {deletedgenes_list}\n")
            finaloutputtxt.write(f"{name_variant2} deleted {remaininggenes_list}\n")
            finaloutputtxt.close()

            print(
                f'\n\n{name_variant2} has completed appended single ko sims, and with no dividing, has finished its final round.\nSee OUTPUT_final/{name_variant2}_finalresult.txt for result.\n')
            log = open(deletionlog, "a+", newline='\n')
            log.write(
                f"\n\n{name_variant2} has completed appended single ko sims, and with no dividing, has finished its final round.\nSee OUTPUT_final/{name_variant2}_finalresult.txt for result.\n")
            log.close()

            continue

        ## if conducted a normal eight group powerset deletion in last simulation round, and none divided, start an appended single kos round
        ## also write a single ko trigger to interpret the results differently in the following round
        elif not os.path.isfile(nthsingleko_trigger) and counter == 0:
            singleko_triggertxt = open(nthsingleko_trigger, "w+", newline='\n')
            singleko_triggertxt.write("\n")
            singleko_triggertxt.close()

            ### use roundnumbername > as future
            roundnumbername = roundnumber_names[x]
            #### use name > as past / current results
            runappendsingleKOS(remaininggenes_list, deletedgenes_list, name, roundnumbername)

            ### use roundnumbername > as future
            # !!!! gets input from {}eightsandgenes.txt / {}singlekopowersets.txt OR list returned from prior function?

            print(
                f'\n{name_variant2} produced 0 divisions and has more than 16 genes remaining, so progressing to appended single ko sims.\n')
            log = open(deletionlog, "a+", newline='\n')
            log.write(
                f"\n{name_variant2} produced 0 divisions and has more than 16 genes remaining, so progressing to appended single ko sims.\n")
            log.close()

            continue

        # if conducted a normal eights deletion in last simulation round, and some divided, continue to another round / final round of eights (depending on n of remaining genes)
        elif counter > 0:
            remaininggenesn = len(remaininggenes_list)
            if remaininggenesn <= 16:
                # write finalround trigger into folder
                finalround_triggertxt = open(nthfinalround_trigger, "w+", newline='\n')
                finalround_triggertxt.write("\n")
                finalround_triggertxt.close()

                # delete singleko_trigger if exists
                if os.path.isfile(nthsingleko_trigger):
                    os.remove(nthsingleko_trigger)

                ### use roundnumbername > as future
                roundnumbername = roundnumber_names[x]
                #### use name > as past / current results
                eightPanelGroupingsGeneration(remaininggenes_list, deletedgenes_list, name, roundnumbername)

                ### use roundnumbername > as future
                # !!!! gets input from {}eightsandgenes.txt / {}singlekopowersets.txt OR list returned from prior function?

                print(f'\n{name_variant2} has less than 16 genes remaining, so progressing to final round.')
                log = open(deletionlog, "a+", newline='\n')
                log.write(f"\n{name_variant2} has less than 16 genes remaining, so progressing to final round.")
                log.close()

                continue

            elif remaininggenesn > 16:
                # delete singleko_trigger if exists
                if os.path.isfile(nthsingleko_trigger):
                    os.remove(nthsingleko_trigger)

                ### use roundnumbername > as future
                roundnumbername = roundnumber_names[x]
                #### use name > as past / current results
                eightPanelGroupingsGeneration(remaininggenes_list, deletedgenes_list, name, roundnumbername)


                ### use roundnumbername > as future
                # !!!! gets input from {}eightsandgenes.txt / {}singlekopowersets.txt OR list returned from prior function?

                print(f'\n{name_variant2} has more than 16 genes remaining, so continuing onto next round.')
                log = open(deletionlog, "a+", newline='\n')
                log.write(f"\n{name_variant2} has more than 16 genes remaining, so continuing onto next round.")
                log.close()

                continue

    # Completion message
    finalyellow = finaloutput.format('yellow')
    finalred = finaloutput.format('red')
    finalblue = finaloutput.format('blue')
    if os.path.isfile(finalyellow) and os.path.isfile(finalred) and os.path.isfile(finalblue):
        print(f'\n\nMinesweeper has finished :) see OUTPUT_final for results.')
        log = open(deletionlog, "a+", newline='\n')
        log.write(f"\n\nMinesweeper has finished :) see OUTPUT_final for results.")
        log.close()

def createKOlist(eightsandgenes_name, combosandgenes_name, roundnumbername):
    '''We use this method to create a pickle file of multiple genes to be
    deleted at once. This pickled file will be read in
    .../wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py'''

    df_names = pd.read_csv("gene_name_id_touse.csv", header=0, sep=';')
    color = combosandgenes_name.split('_')[1]
    to_sim_name = []
    to_sim_ids = []
    c=0
    st = time.time()
    powersetcombos = open(eightsandgenes_name, 'r+').readlines()
    for l1 in powersetcombos[1:]:
        c+=1
        list_kos = l1.split('\t')
        list_kos = ','.join(list_kos[1:])
        l = [str(x.replace('\n', '').replace('\t', '').replace(' ', '').replace(',','')) for x in list_kos.strip().split()]
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

    create_scripts(c, color+'_'+roundnumbername)
    with open("stage4_sims_"+roundnumbername, "wb") as fp:   #Pickling
        pickle.dump(to_sim_ids, fp)
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

def create_scripts(c, combosandgenes_name):
    """ Convert template script to bash script using user input, create exp and ko txt files """
    deletionlog = "OUTPUT_final/deletionlog.txt"

    user_input_values = open(USER_INPUT_TXT, "r+").readlines()[:10]
    value_list = [line.split()[1] for line in user_input_values]
    launchpad_path, project_path, runscripts_path, seeds_number = value_list

    experimentscript = EXPERIMENT_SCRIPT.format(combosandgenes_name)
    exp_script = open(experimentscript, 'w', newline='\n')


    for i in range (0, c, 6):
        first_index=i
        last_index=i+5
        if last_index>=c:
            last_index=c-1

        exp_script.write(f'VARIANT="double_knockout" FIRST_VARIANT_INDEX={first_index} LAST_VARIANT_INDEX={last_index} \\\n')
        exp_script.write(f'CACHED_SIM_DATA=0 PARALLEL_PARCA=1 \\\n')
        exp_script.write(f'SINGLE_DAUGHTERS=0 N_GENS=1 N_INIT_SIMS={seeds_number} \\\n')
        exp_script.write(f'MASS_DISTRIBUTION=1 GROWTH_RATE_NOISE=1 D_PERIOD_DIVISION=1 LOG_TO_DISK_EVERY=60 \\\n')
        exp_script.write(f'LAUNCHPAD_FILE="{launchpad_path}/my_launchpad.yaml" \\\n')
        exp_script.write(f'python {runscripts_path}/fireworks/fw_queue.py\n\n')

    exp_script.write('nohup qlaunch -r -l my_launchpad.yaml -w my_fworker.yaml -q my_qadapter.yaml rapidfire --nlaunches infinite --sleep 40 --maxjobs_queue 100 \n')

    print(f'\nCreated bash script (*.sh), in OUTPUT_script_4.')
    print(f'\n> Now you should copy the bash script to your wcEcoli/runscripts folder.')
    print(f'\n> Remember to read stage4_sims in /wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py.')

    with open(deletionlog, "a+", newline='\n') as log:
        log.write("\n> Created bash script/s in OUTPUT_script_3.\n")
        log.write("\n> Now you should copy the bash script to your wcEcoli/runscripts folder.\n")
        log.write("\n> Remember to read stage3_sims in /wcm2022/wcEcoli/models/ecoli/sim/variants/double_knockout.py.")



def main():
    splashscreen()
    names, roundnumber, previousroundnumber = interpretResults()
    dividingcounter_list = createDividingTxt(names)
    remainingGenes(names, dividingcounter_list, previousroundnumber)
    endingDecision(names, dividingcounter_list, JOB, SIMNAME, previousroundnumber, roundnumber)


main()
