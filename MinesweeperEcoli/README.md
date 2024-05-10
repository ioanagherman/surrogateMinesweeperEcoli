# Minesweeper for *E. coli*
This is an adaptation of the Minesweeper algorithm published in "Designing minimal genomes using whole-cell models" for *M. genitalium*. The algorithm has been changed to work together with the *E. coli* whole-cell model and its corresponding machine learning surrogate. The aim of the algorithm is to remove as many genes as possible from the *in-silico E. coli* genome (modelled using the WCM), while still maintaining division. The paper describing the results obtained using this code is "Accelerated design of *Escherichia coli* genomes with reduced size using a whole-cell model and machine learning surrogate".

The algorithm uses Python and it should be run together with the single-cell *E. coli* WCM published by the Covert Lab and Stanford University. The latest version of the WCM can be found [here](https://github.com/CovertLab/wcEcoli). 

## Algorithm description
The algorithm is split into 4 stages. 
![framework.pdf](https://github.com/ioanagherman/MinesweeperEcoliWCM/files/14747257/framework.pdf)

### Stage 1 (script_1.py):
The algorithm starts by preparing the scripts needed to run all the single gene knockouts using the WCM. To achieve this, a list of all the genes (in the same order as they are present in the genome) modelled by the WCM needs to be provided in INPUT_script_1/genes.txt. Furthermore, you need to provide a csv file that links each gene name to its ID used by the WCM. This file can be generated using the simData object from any simulation of the WCM.

script_1.py will ask for the user to input the location where the *E. coli* WCM is installed, the location of the Fireworks launchpad files (see the WCM GitHub page for more details), the location of the run scripts (see the WCM GitHub page for more details) and the number of repetitions that each simulation should be run for. This information will be saved in INPUT_script_1/user_input.txt and it will be used throughout all stages of Minesweeper.

script_1.py will produce the mineinputko.sh file that will be needed to run all the single gene knockouts using the WCM. This file should be placed inside wcEcoli/runscripts/ under the folder used to run Fireworks scripts.

### Stage 2 (script_2.py):
The script takes as input a list of all genes and their corresponding labels (essential/non-essential) as predicted by the WCM (stored in INPUT_script_1/inputko_endtimes.txt), and it creates a txt file containing only the non-essential genes (nonessential.txt). Then it splits this list of non-essential genes in 56 segments of different sizes such that there will be one segment containing all non-essential genes, two segments containing the first 90% and the last 90% non-essential genes, similarly two segments containing 80%, 70%, 60%, 50% of the non-essential genes, 3 segments containing 33% of the non-essential genes, four 25% segments, eight 12.5% segments, ten 10% segments and twenty 5% segments. These will all be saved in individual txt files called divisionsegmentX.txt. The groups of genes from each of these segments will be knocked out using the WCM.

script_2.py also outputs a mineinput.sh file that will contain the commands needed to run the corresponding knock outs using the WCM. Furthermore, the IDs of the genes knocked out for each of the 56 simulations will be stored in a pickle file called stage2_sims. This file will have to be read in wcEcoli/models/ecoli/sim/variants/double_knockout.py.

### Stage 3 (script_3.py):
The input for this script is a text file containing the division label when each of the 56 groups is removed (INPUT_script_3/divideko_endtimes.txt). The script then takes the genes deleted to obtain the smallest 3 division-producing cells from stage 2 (we call them red, yellow and blue) and it creates powersets between each of these and the genes deleted to obtain the remaining division-producing segments from stage 2. The script outputs a color_0_combosandgenes.txt file that contains the genes to be deleted corresponding to each of these powerset combinations. Furthermore, the script takes into consideration that there might be powersets that are run in two or more of the colours (red, yellow, blue) and it doesn't add these to the output files if they were part of a previous colour. However, it saves these powersets in separate files (color_already_ran.txt). 

script_3.py also outputs a mineinputko_color.sh file that will contain the commands needed to run the corresponding knock outs using the WCM (to be put in the wcEcoli/runscripts folder). Furthermore, the IDs of the genes knocked out for each simulation will be stored in a pickle file called stage3_color_sims. This file will have to be read in wcEcoli/models/ecoli/sim/variants/double_knockout.py. Depending on the size of the genome, this script might produce thousands of simulations and therefore the creation of these pickle files will take several minutes.

### Stage 4 (script_4.py):
This stage is cyclical, in the sense that the script will be run several times given different (updated) outputs from the WCM. For each of the colours from stage 3, we take the largest deletion segment that produced dividing cells, split the remaining genes (not part of the deletion segment) into 8 groups, and then create powersets of the largest segments and these groups. If none of these divide, in the next iteration of script_4 we append individual genes to the largest deletion segment and simulate the corresponding cells. If some of the 8 groups' powersets divide, we update the largest deletion segment producing a dividing cell and the remaining genes accordingly. If when appending individual genes, some cells divide, the list of remaining genes is updated. The algorithm stops either when there are less than 16 genes remaining outside the largest deletion segment or when there is no division following the single individual gene addition to the largest segment.

Note: In the "Accelerated design of *Escherichia coli* genomes with reduced size using a whole-cell model and machine learning surrogate" paper, this loop was stopped before any of the conditions above were met because the genome obtained was considered to be small enough.

script_4.py is run several times, for each iteration given the output from the WCM. The script outputs .sh files containing the running instructions for the WCM as well as the pickle files with the IDs of the genes that have to be knocked out for each colour and each iteration (stage4_sims_color_iteration).


