"""
Name: Azael Zamora

The PSSM matrix should be somewhat similar to the PSSM matrix in the NCBI PSSM Viewer
Implmentation of the PSSM algorithm, also know as PWM
"""

import math                                     # need math to use log
from Bio.SubsMat import MatrixInfo as matlist   # needed to utilize the Blosum62 for pseudocounts

import pandas as panda                          # pandas is utilized to format the PSSM matrix in a file
from tabulate import tabulate                   # tabulate is utilized to write the DataFrame into a grid format

file_name = input("Enter input file name: ")    # input to type in the file name that you want to use

# read_sequences reads the input file and returns the sequences in a 2d list
def read_sequences(fileName):
    sequences = []                      # the list that will store all the sequences from the input fileName
    with open(fileName, 'r') as myfile: # Opens the fileName as readable
        for line in myfile:             # for each line in the input fileName
            sequences.append(line[0:len(line) - 1]) # appends each sequence as a row

    return sequences                    # returns the 2d list that contains the sequences

sequences_list = read_sequences(file_name)   # sequences_list will contain the sequences from the fileName

min_rows = 100000000000                 # min_rows starts with a huge number at first
num_of_sequences = len(sequences_list)  # num_of_sequences is the number of seqeuences that were included in the input file
blosum = matlist.blosum62               # blosum will contain the blosum 62 substitution matrix

for sequence in sequences_list:         # for each sequence in sequences_list
    if min_rows > len(sequence):        # if min_rows is greater than the len(sequence)
        min_rows = len(sequence)        # min_rows = len(sequence), since PSSM needs the matrix to be in a suitable size to prevent from any errors

"""
The first step of the PSSM is to count the number of times each protein appears in
the sequences at each position

The consensus sequence will be created based on the most frequent amino acid in each position from the sequences
"""

consensus = ''                          # the consensus will contain a sequence that is made up of all the most frequent amino acids in per position
amino_acids = 'AGILVMFWPCSTYNQHKRDE'    # this string essentially contains the 20 amino acids
PFM = []                                # PFM is the Position-Frequency Matrix

for i in range(min_rows):               # from position 0 to min_rows for all the sequences in sequences_list
    tmp = []                            # tempory list that contain the occurence of each amino acid in the order of amino_acids
    index = 0
    for j in range(len(amino_acids)):   # this for loop goes through each amino acid in amino_acids
        amino_count = 0                 # counter for each amino acid in each position

        for k in range(len(sequences_list)):            # inner most for loop iterates through each sequence
            if amino_acids[j] == sequences_list[k][i]:  # if the amino_acid[j] == to the amino acid in sequence #k in position i
                amino_count += 1                        # amino_count is incremented

        tmp.append(amino_count)         # tmp adds the amino_count for the amino_acid[j]

        if tmp[index] < tmp[len(tmp) - 1]:
            index = j                   # index will hold the position of the most frequent amino acid from amino_acids

    consensus += amino_acids[index]     # consensus at position i will have the most frequently occuring Amino Acid that occurs
    PFM.append(tmp)                     # PFM adds the tmp as a row with a column of amino_acid frequency occurences

"""
After the PFM is completed, it should result in a min_rows x 20 matrix in which,
the each row is a position from the sequences and,
each column has the number of occurences for each amino acid

If there are no amino acids occurence in a position, then a pseudocount will be used.
The pseudocount will be based from the Blosum 62 substitution matrix

The pseudocount will be s(a, b), where:
    a is the amino acid from the consensus sequence at the current position
    b is the amino from the current column at the current position
"""

PPM = []                            # this is the position-probability matrix

for a in range(len(PFM)):           # number of rows in the PFM matrix
    tmp = []                        # tmp list that will contain the probability of each amino acid for each position

    for b in range(len(PFM[a])):    # column (each column position is the occurence of an amino acid)

        if PFM[a][b] == 0:          # if the position frequency of the amino acid at PFM[a][b] is 0, then include the pseudocount
            pair = (consensus[a], amino_acids[b])   # pair = s(a, b)

            if pair not in blosum:  # is needed since blosum is lower triangular matrix
                denominator = num_of_sequences + (20*blosum[ (tuple(reversed(pair))) ]) # denominator is the number of sequences from the input file + 20 * the pseudocount
                tmp.append( blosum[ (tuple(reversed(pair))) ] / denominator)            # tmp will append the result of the pseudocount from blosum62 matrix divided by the denominator
            else:
                tmp.append( blosum[pair] / (num_of_sequences + 20*blosum[pair]) )       # same thing as in line 88
        else:
            tmp.append(PFM[a][b] / (num_of_sequences + 20) )                            # divides the occurence of PFM[a][b] / 30

    PPM.append(tmp)                                                                     # append the amino acid probability list to the PPM

"""
The next step is to use determine the probability of the amino acid appearing in the Position
ideally, if a certain cell from PPM cotains a value of 0 it will be skipped to prevent any errors from occuring

if there is a value in a cell, then probability will be calculated as:
    log(PPM[a][b] / (the probability that the amino acid will appear, for the purpose of this, assume each amino acid has equal probability which is 0.05))

"""

for a in range(len(PPM)):           # row (the # of positions for the sequence length)
    for b in range(len(PPM[a])):    # column (the 20 amino acids that can occur in each position)
        if PPM[a][b] != 0:          # if the PPM[a][b] does not have the have the value 0 then continue
            PPM[a][b] = math.log10(abs(PPM[a][b]) / 0.05) # use log(PPM[a][b] / 0.05) to find the probability of the amino acid at the ath position

"""
The final step of implementing a similar PSSM from a PSI-Blast:
    add the PPM value to the blosum 62 substitution matrix s(a, b) where:

        a is the amino acid from the current position of the consensus sequence
        b is the amino acid from current column position of the PPM amino acid columns

doing this will yield a similar PSSM matrix to the one that can obtained from a PSI-Blast and viewed in PSSM viewer
"""
PSSM = PPM                          # this will be the PSSM matrix
for a in range(len(PFM)):           # row (the ath position of the sequence)
    for b in range(len(PFM[a])):    # column (the bth column of the amino acid)

        pair = (consensus[a], amino_acids[b])   # pair = s(a, b) from the blosum 62 substitution matrix
        if pair not in blosum:                  # since blosum 62 is lower triangular matrix
            PSSM[a][b] = PPM[a][b] + blosum[(tuple(reversed(pair)))]    # add the PPM[a][b] value to  the substitution score of s(b, a)
        else:
            PSSM[a][b] = PPM[a][b] + blosum[pair]                       # add the PPM[a][b] value to  the substitution score of s(a, b)

        PSSM[a][b] = round(PSSM[a][b])                                  # rounds the amino acid frequency

"""
Once each PSSM value is rounded, it should be a similar pssm matrix from NCBI PSSM viewer but some of the results will vary
Final step of this program is to save the output in a file using pandas package
"""
output_dict = {}                # output_dict will have the columns stored by amino acid as a key value

rows = []                       # this will contain the number of positions 1-nth position

for i in range(1, min_rows + 1):# the for loop iterates from 1 to min_rows + 1
    rows.append(i)

output_dict['Position'] = rows  # output_dict will now contain the header 'Position' with the column indexes from 1 to nth position

for y in range(len(PSSM[0])):   # goes through the 20 amino acids of a position
    amino_output = []           # temporary holder of the yth amino acid
    for z in range(len(PSSM)):              # goes through the positions of the PSSM Matrix
            amino_output.append(PSSM[z][y]) # appends the yth amino acid from the zth position to amino_output

    output_dict[amino_acids[y]] = amino_output  # now contains the header 'yth amino acid' with column values from the zth position

output_file = panda.DataFrame(output_dict)      # output file is now formated in a matrix format by using pandas.DataFrame

with open(file_name[0 : len(file_name) - 4] + "_output.txt", 'w') as myfile:    # saves the output file as file_name_output.txt
    myfile.write(tabulate(output_file, output_dict.keys(), showindex = "never", tablefmt = 'fancy_grid')) # this formats the output into a grid format using output_dict

# end of the program
