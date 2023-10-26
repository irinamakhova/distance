#!/usr/bin/env python
# coding: utf-8

# This script was created to calculate the distances between 
# the specific atoms using coordinates in a PDB file.

# More specifically, I was interested in the elements of the catalytic triad
# of serine hydrolases (aspartate/glutamate, histidine and serine): 
# Oxygen from the hydroxyl group in serine, which forms a hydrogen bond 
# with the ND1 of histidine, and oxygen from the hydroxyl group in aspartate/glutamate, 
# which forms a hydrogen bond with the hydrogen bonded to NE of histidine.

# I wanted to compare the distances for the protein I was interested in, 
# in alpha fold model, in homology model and in crystal structure. 
# Maybe there is a finer method, but since I was interested in the specific atoms 
# and not the residues in general, it was easier for me to write this than to look 
# for existing possibilities. It was important to me that with the known sequence position, 
# I would get the atom positions and atom types suggested 
# so I can choose without having to search in the PDB file myself.

# This script works with PDB files with 11 or 12 columns.

# Usage Example: python3 distance.py -residueposition1 100 -residueposition2 250 -pdb protein.pdb

# # Distance Calculator

import os
import csv
import argparse
from math import sqrt
from tabulate import tabulate

### 1. Defining the working directory containing this script and input file:

working_directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(working_directory)

### 2. Defining Arguments

parser = argparse.ArgumentParser()

### 2.1. Parameter

parser.add_argument("-residueposition1", help="""E.g. "95" - Distance from which molecule need to be calculated""", type=int)
parser.add_argument("-residueposition2", help="""E.g. "95" - Distance to which molecule need to be calculated""", type=int)

### 2.2. Input Data

parser.add_argument("-pdb", help="""PDB-File of the protein""", type=str)

args = parser.parse_args()

### 3. Defining Lists and Variables

residueposition1 = args.residueposition1
residueposition2 = args.residueposition2

pdb = args.pdb

atompositions1 = []
atoms1 = []
atompositions2 = []
atoms2 = []
molecule1 = ""
molecule2 = ""

# The lists "atomposition1" and "atomposition2" will contain all atom positions belonging to chosen residues (3451-3460);
# The lists "atoms1" and "atoms2" will be filled with corresponding atom types (CA, HA, CB, HB etc.);
# The strings "molecule1" and "molecule2" will represent the three-letter code of the respective amino acid 
# from which or to which the distance will be calculated. This serves to check whether the correct residue was selected 
# with the specification of the position in sequence. 


### 4. Reading the file and determination of atom positions

# Asking for the user input from which and to which atom of the specified residue the distance should be calculated

## Some PDB files contain 11 columns, some consist of 12 columns and the values we are interested in are thus "slipped":
##         ANo  Type 3LC    ResNo        X       Y       Z
## ATOM    748  CA   SER    51      -1.164 -14.104  31.891  1.00  0.00           C   
## ATOM    635  CA   PHE A  87       2.592  -5.734  14.792  1.00 98.81           C 

with open(str(pdb), "r+") as pdb_file:
    pdb_all_lines = pdb_file.readlines()
    for line in pdb_all_lines:
        values = line.split()
        if len(values) == 11:
            if values[4] == str(residueposition1):
                molecule1 = str(values[3]) + str(values[4])
                atompositions1.append(int(values[1]))
                atoms1.append(values[2])
            elif values[4] == str(residueposition2):
                molecule2 = str(values[3]) + str(values[4])
                atompositions2.append(int(values[1]))
                atoms2.append(values[2])
        elif len(values) == 12:
            if values[5] == str(residueposition1):
                molecule1 = str(values[3]) + str(values[5])
                atompositions1.append(int(values[1]))
                atoms1.append(values[2])
            elif values[5] == str(residueposition2):
                molecule2 = str(values[3]) + str(values[5])
                atompositions2.append(int(values[1]))
                atoms2.append(values[2])

atoms_in_residue1 = {"Atom Number": atompositions1, "Atom Type": atoms1}
atoms_in_residue2 = {"Atom Number": atompositions2, "Atom Type": atoms2}

print(tabulate(atoms_in_residue1, headers = "keys",  tablefmt = "mixed_grid"))

atomnumber1 = input(f'From which atom (number) of {molecule1} should the distance be calculated? ')

print(tabulate(atoms_in_residue2, headers = "keys", tablefmt = "mixed_grid"))

atomnumber2 = input(f'To which atom (number) of {molecule2} should the distance be calculated? ')

### 5. Getting x, y, z - Coordinates from PDB file

# If the PDB file contain multiple frames, there will be multiple coordinates 
# belonging to the same atomic position. For this case, lists of coordinates are created.
# I used only one-frame PDB, so it still need to be tested.

X_one = []
X_two = []
Y_one = []
Y_two = []
Z_one = []
Z_two = []

distances = []

### Some PDB files contain 11 columns, some consist of 12 columns and the values we are interested in are thus "slipped":
##         ANo  Type 3LC    ResNo        X       Y       Z
## ATOM    748  CA   SER    51      -1.164 -14.104  31.891  1.00  0.00           C   
## ATOM    635  CA   PHE A  87       2.592  -5.734  14.792  1.00 98.81           C 


with open(str(pdb), "r+") as pdb_file:
    pdb_all_lines = pdb_file.readlines()
    for line in pdb_all_lines:
        values = line.split()
        if len(values) == 11: 
            if int(values[1]) == int(atomnumber1):
                X_one.append(float(values[5]))
                Y_one.append(float(values[6]))
                Z_one.append(float(values[7]))
            elif int(values[1]) == int(atomnumber2):
                X_two.append(float(values[5]))
                Y_two.append(float(values[6]))
                Z_two.append(float(values[7]))
        elif len(values) == 12: 
            if int(values[1]) == int(atomnumber1):
                X_one.append(float(values[6]))
                Y_one.append(float(values[7]))
                Z_one.append(float(values[8]))
            elif int(values[1]) == int(atomnumber2):
                X_two.append(float(values[6]))
                Y_two.append(float(values[7]))
                Z_two.append(float(values[8]))
            else: 
                pass

### 6. Distance calculation

for i in range(0, len(X_one)):
    distance = sqrt((X_one[i] - X_two[i])**2 + (Y_one[i] - Y_two[i])**2 + (Z_one[i] - Z_two[i])**2)
    print(distance)
    distances.append(distance)

### 7. Output in terminal or output file generation

if len(distances) > 1:
    with open(f'{working_directory}/distances_{atomnumber1}_to_{atomnumber2}.txt', 'x') as output_file: 
        for dist in distances:
            output_file.write("%s\n" % dist)
    print("\nOutput files generated: \n" + f'distances_{atomnumber1}_to_{atomnumber2}.txt')
else:
    print(distances)
