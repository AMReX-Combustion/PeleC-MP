#!/usr/bin/env python

# This converts data from Canteras write_csv function to Fortran input
# Assumes Cantera data is SI so it converts it CGS
# Assumes Cantera outputs the mass fractions, not the mole fractions
import argparse
import numpy as np
import sys
import csv

parser = argparse.ArgumentParser()

parser.add_argument("-f", help="Input csv file name", type=str, default='')

parser.add_argument("-out", help="Output file name", type=str, default='')

args = parser.parse_args()

infile = args.f
outfile = args.out

if (infile == ''):
    inputfile = raw_input("Give input csv file name: ")
    infile = inputfile
if (outfile == ''):
    inputfile = raw_input("Give desired output file name: ")
    outfile = inputfile

numrow = 0
numcol = 0
with open(infile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    numcol = len(next(csv_reader))
    csv_file.seek(0)
    for row in csv_reader:
        numrow += 1
# Find number of species by removing z, u, V, T, and rho
numextra = 5
numspec = numcol - numextra
numrow -= 1 # Remove the header row
# Specify which is the x, u, T, rho, and cn columns
xcol = 0
ucol = 1
tcol = 3
rhocol = 4
cncol = 5
# Extract and correct units (SI to CGS) for x, u, T, rho, and c_n
xdata = np.zeros(numrow)
udata = np.zeros(numrow)
tdata = np.zeros(numrow)
rhodata = np.zeros(numrow)
cndata = np.zeros((numrow, numspec))
convx = 100.
convu = 100.
convrho = 1000./(100**3)
spec_names = []
with open(infile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    k = 0
    for row in csv_reader:
        # Extract names of species
        if (k == 0):
            for n in range(numextra, numspec+numextra):
                spec_names.append(row[n])
        else:
            xdata[k-1] = float(row[xcol])*convx
            udata[k-1] = float(row[ucol])*convu
            tdata[k-1] = float(row[tcol])
            rhodata[k-1] = float(row[rhocol])*convrho
            for n in range(0,numspec):
                curcol = cncol + n
                cndata[k-1,n] = float(row[curcol])
        k += 1
finalx = xdata[numrow-1]
startx = 0.
with open(outfile, mode='w') as csv_file:
    headerstring = 'VARIABLES = \"X\" \"U\" \"T\" \"rho\"'
    for n in range(0,numspec):
        headerstring += ' \"' + spec_names[n] + '\"'
    headerstring += '\n'
    csv_file.write(headerstring)
    csv_writer = csv.writer(csv_file, delimiter=' ',lineterminator='\n')
    k = 0
    for row in range(1,numrow):
        outval = np.array([startx + xdata[k], udata[k], tdata[k], rhodata[k]])
        outval = np.append(outval, cndata[k,:])
        csv_writer.writerow(outval)
        k += 1
