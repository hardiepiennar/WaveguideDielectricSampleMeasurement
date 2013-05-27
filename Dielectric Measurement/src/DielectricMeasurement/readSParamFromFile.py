'''
Created on May 12, 2013

Read S-parameter files into S11,S21,S12,S22 objects
returned into as a S array [S11,S21,S12,S22]
TODO:
Extend functionality to detect parameter types
@author: Hardie Pienaar
'''
from numpy import *
import csv

def readSParamFromFile(filename):
    startRow = 3;
    
    with open(filename) as csvfile:
        reader = csv.reader(csvfile,delimiter=' ')
        #Initiate S-Parameter objects
        freq = []
        S11 = []
        S21 = []
        S12 = []
        S22 = []
        
        rowNo = 0    
        for row in reader:
            #skip header rows        
            if rowNo < startRow:
                rowNo = rowNo + 1
            else:
                #Store objects
                string = '    '.join(row)
                elements = string.split()
                freq.append(float(elements[0]))
                S11.append(float(elements[1]) +1j*float(elements[2]))
                S21.append(float(elements[3]) +1j*float(elements[4]))
                S12.append(float(elements[5]) +1j*float(elements[6]))
                S22.append(float(elements[7]) +1j*float(elements[8]))
        freq = array(freq)
        S11 = array(S11)
        S11 = array(S11)
        S11 = array(S11)
        S11 = array(S11)
    return [freq, S11,S21,S12,S22]

        