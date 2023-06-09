#  Usage: python3 removeGates.py <input.qasm> <output.qasm> <#remove gates> <seed>

import sys
import os
import re
import random as rd
from random import randrange


rd.seed(int(sys.argv[4]))
nRemove = int(sys.argv[3])
prefixLines = []
gateLines = []
with open(sys.argv[1], 'r') as inFile:
    lines = inFile.readlines()
    for line in lines:
        words = re.split(r' ', line)
        while(' ' in words):
            words.remove(' ')
        if(words[0] == 'OPENQASM'):
            prefixLines.append(line)
        elif(words[0] == 'qreg'):
            prefixLines.append(line)
        elif(words[0] == 'creg'):
            prefixLines.append(line)
        elif(words[0] == 'include'):
            prefixLines.append(line)
        else:
            gateLines.append(line)

    for _ in range(nRemove):
        gateLines.pop(randrange(len(gateLines)))

with open(sys.argv[2], 'w') as outFile:
    for line in prefixLines:
        outFile.write(line)
    for line in gateLines:
        outFile.write(line)
