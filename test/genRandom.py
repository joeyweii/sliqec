# Usage: python3 genRandom.py <#qubits> <#gates> <output_file> <seed>

import sys
import random as rd
import os

nQubits = int(sys.argv[1])
nGates = int(sys.argv[2])

rd.seed(int(sys.argv[4]))

with open(sys.argv[3], "wt") as fout:
    fout.write('OPENQASM 2.0;\ninclude \"qelib1.inc\";\n')
    fout.write('qreg q[' + str(nQubits) + '];\n')

    for j in range(nGates):
        gatetype = rd.randint(0, 4)
        if gatetype == 0: # H
            qubit = rd.randint(0, nQubits-1)
            fout.write('h q[' + str(qubit) +'];\n')
        elif gatetype == 1: # S
            qubit = rd.randint(0, nQubits-1)
            fout.write('s q[' + str(qubit) +'];\n')
        elif gatetype == 2: # T
            qubit = rd.randint(0, nQubits-1)
            fout.write('t q[' + str(qubit) +'];\n')
        elif gatetype == 3: # CX
            qubit1 = rd.randint(0, nQubits-1)
            qubit2 = qubit1
            while qubit2 == qubit1:
                qubit2 = rd.randint(0, nQubits-1)
            fout.write('cx q[' + str(qubit1) +'], q[' + str(qubit2) +'];\n')
        elif gatetype == 4: # CCX
            qubit1 = rd.randint(0, nQubits-1)
            qubit2 = rd.randint(0, nQubits-1)
            while qubit2 == qubit1:
                qubit2 = rd.randint(0, nQubits-1)
            qubit3 = rd.randint(0, nQubits-1)
            while qubit3 == qubit1 or qubit3 == qubit2:
                qubit3 = rd.randint(0, nQubits-1)

            fout.write('ccx q[' + str(qubit1) +'], q[' + str(qubit2) +'], q['+str(qubit3)+'];\n')
        elif gatetype == 5: # Rx
            qubit = rd.randint(0, nQubits-1)
            fout.write('rx(pi/2) q[' + str(qubit) +'];\n')
        elif gatetype == 6: # Ry
            qubit = rd.randint(0, nQubits-1)
            fout.write('ry(pi/2) q[' + str(qubit) +'];\n')
        elif gatetype == 7: # X
            qubit = rd.randint(0, nQubits-1)
            fout.write('x q[' + str(qubit) +'];\n')
        elif gatetype == 8: # Y
            qubit = rd.randint(0, nQubits-1)
            fout.write('y q[' + str(qubit) +'];\n')
        elif gatetype == 9: # Z
            qubit = rd.randint(0, nQubits-1)
            fout.write('z q[' + str(qubit) +'];\n')
