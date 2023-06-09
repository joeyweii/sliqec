"Replace CCX in random Clifford+T+CCX circuits with Cliffor+T gates"
"Usage: python3 replaceCCX.py <Input QASM> <Output QASM>"

import sys
import os
import re

with open(sys.argv[1], "rt") as fin:
    with open(sys.argv[2], "wt") as fout:
        for line in fin:
            if 'ccx' in line:
                strlist = re.split("\[|\]", line)
                qubit1 = strlist[1]
                qubit2 = strlist[3]
                qubit3 = strlist[5]
                fout.write('h q[' + qubit3 + '];\n')
                fout.write('cx q[' + qubit3 + '], q[' + qubit2 + '];\n')
                fout.write('tdg q[' + qubit2 + '];\n')
                fout.write('cx q[' + qubit1 + '], q[' + qubit2 + '];\n')
                fout.write('t q[' + qubit2 + '];\n')
                fout.write('cx q[' + qubit3 + '], q[' + qubit2 + '];\n')
                fout.write('tdg q[' + qubit2 + '];\n')
                fout.write('cx q[' + qubit1 + '], q[' + qubit2 + '];\n')
                fout.write('t q[' + qubit2 + '];\n')
                fout.write('cx q[' + qubit1 + '], q[' + qubit3 + '];\n')
                fout.write('tdg q[' + qubit3 + '];\n')
                fout.write('cx q[' + qubit1 + '], q[' + qubit3 + '];\n')
                fout.write('t q[' + qubit1 + '];\n')
                fout.write('t q[' + qubit3 + '];\n')
                fout.write('h q[' + qubit3 + '];\n')
            else:
                fout.write(line)
