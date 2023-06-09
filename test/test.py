import sys
import subprocess
import json
from tqdm import tqdm

if __name__ == '__main__':
    nQubit = 5
    nGate = 100

    #  print("Start checking the equivalent cases.")
    for seed in tqdm(range(100), desc='Checking the equivalent cases...'):
        genUCommand = "python3 ./genRandom.py {} {} U.qasm {}".format(nQubit, nGate, seed)
        subprocess.run(genUCommand, shell=True)
        genVCommand = "python3 ./replaceCCX.py U.qasm V.qasm"
        subprocess.run(genVCommand, shell=True)

        checkCommand = "../SliQEC_v2 --circuit1 U.qasm --circuit2 V.qasm"
        checkOutput = subprocess.getoutput(checkCommand)
        checkOutputJSON = json.loads(checkOutput)
        if(checkOutputJSON["equivalence"] == "not_equivalent"):
            print("Test failed.")
            exit()

        subprocess.run("rm U.qasm V.qasm", shell=True)

    for seed in tqdm(range(100), desc='Checking the non-equivalent cases...'):
        genUCommand = "python3 ./genRandom.py {} {} U.qasm {}".format(nQubit, nGate, seed)
        subprocess.run(genUCommand, shell=True)
        genVtemCommand = "python3 ./replaceCCX.py U.qasm V_tem.qasm"
        subprocess.run(genVtemCommand, shell=True)
        genVCommand = "python3 ./removeGates.py V_tem.qasm V.qasm 5 {}".format(seed)
        subprocess.run(genVCommand, shell=True)

        checkCommand = "../SliQEC_v2 --circuit1 U.qasm --circuit2 V.qasm"
        checkOutput = subprocess.getoutput(checkCommand)
        checkOutputJSON = json.loads(checkOutput)
        if(checkOutputJSON["equivalence"] == "equivalent"):
            print("Test failed.")
            exit()

        subprocess.run("rm U.qasm V.qasm V_tem.qasm", shell=True)

    print("Test passed.")
