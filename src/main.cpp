#include <string>
#include <sys/time.h>
#include <fstream>

#include "qcCheck.h"
#include "memMeasure.h"

int main(int argc, char *argv[]) {
    // Program options
    bool fReorder = true;
    int fInitBitWidth = 4;
    std::string fApproach = "miter";
    std::string fBitWidthControl = "extend_bitwidth";
    std::string circuit1Name, circuit2Name;

    for (int i = 1; i < argc;) {
        if (std::string(argv[i]) == "--circuit1") {
            circuit1Name = std::string(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]) == "--circuit2") {
            circuit2Name = std::string(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]) == "--reorder") {
            if (std::string(argv[i + 1]) == "0")
                fReorder = false;
            else if (std::string(argv[i + 1]) == "1")
                fReorder = true;
            else
                assert(0 && "Parse -reorder option fail");
            i += 2;
        } else if (std::string(argv[i]) == "--init_bitwidth") {
            fInitBitWidth = std::stoi(std::string(argv[i + 1]));
            i += 2;
        } else if (std::string(argv[i]) == "--approach") {
            fApproach = std::string(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]) == "--bitwidth_control") {
            fBitWidthControl = std::string(argv[i + 1]);
            i += 2;
        } else
            assert(0 && "Undefined options.");
    }

    assert(circuit1Name != "" && "Circuit1 not specified.");
    assert(circuit2Name != "" && "Circuit2 not specified.");
    Circuit *circuitU = parseQASM(circuit1Name);
    Circuit *circuitV = parseQASM(circuit2Name);

    int nQubits = circuitU->getNumberQubits();
    assert(circuitV->getNumberQubits() == nQubits);

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    double memPeak;

    gettimeofday(&tStart, NULL);

    Checker checker(nQubits);
    checker.setInitBitWidth(fInitBitWidth);
    checker.setAutoReorder(fReorder);
    if (fBitWidthControl == "extend_bitwidth")
        checker.setBitWidthControl(BitWidthControl::ExtendBitWidth);
    else if (fBitWidthControl == "drop_lsb")
        checker.setBitWidthControl(BitWidthControl::DropLSB);
    else
        assert(0 && "Undefined argument --bitwidth_control");

    if (fApproach == "miter")
        checker.checkByConstructMiter(circuitU, circuitV);
    else if (fApproach == "construct")
        checker.checkByConstructFunctionality(circuitU, circuitV);
    else if (fApproach == "simulation")
        checker.checkBySimulation(circuitU, circuitV);
    else
        assert(0 && "Undefined argument --approach");

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS() / 1024.0 / 1024.0 / 1024.0;

    checker.addElementToOutputJSON("num_qubits", std::to_string(nQubits));
    checker.addElementToOutputJSON("num_gates_circuit1",
                                   std::to_string(circuitU->getGateCount()));
    checker.addElementToOutputJSON("num_gates_circuit2",
                                   std::to_string(circuitV->getGateCount()));
    checker.addElementToOutputJSON("runtime", std::to_string(runtime));
    checker.addElementToOutputJSON("memory", std::to_string(memPeak));
    checker.printOutputJSON();

    return 0;
}
