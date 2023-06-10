#include <boost/program_options.hpp>
#include <sys/time.h>
#include <fstream>

#include "qcCheck.h"
#include "memMeasure.h"

extern Circuit *qasmParser(const std::string &filename);

int main(int argc, char **argv)
{
    // Program options
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()("help", "produce help message.")(
        "reorder",
        po::value<bool>()->default_value(true),
        "allow variable reordering or not.\n"
        "0: disable 1: enable default: 1")(
        "bitwidth_control",
        po::value<std::string>()->default_value("extend_bitwidth"),
        "bitwidth control when overflowing\n"
        "extend_bitwidth: extend r when overflows\n"
        "drop_lsb: drop lsb when lsb are all 0.")(
        "init_bitwidth",
        po::value<int>()->default_value(4),
        "initial bitwidth r\n"
        "default: 4")("circuit1",
                      po::value<std::string>()->implicit_value(""),
                      "circuit1 under equivalence checking.\n")(
        "circuit2",
        po::value<std::string>()->implicit_value(""),
        "circuit2 under equivalence checking.\n")(
        "approach",
        po::value<std::string>()->default_value("miter"),
        "approach of equivalence checking\n"
        "miter: construct miter\n"
        "construct: construct fucntionality\n"
        "simulation: simulation\n");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
    {
        std::cout << description;
        return 0;
    }

    bool fReorder = vm["reorder"].as<bool>();
    int fInitBitWidth = vm["init_bitwidth"].as<int>();
    std::string fApproach = vm["approach"].as<std::string>();
    std::string fBitWidthControl = vm["approach"].as<std::string>();

    Circuit *circuitU = qasmParser(vm["circuit1"].as<std::string>());
    Circuit *circuitV = qasmParser(vm["circuit2"].as<std::string>());

    int nQubits = circuitU->getNumberQubits();
    assert(circuitV->getNumberQubits() == nQubits);

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    Checker checker;
    checker.setInitBitWidth(fInitBitWidth);
    checker.setAutoReorder(fReorder);
    if (fBitWidthControl == "extend_bitwidth")
        checker.setBitWidthControl(BitWidthControl::ExtendBitWidth);
    else if (fBitWidthControl == "drop_lsb")
        checker.setBitWidthControl(BitWidthControl::DropLSB);
    else
        assert(0 && "Undefined argument bitwidth_control");

    if (fApproach == "miter")
        checker.checkByConstructMiter(circuitU, circuitV);
    else if (fApproach == "construct")
        checker.checkByConstructFunctionality(circuitU, circuitV);
    else if (fApproach == "simulation")
        checker.checkBySimulation(circuitU, circuitV);
    else
        assert(0 && "Undefined argument approach");

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS();

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
