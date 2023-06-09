#include <boost/program_options.hpp>
#include <sys/time.h> 
#include <fstream>

#include "qcCheck.h"
#include "memMeasure.h"

extern Circuit* qasmParser(const std::string& filename);

int main(int argc, char **argv)
{
    // Program options
    namespace po = boost::program_options;
    po::options_description description("Options");
    description.add_options()
    ("help", "produce help message.")
    ("reorder", po::value<bool>()->default_value(true), "allow variable reordering or not.\n"
                                                             "0: disable 1: enable default: 1") 
	("bitwidth_control", po::value<int>()->default_value(0), "bitwidth control when overflowing\n"
												 "0: extendBitwidth 1: dropLSB")
	("init_bitwidth", po::value<int>()->default_value(32), "initial bitwidth r\n"
															"default: 32")
    ("circuit1", po::value<std::string>()->implicit_value(""), "circuit1 under equivalence checking")
    ("circuit2", po::value<std::string>()->implicit_value(""), "circuit2 under equivalence checking.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, description), vm);
    po::notify(vm);

    if (vm.count("help") || argc == 1)
    {
	    std::cout << description;
	    return 0;
	}

    bool fReorder = vm["reorder"].as<bool>();
	int	 fBitWidthControl = vm["bitwidth_control"].as<int>();
	int	 fInitBitWidth = vm["init_bitwidth"].as<int>();

    Circuit *circuitU = qasmParser(vm["circuit1"].as<std::string>());
    Circuit *circuitV = qasmParser(vm["circuit2"].as<std::string>());

	int nQubits = circuitU->getNumberQubits();
    assert(circuitV->getNumberQubits() == nQubits);

    struct timeval tStart, tFinish;
    double elapsedTime;
    double runtime;
    size_t memPeak;

    gettimeofday(&tStart, NULL);

    Checker checker(nQubits, fInitBitWidth, fBitWidthControl, fReorder);

    checker.check(circuitU, circuitV);

    gettimeofday(&tFinish, NULL);
    elapsedTime = (tFinish.tv_sec - tStart.tv_sec) * 1000.0;
    elapsedTime += (tFinish.tv_usec - tStart.tv_usec) / 1000.0;

    runtime = elapsedTime / 1000.0;
    memPeak = getPeakRSS();
    
    std::cout << "----- Circuit Info. -----\n";
    std::cout << "#Qubits: " << nQubits << '\n';
    std::cout << "#Gates in circuit1: " << circuitU->getGateCount() << '\n';
    std::cout << "#Gates in circuit2: " << circuitV->getGateCount() << '\n';

    std::cout << "----- Resource Usage -----\n";
    std::cout << "Runtime: " << runtime << " seconds\n";
    std::cout << "Peak memory usage: " << memPeak << " bytes\n"; 
    std::cout << "--------------------------\n";

    return 0;
}
