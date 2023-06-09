#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

#include "qCkt.h"

/**Function*************************************************************

  Synopsis    [parse QASM file and store information into gate/qubit/n]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

Circuit* qasmParser(const std::string &filename)
{
    
    // Parse QASM files
    std::ifstream inFile;

    inFile.open(filename);
    if (!inFile)
    {
        std::cerr << "File \"" << filename << "\" cannot be opened." << std::endl;
        return 0;
    }

	Circuit* circuit = nullptr;
    std::string inStr;
    
	int nQubits = -1;
    while(getline(inFile, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        if (inStr.find_first_not_of("\t\n ") == std::string::npos) continue;

		std::stringstream inStr_ss(inStr);
		getline(inStr_ss, inStr, ' ');
		
		if (inStr == "qreg")
		{
			getline(inStr_ss, inStr, '[');
			getline(inStr_ss, inStr, ']');
			nQubits = stoi(inStr);
			circuit = new Circuit(nQubits);
		}
		else if (inStr == "creg"){;}
		else if (inStr == "OPENQASM"){;}
		else if (inStr == "include"){;}
		else if (inStr == "measure"){;}
		else
		{
			assert(circuit);
			if (inStr == "x")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::X, std::vector<int>(1, iQubit));
			}
			else if (inStr == "y")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::Y, std::vector<int>(1, iQubit));
			}
			else if (inStr == "z")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::Z, std::vector<int>(1, iQubit));
			}
			else if (inStr == "h")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::H, std::vector<int>(1, iQubit));
			}
			else if (inStr == "s")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::S, std::vector<int>(1, iQubit));
			}
			else if (inStr == "sdg")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::SDG, std::vector<int>(1, iQubit));
			}
			else if (inStr == "t")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::T, std::vector<int>(1, iQubit));
			}
			else if (inStr == "tdg")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::TDG, std::vector<int>(1, iQubit));
			}
			else if (inStr == "rx(pi/2)")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::RX_PI_2, std::vector<int>(1, iQubit));
			}
			else if (inStr == "rx(-pi/2)")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::RX_PI_2_DG, std::vector<int>(1, iQubit));
			}
			else if (inStr == "ry(pi/2)")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::RY_PI_2, std::vector<int>(1, iQubit));
			}
			else if (inStr == "ry(-pi/2)")
			{
				getline(inStr_ss, inStr, '[');
				getline(inStr_ss, inStr, ']');
				int iQubit = stoi(inStr);
				assert(iQubit < nQubits);
				circuit->addGate(GateType::RY_PI_2_DG, std::vector<int>(1, iQubit));
			}
			else if (inStr == "cx")
			{
				std::vector<int> iQubitsList(2);
				for (int i = 0; i < 2; i++)
				{
					getline(inStr_ss, inStr, '[');
					getline(inStr_ss, inStr, ']');
					iQubitsList[i] = stoi(inStr);
					assert(iQubitsList[i] < nQubits);
				}
				circuit->addGate(GateType::CX, iQubitsList);
			}
			else if (inStr == "cz")
			{
				std::vector<int> iQubitsList(2);
				for (int i = 0; i < 2; i++)
				{
					getline(inStr_ss, inStr, '[');
					getline(inStr_ss, inStr, ']');
					iQubitsList[i] = stoi(inStr);
					assert(iQubitsList[i] < nQubits);
				}
				circuit->addGate(GateType::CZ, iQubitsList);
			}
			else if (inStr == "swap")
			{
				std::vector<int> iQubitsList(2);
				for (int i = 0; i < 2; i++)
				{
					getline(inStr_ss, inStr, '[');
					getline(inStr_ss, inStr, ']');
					iQubitsList[i] = stoi(inStr);
					assert(iQubitsList[i] < nQubits);
				}
				circuit->addGate(GateType::SWAP, iQubitsList);
			}
			else if (inStr == "cswap")
			{
				std::vector<int> iQubitsList(3);
				for (int i = 0; i < 3; i++)
				{
					getline(inStr_ss, inStr, '[');
					getline(inStr_ss, inStr, ']');
					iQubitsList[i] = stoi(inStr);
					assert(iQubitsList[i] < nQubits);
				}
				circuit->addGate(GateType::CSWAP, iQubitsList);
			}
			else if (inStr == "ccx" || inStr == "mcx")
			{
				std::vector<int> iQubitsList(0);
				getline(inStr_ss, inStr, '[');
				while(getline(inStr_ss, inStr, ']'))
				{
					iQubitsList.push_back(stoi(inStr));
					assert(iQubitsList.back() < nQubits);
					getline(inStr_ss, inStr, '[');
				}
				circuit->addGate(GateType::CCX, iQubitsList);
			}
			else
			{
				std::cerr << std::endl
						<< "[Warning]: Syntax \'" << inStr << "\' is not supported. The line is ignored ..." << std::endl;
			}
        }
    }

    inFile.close();

	return circuit;
}
