#ifndef _QUANTUM_CIRCUIT_H_
#define _QUANTUM_CIRCUIT_H_

#include <vector>

enum class GateType
{
	X,              // Pauli-X
	Y,              // Pauli-Y
	Z,              // Pauli-Z
	H,              // Hadamard
	S,              // Phase 
	SDG,            // Phase inverse
	T,              // Pi/8
	TDG,            // Pi/8 inverse
	RX_PI_2,        // Rotation-X Pi/2
	RX_PI_2_DG,     // Rotation-X Pi/2 inverse
	RY_PI_2,        // Rotation-Y Pi/2 
	RY_PI_2_DG,     // Rotation-Y Pi/2 inverse
	CX,             // Controlled-NOT
	CZ,             // Controlled-Z
	CCX,            // Toffoli
	SWAP,           // SWAP
	CSWAP           // Fredkin
};

class Gate
{
public:
	explicit Gate(const GateType gateType, const std::vector<int> &qubits)
	:_gateType(gateType), _qubits(qubits)
	{}

	const GateType getType() const { return _gateType; }
	const std::vector<int>& getQubits() const { return _qubits; }

private:
	GateType _gateType;
	std::vector<int> _qubits;
};

class Circuit 
{
public:
	explicit Circuit(const int nQubits)
	: _nQubits(nQubits)
	{}

	const Gate* getGate(const int index) const
	{
        assert(index < static_cast<int>(_vGates.size()));
		return _vGates[index];
	}

	void addGate(const GateType gateType, const std::vector<int> &qubits)
	{
		for(const auto qubit: qubits)
			assert(qubit < _nQubits);

		_vGates.push_back(new Gate(gateType, qubits));
	}

	int getNumberQubits() const { return _nQubits; }
	int getGateCount() const { return _vGates.size(); }

private:
	std::vector<Gate*> _vGates;
	int  _nQubits;
};
#endif
