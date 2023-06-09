#ifndef _QUANTUM_CIRCUIT_H_
#define _QUANTUM_CIRCUIT_H_

#include <vector>

enum class GateType
{
    X,          // Pauli-X
    Y,          // Pauli-Y
    Z,          // Pauli-Z
    H,          // Hadamard
    S,          // Phase
    SDG,        // Phase inverse
    T,          // Pi/8
    TDG,        // Pi/8 inverse
    RX_PI_2,    // Rotation-X Pi/2
    RX_PI_2_DG, // Rotation-X Pi/2 inverse
    RY_PI_2,    // Rotation-Y Pi/2
    RY_PI_2_DG, // Rotation-Y Pi/2 inverse
    CX,         // Controlled-NOT
    CZ,         // Controlled-Z
    CCX,        // Toffoli
    SWAP,       // SWAP
    CSWAP       // Fredkin
};

class Gate
{
public:
    explicit Gate(const GateType gateType, const std::vector<int> &qubits)
        : _gateType(gateType)
        , _qubits(qubits)
    {
    }

    const GateType getType() const { return _gateType; }
    const std::vector<int> &getQubits() const { return _qubits; }
    void setType(const GateType gateType) { _gateType = gateType; }
    void liftQubits(int numLift)
    {
        for (int &qubit : _qubits)
            qubit += numLift;
    }

    void unliftQubits(int numLift)
    {
        for (int &qubit : _qubits)
            qubit -= numLift;
    }

private:
    GateType _gateType;
    std::vector<int> _qubits;
};

class Circuit
{
public:
    explicit Circuit(const int nQubits)
        : _nQubits(nQubits)
    {
    }

    const Gate *getGate(const int index) const { return _vGates[index]; }

    void addGate(const GateType gateType, const std::vector<int> &qubits)
    {
        _vGates.push_back(new Gate(gateType, qubits));
    }

    int getNumberQubits() const { return _nQubits; }
    int getGateCount() const { return _vGates.size(); }

    void liftAllGateQubits(int numLift)
    {
        for (Gate *gate : _vGates)
            gate->liftQubits(numLift);
    }

    void unliftAllGateQubits(int numLift)
    {
        for (Gate *gate : _vGates)
            gate->unliftQubits(numLift);
    }

    void daggerAllGate()
    {
        for (Gate *gate : _vGates)
        {
            const GateType &gateType = gate->getType();
            if (gateType == GateType::S)
                gate->setType(GateType::SDG);
            else if (gateType == GateType::SDG)
                gate->setType(GateType::S);
            else if (gateType == GateType::T)
                gate->setType(GateType::TDG);
            else if (gateType == GateType::TDG)
                gate->setType(GateType::T);
            else if (gateType == GateType::RX_PI_2)
                gate->setType(GateType::RX_PI_2_DG);
            else if (gateType == GateType::RX_PI_2_DG)
                gate->setType(GateType::RX_PI_2);
        }
    }

private:
    std::vector<Gate *> _vGates;
    int _nQubits;
};
#endif
