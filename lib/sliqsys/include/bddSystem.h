#pragma once

#include "qCkt.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <deque>

#include "../cudd/cudd/cudd.h"
#include "../cudd/cudd/cuddInt.h"
#include "../cudd/util/util.h"

class Tensor
{
public:
    Tensor(int r, int rank)
        : _k(0)
        , _r(r)
        , _rank(rank)
        , _allBDD(std::vector<std::deque<DdNode *>>(
              4, std::deque<DdNode *>(_r, nullptr)))
    {
    }

    int _k;
    int _r;
    int _rank;
    std::vector<std::deque<DdNode *>> _allBDD;
};

class BDDSystem
{
public:
    enum class BitWidthControl
    {
        ExtendBitWidth,
        DropLSB
    };

    explicit BDDSystem()
        : _ddManager(Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0))
        , _w(4)
        , _initBitWidth(4)
        , _bitWidthControl(BitWidthControl::ExtendBitWidth)
    {
    }

    void setInitBitWidth(const int initBitWidth)
    {
        _initBitWidth = initBitWidth;
    }

    void setAutoReorder(const bool isReorder)
    {
        if (isReorder)
            Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);
        else
            Cudd_AutodynDisable(_ddManager);
    }

    void setBitWidthControl(BitWidthControl bitWidthControl)
    {
        _bitWidthControl = bitWidthControl;
    }

protected:
    /* gateOpe.cpp */
    void Toffoli(Tensor *tensor, const std::vector<int> &qubits);
    void Fredkin(Tensor *tensor, const std::vector<int> &qubits);
    void Hadamard(Tensor *tensor, const std::vector<int> &qubits);
    void Rx_pi_2(Tensor *tensor, const std::vector<int> &qubits);
    void Rx_pi_2_dagger(Tensor *tensor, const std::vector<int> &qubits);
    void Ry_pi_2(Tensor *tensor,
                 const std::vector<int> &qubits,
                 const bool fTranspose);
    void Ry_pi_2_dagger(Tensor *tensor,
                        const std::vector<int> &qubits,
                        const bool fTranspose);
    void Phase_shift(
        Tensor *tensor,
        const std::vector<int> &qubits,
        const int phase); // phase can only be 2 to the power of an integer
    void Phase_shift_dagger(Tensor *tensor,
                            const std::vector<int> &qubits,
                            const int phase);
    void PauliX(Tensor *tensor, const std::vector<int> &qubits);
    void PauliY(Tensor *tensor,
                const std::vector<int> &qubits,
                const bool fTranspose);
    void PauliZ(Tensor *tensor, const std::vector<int> &qubits);
    void applyGate(const Gate *gate, Tensor *tensor, const bool fTranspose);

    /* propCheck.cpp */
    DdNode *sparsityDD(Tensor *tensor) const;
    double sparsity(Tensor *tensor) const;
    bool eqCheckTwoTensor(Tensor *tensor1, Tensor *tensor2);

    /* misc.cpp */
    Tensor *newTensor(int rank);
    void deleteTensor(Tensor *tensor);
    void printTensor(Tensor *tensor) const;
    bool isTensorLSBZero(Tensor *tensor) const;
    void incBDDsBitWidth(std::vector<std::deque<DdNode *>> &allBDD);
    void dropTensorLSB(Tensor *tensor);
    void dropTensorMSB(Tensor *tensor);
    void dropTensorBits(Tensor *tensor);
    void increaseTensorKByOne(Tensor *tensor);
    bool checkAdderOverflow(DdNode *g, DdNode *h, DdNode *crin) const;

    DdManager *_ddManager;            // BDD manager.
    int _w;                           // # of integers = 4.
    int _initBitWidth;                // initial bitwidth when new a tensor
    BitWidthControl _bitWidthControl; // mode of bits' control
};

using BitWidthControl = BDDSystem::BitWidthControl;
