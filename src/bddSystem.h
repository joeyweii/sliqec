#ifndef _BDDSYSTEN_H_
#define _BDDSYSTEM_H_

#include <iostream>
#include <cstdlib> 
#include <string> 
#include <vector>
#include <deque>

#include "../cudd/cudd/cudd.h"
#include "../cudd/cudd/cuddInt.h"
#include "../cudd/util/util.h"
#include "qCkt.h"

class Tensor
{
public:
	int		_k;
	int		_r;
	int		_rank;
	std::vector<std::deque<DdNode*>> _allBDD;

	Tensor(int r, int rank)
	:_k(0), _r(r), _rank(rank),
	_allBDD(std::vector<std::deque<DdNode*>>(4, std::deque<DdNode*>(_r, nullptr)))
	{}
};

class BDDSystem
{
protected:

	enum class BitWidthMode
	{
		ExtendBitWidth,
		DropLSB,
		DropMSB
	};

    explicit BDDSystem(int maxRank, int fBitWidthMode, bool fReorder)
    :   _ddManager(nullptr),
        _w(4), _maxNodeCount(0)
    {
		_ddManager = Cudd_Init(maxRank, maxRank, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
		switch(fBitWidthMode)
		{
			case 1: _bitWidthMode = BitWidthMode::DropLSB; break;
			case 2: _bitWidthMode = BitWidthMode::DropMSB; break;
			default: _bitWidthMode = BitWidthMode::ExtendBitWidth; 
		}
		if (fReorder) Cudd_AutodynEnable(_ddManager, CUDD_REORDER_SYMM_SIFT);
	}

    virtual ~BDDSystem()  
    {
		if(Cudd_CheckZeroRef(_ddManager) > 10)
			std::cout << "\nThe number of referenced nodes = " << Cudd_CheckZeroRef(_ddManager) << std::endl;
		Cudd_Quit(_ddManager);
    }

    /* gateOpe.cpp */
    void Toffoli(Tensor *tensor, const std::vector<int> &qubits);
    void Fredkin(Tensor *tensor, const std::vector<int> &qubits);
    void Hadamard(Tensor *tensor, const std::vector<int> &qubits);
    void Rx_pi_2(Tensor *tensor, const std::vector<int> &qubits);
    void Rx_pi_2_dagger(Tensor *tensor, const std::vector<int> &qubits);
    void Ry_pi_2(Tensor *tensor, const std::vector<int> &qubits, const bool fTranspose);
    void Ry_pi_2_dagger(Tensor *tensor, const std::vector<int> &qubits, const bool fTranspose);
    void Phase_shift(Tensor *tensor, const std::vector<int> &qubits, const int phase); // phase can only be 2 to the power of an integer
    void Phase_shift_dagger(Tensor *tensor, const std::vector<int> &qubits, const int phase);
    void PauliX(Tensor *tensor, const std::vector<int> &qubits);
    void PauliY(Tensor *tensor, const std::vector<int> &qubits, const bool fTranspose);
    void PauliZ(Tensor *tensor, const std::vector<int> &qubits); 
	void applyGate(const Gate* gate, Tensor *tensor, const bool fTranspose);

	/* propCheck.cpp */
	DdNode* sparsityDD(Tensor *tensor) const;
	double sparsity(Tensor *tensor) const;
	bool eqCheckTwoTensor(Tensor *tensor1, Tensor *tensor2);

    /* misc.cpp */
	Tensor* newTensor(int r, int rank);
	void deleteTensor(Tensor* tensor);
	void printTensor(Tensor* tensor) const;
	bool isTensorLSBZero(Tensor* tensor) const;
    void incBDDsBitWidth(std::vector<std::deque<DdNode*>> &allBDD);
    void dropTensorLSB(Tensor* tensor);
    void dropTensorMSB(Tensor* tensor);
	void dropTensorBits(Tensor *tensor);
	void increaseTensorKByOne(Tensor *tensor);
    bool checkAdderOverflow(DdNode *g, DdNode *h, DdNode *crin) const;
    void updateMaxNodeCount();

    DdManager *_ddManager;			// BDD manager.
    int _w;							// # of integers = 4.
    unsigned long _maxNodeCount;	// node count.
	BitWidthMode	_bitWidthMode;			// mode of bits' control
};

#endif
