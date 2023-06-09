#include "qcCheck.h"

// Constructor
Checker::Checker
(
    int nQubits,
	int fInitBitWidth,
	int fBitWidthControl,
    bool fReorder
)
:   
	BDDSystem
    ( 
		nQubits,
		fBitWidthControl,
        fReorder
    )
{
	_U = newTensor(fInitBitWidth, nQubits*2);
	_V = newTensor(fInitBitWidth, nQubits*2);
}

// Destructor
Checker::~Checker()
{
	deleteTensor(_U);
	deleteTensor(_V);
}

void Checker::check(const Circuit *circuitU, const Circuit * circuitV)
{
    constructUandV(circuitU, circuitV);
}

/**Function*************************************************************

  Synopsis    [Construct the unitary matrix of U and V.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Checker::constructUandV(const Circuit *circuitU, const Circuit *circuitV)
{
    bool fTranspose = false;

    initTensorToIdentityMatrix(_U);
	for(int i = 0; i < circuitU->getGateCount(); ++i)
		applyGate(circuitU->getGate(i), _U, fTranspose);

    initTensorToIdentityMatrix(_V);
	for(int i = 0; i < circuitV->getGateCount(); ++i)
		applyGate(circuitV->getGate(i), _V, fTranspose);

    bool checkResult = eqCheckTwoTensor(_U, _V);
    if(checkResult)
        std::cout << "Equivalent." << std::endl;
    else
        std::cout << "Not equivalent." << std::endl;

    std::cout << "Max node count: " << _maxNodeCount << std::endl;
    std::cout << "r: "<< _U->_r << std::endl;
}

/**Function*************************************************************

  Synopsis    [Initialize a identity matrix to a tensor.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void Checker::initTensorToIdentityMatrix(Tensor *tensor)
{
	DdNode *tmp1, *tmp2, *tmp3;
    int n = tensor->_rank / 2;

    for (int i = 0; i < tensor->_r; ++i)
    {
        for(int j = 0; j < _w; ++j)
        {
            if(i == 0 && j == _w-1)
            {
                tensor->_allBDD[j][i] = Cudd_ReadOne(_ddManager);
                Cudd_Ref(tensor->_allBDD[j][i]);
                for(int k = 0; k < n; ++k)
                {
                    tmp1 = Cudd_bddAnd(_ddManager, Cudd_Not(Cudd_bddIthVar(_ddManager, k)), Cudd_Not(Cudd_bddIthVar(_ddManager, k + n)));
                    Cudd_Ref(tmp1);
                    tmp2 = Cudd_bddAnd(_ddManager, Cudd_bddIthVar(_ddManager, k), Cudd_bddIthVar(_ddManager, k + n));
                    Cudd_Ref(tmp2);
                    tmp3 = Cudd_bddOr(_ddManager, tmp1, tmp2);
                    Cudd_Ref(tmp3);
                    Cudd_RecursiveDeref(_ddManager, tmp1);
                    Cudd_RecursiveDeref(_ddManager, tmp2);
                    tmp1 = Cudd_bddAnd(_ddManager, tensor->_allBDD[j][i], tmp3);
                    Cudd_Ref(tmp1);
                    Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[j][i]);
                    Cudd_RecursiveDeref(_ddManager, tmp3);
                    tensor->_allBDD[j][i] = tmp1;
                }
            }
            else
            {
                tensor->_allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(tensor->_allBDD[j][i]);
            }
        }
    }
}
