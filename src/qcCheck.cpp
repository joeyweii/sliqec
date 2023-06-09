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
     fInitBitWidth,
     fBitWidthControl,
     fReorder
    ),
    _nQubits(nQubits)
{
}

/**Function*************************************************************

  Synopsis    [Construct the unitary matrix of U and V and check their equivalence.]

  Description []

  SideEffects []

  SeeAlso     []

 ***********************************************************************/
void Checker::checkByConstructFunctionality(const Circuit *circuitU, const Circuit *circuitV)
{
    Tensor* _U = newTensor(_nQubits*2);
    Tensor* _V = newTensor(_nQubits*2);

    bool fTranspose = false;

    initTensorToIdentityMatrix(_U);
    for(int i = 0; i < circuitU->getGateCount(); ++i)
        applyGate(circuitU->getGate(i), _U, fTranspose);

    initTensorToIdentityMatrix(_V);
    for(int i = 0; i < circuitV->getGateCount(); ++i)
        applyGate(circuitV->getGate(i), _V, fTranspose);

    bool checkResult = eqCheckTwoTensor(_U, _V);
    if(checkResult)
        addElementToOutputJSON("equivalence", "equivalent");
    else
        addElementToOutputJSON("equivalence", "not_equivalent");

    addElementToOutputJSON("num_nodes", std::to_string(_maxNodeCount));

    deleteTensor(_U);
    deleteTensor(_V);
}

/**Function*************************************************************

  Synopsis    [Construct the miter of U and V and check if the miter equals to I.]

  Description []

  SideEffects []

  SeeAlso     []

 ***********************************************************************/
void Checker::checkByConstructMiter(Circuit *circuitU, const Circuit *circuitV)
{
    circuitU->daggerAllGate();
    circuitU->liftAllGateQubits(_nQubits);

    Tensor* miter = newTensor(_nQubits*2);
    Tensor* identityMatrix = newTensor(_nQubits*2);
    initTensorToIdentityMatrix(miter);
    initTensorToIdentityMatrix(identityMatrix);

    int sizeU = circuitU->getGateCount(), sizeV = circuitV->getGateCount();
    int ratio = sizeV / sizeU;
    int idxU = 0, idxV = 0;

    while(idxU < sizeU || idxV < sizeV)
    {
        if(idxU < sizeU)
        {
            applyGate(circuitU->getGate(idxU), miter, true);
            ++idxU;
        }

        int count = 0;
        while(idxV < sizeV && count < ratio)
        {
            applyGate(circuitV->getGate(idxV), miter, false);
            ++count;
            ++idxV;
        }
    }

    bool checkResult = eqCheckTwoTensor(miter, identityMatrix);
    if(checkResult)
        addElementToOutputJSON("equivalence", "equivalent");
    else
        addElementToOutputJSON("equivalence", "not_equivalent");

    addElementToOutputJSON("num_nodes", std::to_string(_maxNodeCount));

    circuitU->daggerAllGate();
    circuitU->unliftAllGateQubits(_nQubits);
    deleteTensor(miter);
    deleteTensor(identityMatrix);
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

/**Function*************************************************************

  Synopsis    [Add a pair of key/value strings into _outputJSON.]

  Description []

  SideEffects []

  SeeAlso     []

 ***********************************************************************/

void Checker::addElementToOutputJSON(const std::string key, const std::string value)
{
    _outputJSON.push_back(std::make_pair(key, value));
}

/**Function*************************************************************

  Synopsis    [Print the content of _outputJSON to the stdout.]

  Description []

  SideEffects []

  SeeAlso     []

 ***********************************************************************/

void Checker::printOutputJSON() const
{
    std::cout << "{\n";
    for(int i = 0; i < _outputJSON.size(); ++i)
    {
        auto &p = _outputJSON[i];
        std::cout << "  ";
        std::cout << "\"" << p.first << "\"";
        std::cout << ": ";
        std::cout << "\"" << p.second << "\"";
        if(i != _outputJSON.size()-1) std::cout << ',';
        std::cout << '\n';
    }
    std::cout << "}\n";
}
