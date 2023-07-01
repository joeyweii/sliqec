#include "qcCheck.h"
#include "bddSystem.h"
#include "cudd.h"
#include <unistd.h>

void Checker::checkByConstructFunctionality(const Circuit *circuitU,
                                            const Circuit *circuitV) {
    int nQubits = circuitU->getNumberQubits();
    assert(circuitV->getNumberQubits == nQubits);

    initIdentityBDD(nQubits);
    initZeroBDD();

    Tensor *U = newTensor(nQubits * 2);
    Tensor *V = newTensor(nQubits * 2);

    bool fTranspose = false;

    initTensorToIdentityMatrix(U);
    for (int i = 0; i < circuitU->getGateCount(); ++i)
        applyGate(circuitU->getGate(i), U, fTranspose);

    initTensorToIdentityMatrix(V);
    for (int i = 0; i < circuitV->getGateCount(); ++i)
        applyGate(circuitV->getGate(i), V, fTranspose);

    if (eqCheckTwoTensor(U, V)) {
        addElementToOutputJSON("equivalence", "equivalent");
        deleteTensor(U);
        deleteTensor(V);
    } else {
        DdNode *pos = pickOneNonzeroCommonEntryPosition(U);

        Tensor *Un = constructAlphaIdentity(V, pos);
        deleteTensor(V);
        Tensor *Vn = constructAlphaIdentity(U, pos);
        deleteTensor(U);
        Cudd_RecursiveDeref(_ddManager, pos);

        for (int i = 0; i < circuitU->getGateCount(); ++i)
            applyGate(circuitU->getGate(i), Un, fTranspose);
        for (int i = 0; i < circuitV->getGateCount(); ++i)
            applyGate(circuitV->getGate(i), Vn, fTranspose);

        if (eqCheckTwoTensor(Un, Vn))
            addElementToOutputJSON("equivalence",
                                   "equivalent_up_to_global_phase");
        else
            addElementToOutputJSON("equivalence", "not_equivalent");

        deleteTensor(Un);
        deleteTensor(Vn);
    }

    addElementToOutputJSON("max_num_nodes",
                           std::to_string(Cudd_ReadPeakNodeCount(_ddManager)));

    Cudd_RecursiveDeref(_ddManager, _zeroBDD);
    Cudd_RecursiveDeref(_ddManager, _identityBDD);
}

void Checker::checkBySimulation(const Circuit *circuitU,
                                const Circuit *circuitV) {
    int nQubits = circuitU->getNumberQubits();
    assert(circuitV->getNumberQubits == nQubits);

    Tensor *U = newTensor(nQubits);
    Tensor *V = newTensor(nQubits);

    bool fTranspose = false;

    initTensorToBasisState(U);
    for (int i = 0; i < circuitU->getGateCount(); ++i)
        applyGate(circuitU->getGate(i), U, fTranspose);

    initTensorToBasisState(V);
    for (int i = 0; i < circuitV->getGateCount(); ++i)
        applyGate(circuitV->getGate(i), V, fTranspose);

    bool checkResult = eqCheckTwoTensor(U, V);
    if (checkResult)
        addElementToOutputJSON("equivalence", "probably_equivalent");
    else
        addElementToOutputJSON("equivalence", "not_equivalent");

    addElementToOutputJSON("max_num_nodes",
                           std::to_string(Cudd_ReadPeakNodeCount(_ddManager)));

    deleteTensor(U);
    deleteTensor(V);
}

void Checker::checkByConstructMiter(Circuit *circuitU,
                                    const Circuit *circuitV) {
    int nQubits = circuitU->getNumberQubits();
    assert(circuitV->getNumberQubits == nQubits);

    initIdentityBDD(nQubits);
    initZeroBDD();

    circuitU->daggerAllGate();
    circuitU->liftAllGateQubits(nQubits);

    Tensor *miter = newTensor(nQubits * 2);
    Tensor *identityMatrix = newTensor(nQubits * 2);
    initTensorToIdentityMatrix(miter);
    initTensorToIdentityMatrix(identityMatrix);

    int sizeU = circuitU->getGateCount(), sizeV = circuitV->getGateCount();
    int idxU = 0, idxV = 0;

    while (idxU < sizeU || idxV < sizeV) {
        if (idxU < sizeU) {
            applyGate(circuitU->getGate(idxU), miter, true);
            ++idxU;
        }

        while (idxV < sizeV && idxV * sizeU < idxU * sizeV) {
            applyGate(circuitV->getGate(idxV), miter, false);
            ++idxV;
        }
    }

    if (eqCheckTwoTensor(miter, identityMatrix))
        addElementToOutputJSON("equivalence", "equivalent");
    else if (checkIsTensorIdentityGlobalPhase(miter))
        addElementToOutputJSON("equivalence", "equivalent_up_to_global_phase");
    else
        addElementToOutputJSON("equivalence", "not_equivalent");

    addElementToOutputJSON("max_num_nodes",
                           std::to_string(Cudd_ReadPeakNodeCount(_ddManager)));

    circuitU->daggerAllGate();
    circuitU->unliftAllGateQubits(nQubits);
    deleteTensor(miter);
    deleteTensor(identityMatrix);
    Cudd_RecursiveDeref(_ddManager, _identityBDD);
    Cudd_RecursiveDeref(_ddManager, _zeroBDD);
}

void Checker::initTensorToIdentityMatrix(Tensor *tensor) {
    for (int i = 0; i < tensor->_r; ++i) {
        for (int j = 0; j < _w; ++j) {
            if (i == 0 && j == _w - 1) {
                tensor->_allBDD[j][i] = _identityBDD;
                Cudd_Ref(tensor->_allBDD[j][i]);
            } else {
                tensor->_allBDD[j][i] = _zeroBDD;
                Cudd_Ref(tensor->_allBDD[j][i]);
            }
        }
    }
}

void Checker::initTensorToBasisState(Tensor *tensor) {
    std::vector<bool> basisState(tensor->_rank, false);

    DdNode *var, *tmp;

    for (int i = 0; i < tensor->_r; i++) {
        if (i == 0) {
            for (int j = 0; j < _w - 1; j++) {
                tensor->_allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(tensor->_allBDD[j][i]);
            }
            tensor->_allBDD[_w - 1][i] = Cudd_ReadOne(_ddManager);
            Cudd_Ref(tensor->_allBDD[_w - 1][i]);
            for (int j = tensor->_rank - 1; j >= 0; j--)

            {
                var = Cudd_bddIthVar(_ddManager, j);
                if (!basisState[j])
                    tmp = Cudd_bddAnd(
                        _ddManager, Cudd_Not(var), tensor->_allBDD[_w - 1][i]);
                else
                    tmp = Cudd_bddAnd(
                        _ddManager, var, tensor->_allBDD[_w - 1][i]);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[_w - 1][i]);
                tensor->_allBDD[_w - 1][i] = tmp;
            }
        } else {
            for (int j = 0; j < _w; j++) {
                tensor->_allBDD[j][i] = Cudd_Not(Cudd_ReadOne(_ddManager));
                Cudd_Ref(tensor->_allBDD[j][i]);
            }
        }
    }
}

bool Checker::checkIsTensorIdentityGlobalPhase(Tensor *tensor) {
    for (int i = 0; i < _w; ++i) {
        for (int j = 0; j < tensor->_r; ++j) {
            if ((tensor->_allBDD[i][j] != _identityBDD) &&
                (tensor->_allBDD[i][j] != _zeroBDD))
                return false;
        }
    }
    return true;
}

void Checker::addElementToOutputJSON(const std::string key,
                                     const std::string value) {
    _outputJSON.push_back(std::make_pair(key, value));
}

void Checker::printOutputJSON() const {
    std::cout << "{\n";
    for (size_t i = 0; i < _outputJSON.size(); ++i) {
        auto &p = _outputJSON[i];
        std::cout << "  ";
        std::cout << "\"" << p.first << "\"";
        std::cout << ": ";
        std::cout << "\"" << p.second << "\"";
        if (i != _outputJSON.size() - 1) std::cout << ',';
        std::cout << '\n';
    }
    std::cout << "}\n";
}

void Checker::initIdentityBDD(int nQubits) {
    // shuffle BDD into a better order for identityBDD
    int *permut = new int[2 * nQubits];
    for (int i = 0; i < nQubits; i++) {
        permut[2 * i] = i;
        permut[2 * i + 1] = i + nQubits;
    }
    Cudd_ShuffleHeap(_ddManager, permut);
    delete[] permut;

    DdNode *tmp1, *tmp2, *tmp3;
    _identityBDD = Cudd_ReadOne(_ddManager);
    Cudd_Ref(_identityBDD);
    for (int i = nQubits - 1; i >= 0; --i) {
        tmp1 = Cudd_bddAnd(_ddManager,
                           Cudd_Not(Cudd_bddIthVar(_ddManager, i)),
                           Cudd_Not(Cudd_bddIthVar(_ddManager, i + nQubits)));
        Cudd_Ref(tmp1);
        tmp2 = Cudd_bddAnd(_ddManager,
                           Cudd_bddIthVar(_ddManager, i),
                           Cudd_bddIthVar(_ddManager, i + nQubits));
        Cudd_Ref(tmp2);
        tmp3 = Cudd_bddOr(_ddManager, tmp1, tmp2);
        Cudd_Ref(tmp3);
        Cudd_RecursiveDeref(_ddManager, tmp1);
        Cudd_RecursiveDeref(_ddManager, tmp2);
        tmp1 = Cudd_bddAnd(_ddManager, _identityBDD, tmp3);
        Cudd_Ref(tmp1);
        Cudd_RecursiveDeref(_ddManager, _identityBDD);
        Cudd_RecursiveDeref(_ddManager, tmp3);
        _identityBDD = tmp1;
    }
}

void Checker::initZeroBDD() {
    _zeroBDD = Cudd_Not(Cudd_ReadOne(_ddManager));
}

DdNode *Checker::pickOneNonzeroCommonEntryPosition(Tensor *tensor) {
    DdNode **vars = new DdNode *[tensor->_rank];
    for (int k = 0; k < tensor->_rank; ++k) {
        vars[k] = Cudd_bddIthVar(_ddManager, k);
    }

    DdNode *pos = nullptr;

    for (int i = 0; i < _w; ++i) {
        for (int j = 0; j < tensor->_r; ++j) {
            if (tensor->_allBDD[i][j] != _zeroBDD) {
                pos = Cudd_bddPickOneMinterm(
                    _ddManager, tensor->_allBDD[i][j], vars, tensor->_rank);
                Cudd_Ref(pos);
                break;
            }
        }
    }
    assert(pos != nullptr);
    delete[] vars;
    return pos;
}

Tensor *Checker::constructAlphaIdentity(Tensor *tensor, DdNode *pos) {
    Tensor *alphaI = newTensor(tensor->_rank);
    alphaI->_r = tensor->_r;
    alphaI->_k = tensor->_k;
    alphaI->_allBDD = std::vector(4, std::deque<DdNode *>(alphaI->_r));
    for (int i = 0; i < _w; ++i) {
        for (int j = 0; j < tensor->_r; ++j) {
            DdNode *curBDD = tensor->_allBDD[i][j];
            if (Cudd_Cofactor(_ddManager, curBDD, pos) != _zeroBDD)
                alphaI->_allBDD[i][j] = _identityBDD;
            else
                alphaI->_allBDD[i][j] = _zeroBDD;
            Cudd_Ref(alphaI->_allBDD[i][j]);
        }
    }
    return alphaI;
}
