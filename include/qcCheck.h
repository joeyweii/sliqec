#pragma once

#include "bddSystem.h"

class Checker : public BDDSystem {
public:
    // Constructor and Destructor
    Checker(int nQubits = 0)
        : BDDSystem(2 * nQubits)
        , _identityBDD(nullptr)
        , _zeroBDD(nullptr) {}
    ~Checker() {}

    void addElementToOutputJSON(const std::string key, const std::string value);
    void printOutputJSON() const;

    void checkByConstructFunctionality(const Circuit *circuitU,
                                       const Circuit *circuitV);
    void checkByConstructMiter(Circuit *circuitU, const Circuit *circuitV);
    void checkBySimulation(const Circuit *circuitU, const Circuit *circuitV);

private:
    std::vector<std::pair<std::string, std::string>> _outputJSON;
    DdNode *_identityBDD;
    DdNode *_zeroBDD;

    void initIdentityBDD(int nQubits);
    void initZeroBDD();
    void initTensorToIdentityMatrix(Tensor *tensor);
    void initTensorToBasisState(Tensor *tensor);
    bool checkIsTensorIdentityGlobalPhase(Tensor *tensor);
    DdNode *pickOneNonzeroCommonEntryPosition(Tensor *tensor);
    Tensor *constructAlphaIdentity(Tensor *tensor, DdNode *minterm);
};
