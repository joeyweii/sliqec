#ifndef _CHECKER_H_
#define _CHECKER_H_

#include "bddSystem.h"

class Checker : public BDDSystem
{
public:
    // Constructor and Destructor
    explicit Checker(int nQubits,
                     int fInitBitWidth,
                     int fBitWidthMode,
                     bool fReorder);
    ~Checker() {}

    void addElementToOutputJSON(const std::string key, const std::string value);
    void printOutputJSON() const;

    void checkByConstructFunctionality(const Circuit *circuitU,
                                       const Circuit *circuitV);
    void checkByConstructMiter(Circuit *circuitU, const Circuit *circuitV);
    void checkBySimulation(const Circuit *circuitU, const Circuit *circuitV);

private:
    std::vector<std::pair<std::string, std::string>> _outputJSON;

    void initTensorToIdentityMatrix(Tensor *tensor);
    void initTensorToBasisState(Tensor *tensor);
};

#endif
