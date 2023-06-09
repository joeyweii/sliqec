#ifndef _CHECKER_H_
#define _CHECKER_H_

#include "bddSystem.h"

class Checker : public BDDSystem
{
    public:
        // Constructor and Destructor
        explicit Checker
            (
             int nQubits,
             int fInitBitWidth,
             int fBitWidthMode,
             bool fReorder
            );

        ~Checker();
        void addElementToOutputJSON(const std::string key, const std::string value);
        void printOutputJSON() const;
        void check(const Circuit *circuitU, const Circuit * circuitV);

    private:
        Tensor *_U;
        Tensor *_V;
        std::vector<std::pair<std::string, std::string>> _outputJSON;
        void constructUandV(const Circuit *circuitU, const Circuit *circuitV);
        void initTensorToIdentityMatrix(Tensor *tensor);
};

#endif
