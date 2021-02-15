#pragma once

#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/helper/map.h>

namespace sofa::msofaplugin::solver {

/// Empty class used for default solver implementation without multi-threading support
class NoThreadManager {
public:
    static bool isAsyncSolver() { return false; }

    static inline std::string Name() { return ""; }
};

using namespace sofa::msofaplugin::matrix;

class MBaseLinearSolver : public sofa::core::behavior::LinearSolver {
public:

    SOFA_ABSTRACT_CLASS(MBaseLinearSolver,sofa::core::behavior::LinearSolver);

    Data<unsigned> d_updateStep;

    MBaseLinearSolver()
    : d_updateStep(initData(&d_updateStep,(unsigned) 1,"updateStep","Number of step before next update of the matrix (0 is computed only at initial step)")) {
        m_updateStep = 0;
    }

    virtual void resetSystem() override {}

    // set the incoming matrix and inverse the system to prepare for the solve
    virtual void setSystemMBKMatrix(const core::MechanicalParams* mparams) override {
        m_params=*mparams;

        if (m_updateStep == 0) {
            buildSystemMatrix();
        } else if(d_updateStep.getValue() != 0 && m_updateStep>=d_updateStep.getValue()) {
            buildSystemMatrix();
            m_updateStep = 0;
        }

        m_updateStep++;
    }

    void setSystemRHVector(core::MultiVecDerivId v) override {
        m_rh = v;
    }

    void setSystemLHVector(core::MultiVecDerivId v) override {
        m_lh = v;
    }

    void solveSystem() override {
        BaseSystemMatrix::MechanicalVectorId x = getMatrix()->getMVecId(this->m_lh,false);
        BaseSystemMatrix::MechanicalVectorId b = getMatrix()->getMVecId(this->m_rh);

        solve(x,b);
    }

    virtual void mult(core::MechanicalParams * mparams, bool acc = false) {
        BaseSystemMatrix::MechanicalVectorId x = getMatrix()->getMVecId(this->m_lh, acc);
        BaseSystemMatrix::MechanicalVectorId b = getMatrix()->getMVecId(this->m_rh);

        getMatrix()->mult(mparams,x,b,acc);
    }

    virtual void solve(BaseSystemMatrix::MechanicalVectorId & x,BaseSystemMatrix::MechanicalVectorId & b) = 0;

    virtual void buildSystemMatrix() = 0;

    virtual BaseSystemMatrix * getMatrix() = 0;

    const core::MechanicalParams * getMParams() const {
        return &m_params;
    }

//    virtual void setWi(defaulttype::BaseMatrix * res) = 0;

protected:
    core::MultiVecDerivId m_rh;
    core::MultiVecDerivId m_lh;
    core::MechanicalParams m_params;
    unsigned m_updateStep;


};

}
