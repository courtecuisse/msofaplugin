#pragma once

#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/vector.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <map>

namespace sofa::msofaplugin::matrix {

class MUnbuiltMatrix : public BaseSystemMatrix {
public:

    SOFA_CLASS(MUnbuiltMatrix,BaseSystemMatrix);

    virtual void clear() override {}

    virtual void buildMatrix() override {}

    virtual void mult(core::MechanicalParams * mparams, MechanicalVectorId & x, const MechanicalVectorId & b, bool acc = false) override {
//        TIMER_START(Mult);

        mparams->setDx(b.id());
        mparams->setDf(x);

        if (! acc) {
            simulation::MechanicalResetForceVisitor(mparams, x, false).execute(this->getContext());
//            simulation::MechanicalPropagateDxAndResetForceVisitor(&m_params, b, x, false).execute(m_context);
        }

        simulation::MechanicalAddMBKdxVisitor(mparams, x, acc).execute(this->getContext());
        simulation::MechanicalApplyConstraintsVisitor(mparams, x).execute(this->getContext());

//        TIMER_END(Mult);
//        TIMER_PRINT("Mult " << Mult << " ms");
    }

    double dot(MechanicalVectorId & a,MechanicalVectorId & b) override {
//        sofa::core::behavior::MultiVecDeriv A(&this->m_vop,a.getId()) ;
//        sofa::core::behavior::MultiVecDeriv B(&this->m_vop,b.getId()) ;

        return a.dot(b);
    }

    void eq(MechanicalVectorId & a,MechanicalVectorId & b) override {
//        sofa::core::behavior::MultiVecDeriv A(&this->m_vop,a.getId()) ;
//        sofa::core::behavior::MultiVecDeriv B(&this->m_vop,b.getId()) ;

        return a.eq(b);
    }

    void eq(MechanicalVectorId & a, MechanicalVectorId & b,MechanicalVectorId & c,double f) override {
//        sofa::core::behavior::MultiVecDeriv A(&this->m_vop,a.getId()) ;
//        sofa::core::behavior::MultiVecDeriv B(&this->m_vop,b.getId()) ;
//        sofa::core::behavior::MultiVecDeriv C(&this->m_vop,c.getId()) ;

        return a.eq(b,c,f);
    }

    void peq(MechanicalVectorId & a, MechanicalVectorId & b,double f) override {
//        sofa::core::behavior::MultiVecDeriv A(&this->m_vop,a.getId()) ;
//        sofa::core::behavior::MultiVecDeriv B(&this->m_vop,b.getId()) ;

        return a.peq(b,f);
    }

    virtual void freeMVector(MechanicalVectorId & v) {}

    virtual void syncMVector(MechanicalVectorId & v,bool copy) {}

    virtual unsigned colSize() const {
        return 0 ;
    }

    virtual unsigned rowSize() const {
        return 0 ;
    }

    virtual unsigned char getDomain (int x, int y) const {
        return 0 ;
    }
protected:
    core::MechanicalParams m_params;
};

}
