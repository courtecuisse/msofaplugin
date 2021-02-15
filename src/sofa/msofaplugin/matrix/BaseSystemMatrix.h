#pragma once

#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/vector.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/helper/AdvancedTimer.h>

namespace sofa::msofaplugin::matrix {

class BaseSystemMatrix : public core::objectmodel::BaseObject {
public:

    class MechanicalVectorId : public sofa::core::behavior::MultiVecDeriv {
        friend class BaseSystemMatrix;
    public:
        virtual ~MechanicalVectorId() {
            if (! dynamic) m_matrix->freeMVector(*this);
            delete this->vop;
        }

    protected:
        MechanicalVectorId(BaseSystemMatrix * m, core::MultiVecDerivId & v,bool sync = true)
        : sofa::core::behavior::MultiVecDeriv(new simulation::common::VectorOperations(core::ExecParams::defaultInstance(),m->getContext()),v)
        , m_matrix(m) {
            m_matrix->syncMVector(*this,sync);
        }

        MechanicalVectorId(BaseSystemMatrix * m)
        : sofa::core::behavior::MultiVecDeriv(new simulation::common::VectorOperations(core::ExecParams::defaultInstance(),m->getContext()))
        , m_matrix(m) {
            m_matrix->syncMVector(*this,false);            
        }

        BaseSystemMatrix * m_matrix;
    };

    SOFA_ABSTRACT_CLASS(BaseSystemMatrix, core::objectmodel::BaseObject);

    void buildSystemMatrix() {
        preBuildMatrix();
        clear();
        buildMatrix();
        postBuildMatrix();
    }

    //Default: Use addDforce implementation
    virtual void mult(core::MechanicalParams * mparams, MechanicalVectorId & x, const MechanicalVectorId & b, bool acc = false)  = 0;

    virtual double dot(MechanicalVectorId & a, MechanicalVectorId & b) = 0;

    virtual void eq(MechanicalVectorId & a, MechanicalVectorId & b) = 0;

    virtual void eq(MechanicalVectorId & a, MechanicalVectorId & b, MechanicalVectorId & c,double f) = 0;

    virtual void peq(MechanicalVectorId & a, MechanicalVectorId & b,double f) = 0;

    virtual void freeMVector(MechanicalVectorId & v) = 0;

    virtual void syncMVector(MechanicalVectorId & v,bool copy) = 0;

    virtual unsigned colSize() const = 0 ;

    virtual unsigned rowSize() const = 0 ;

    MechanicalVectorId getMVecId(core::MultiVecDerivId & v,bool sync = true) {
        return MechanicalVectorId(this,v,sync);
    }

    MechanicalVectorId createMVecId() {
        return MechanicalVectorId(this);
    }

protected:
    virtual void clear() = 0;

    virtual void buildMatrix() = 0;

    virtual void preBuildMatrix() {}

    virtual void postBuildMatrix() {}
};

}
