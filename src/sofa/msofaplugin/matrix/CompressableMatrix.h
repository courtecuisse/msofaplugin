#pragma once

#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>

namespace sofa::msofaplugin::matrix {

template<class Real>
class CompressedMatrix {
public:
    typedef std::shared_ptr<CompressedMatrix<Real> > SPtr;

    virtual unsigned colSize() const = 0;

    virtual unsigned rowSize() const = 0;

    virtual const int * getColptr() const = 0;

    virtual const int * getRowind() const = 0;

    virtual const Real * getValues() const = 0;
};


template<class TVecReal,class TVecInt>
class CompressableMatrix : public BaseSystemMatrix {
public:

    typedef typename RealType<TVecReal>::Real Real;
    typedef TVecReal VecReal;
    typedef TVecInt VecInt;

    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE2(CompressableMatrix,VecReal,VecInt), BaseSystemMatrix);

    virtual void preBuildMatrix() override {
        //Visitor to retreive the list of mstates
        class GetMechanicalStateVisitor : public sofa::simulation::BaseMechanicalVisitor {
        public:
            GetMechanicalStateVisitor(sofa::type::vector<core::behavior::BaseMechanicalState*> & v)
            : sofa::simulation::BaseMechanicalVisitor(core::ExecParams::defaultInstance())
            , m_states(v) { m_states.clear(); }

            Result fwdMechanicalState(simulation::Node* /*node*/, core::behavior::BaseMechanicalState* mm) override {
                m_states.push_back(mm);
                return RESULT_CONTINUE;
            }

            std::vector<core::behavior::BaseMechanicalState*> & m_states;
        };
        //Get all the mstate concerned by the solver
        sofa::type::vector<core::behavior::BaseMechanicalState*> states;
        GetMechanicalStateVisitor(states).execute(getContext());

        m_globalSize = 0;
        m_stateAccessor.clear();
        for (unsigned i=0;i<states.size();i++) {
            auto acc = TStateAccessor<VecReal>::create(states[i]);
            acc->f_printLog.setValue(this->f_printLog.getValue());
            m_stateAccessor.push_back(acc);
            m_globalSize += states[i]->getSize() * states[i]->getDerivDimension();
        }

        for (unsigned i=0;i<m_mechanicalVectors.size();i++) {
            if (m_mechanicalVectors[i] != NULL)
                m_mechanicalVectors[i]->update(m_stateAccessor);
        }
    }

    virtual void postBuildMatrix() override {}

    virtual typename CompressedMatrix<Real>::SPtr getCompressedMatrix(const core::MechanicalParams * mparams) = 0;

    double dot(MechanicalVectorId & a,MechanicalVectorId & b) override {
        auto Amv = getMechanicalVector(a);
        auto Bmv = getMechanicalVector(b);

        const Real * A = Amv->read()->data();
        const Real * B = Bmv->read()->data();

        double res = 0.0;
        for (unsigned i=0;i<m_globalSize;i++) res += A[i]*B[i];
        return res;
    }

    virtual void eq(MechanicalVectorId & a,MechanicalVectorId & b) override {
        auto Amv = getMechanicalVector(a);
        auto Bmv = getMechanicalVector(b);

        Real * A = Amv->write()->data();
        const Real * B = Bmv->read()->data();

        for (unsigned i=0;i<this->m_globalSize;i++) A[i] = B[i];
    }

    virtual void eq(MechanicalVectorId & a, MechanicalVectorId & b,MechanicalVectorId & c,double f) override {
        auto Amv = getMechanicalVector(a);
        auto Bmv = getMechanicalVector(b);
        auto Cmv = getMechanicalVector(c);

        Real * A = Amv->write()->data();
        const Real * B = Bmv->read()->data();
        const Real * C = Cmv->read()->data();

        for (unsigned i=0;i<this->m_globalSize;i++) A[i] = B[i] + C[i]*f;
    }

    virtual void peq(MechanicalVectorId & a, MechanicalVectorId & b,double f) override {
        auto Amv = getMechanicalVector(a);
        auto Bmv = getMechanicalVector(b);

        Real * A = Amv->write()->data();
        const Real * B = Bmv->read()->data();

        for (unsigned i=0;i<this->m_globalSize;i++) A[i] += B[i]*f;
    }

    virtual void freeMVector(MechanicalVectorId & v) {
        auto mv = getMechanicalVector(v);
        mv->syncToState();
        mv->invalidateLocal();
    }

    virtual void syncMVector(MechanicalVectorId & v,bool copy) {
        auto mv = getMechanicalVector(v);
        if (copy) mv->syncToLocal();
        mv->validateLocal();
    }

    typename MechanicalVector<VecReal>::SPtr getMechanicalVector(const BaseSystemMatrix::MechanicalVectorId & id) const {
        auto vid = id.id().getId(m_stateAccessor[0]->getState());
        unsigned index = vid.getIndex();

        if (index >= m_mechanicalVectors.size()) m_mechanicalVectors.resize(index+1);

        if (m_mechanicalVectors[index] != NULL) return m_mechanicalVectors[index];

        if (m_stateAccessor.size() == 1) {
            if (auto vecAcc = core::objectmodel::SPtr_dynamic_cast<TStateAccessor<VecReal>>(m_stateAccessor.front())) {
                m_mechanicalVectors[index] = typename MechanicalVector<VecReal>::SPtr(new RawMechanicalVector<VecReal>(index,vecAcc));
                return m_mechanicalVectors[index];
            }
        }

//        m_mechanicalVectors[vid] = typename MechanicalVector<VecReal>::SPtr(new GenericSynchronizer<VecReal>(id,m_stateAccessor));
        m_mechanicalVectors[index] = typename MechanicalVector<VecReal>::SPtr(new SyncMechanicalVector<VecReal>(index,m_stateAccessor));
        return m_mechanicalVectors[index];
    }

protected:
    mutable std::vector<typename MechanicalVector<VecReal>::SPtr> m_mechanicalVectors;
    std::vector<BaseStateAccessor::SPtr> m_stateAccessor;
    unsigned m_globalSize;


};

}
