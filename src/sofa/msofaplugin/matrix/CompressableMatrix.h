#pragma once

#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/msofaplugin/matrix/StateSynchronizer.h>
#include <sofa/msofaplugin/matrix/StateVectorOperator.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::msofaplugin::matrix {

template<class Real>
class CompressedMatrix {
public:
    virtual unsigned colSize() const = 0;

    virtual unsigned rowSize() const = 0;

    virtual const int * getRowBegin() const = 0;

    virtual const int * getColsIndex() const = 0;

    virtual const Real * getColsValue() const = 0;
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
            GetMechanicalStateVisitor(helper::vector<core::behavior::BaseMechanicalState*> & v)
            : sofa::simulation::BaseMechanicalVisitor(core::ExecParams::defaultInstance())
            , m_states(v) { m_states.clear(); }

            Result fwdMechanicalState(simulation::Node* /*node*/, core::behavior::BaseMechanicalState* mm) override {
                m_states.push_back(mm);
                return RESULT_CONTINUE;
            }

            std::vector<core::behavior::BaseMechanicalState*> & m_states;
        };
        //Get all the mstate concerned by the solver
        helper::vector<core::behavior::BaseMechanicalState*> states;
        GetMechanicalStateVisitor(states).execute(getContext());

        m_globalSize = 0;
        m_vectorOperations.clear();
        for (unsigned i=0;i<states.size();i++) {
            m_vectorOperations.push_back(createVectorOperations(states[i]));
            m_globalSize += states[i]->getSize() * states[i]->getDerivDimension();
        }

        for (unsigned i=0;i<m_mechanicalVectors.size();i++) {
            if (m_mechanicalVectors[i] != NULL)
                m_mechanicalVectors[i]->update(m_vectorOperations);
        }
    }

    virtual void postBuildMatrix() override {}

    virtual std::shared_ptr<CompressedMatrix<Real>> getCompressedMatrix(const core::MechanicalParams * mparams) = 0;

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
        auto vid = id.id().getId(m_vectorOperations[0]->getState());
        unsigned index = vid.getIndex();

        if (index >= m_mechanicalVectors.size()) m_mechanicalVectors.resize(index+1);

        if (m_mechanicalVectors[index] != NULL) return m_mechanicalVectors[index];

        if (m_vectorOperations.size() == 1) {
            if (auto vecAcc = dynamic_cast<StateVectorAccessor<VecReal> *>(m_vectorOperations.front().get())) {
                m_mechanicalVectors[index] = typename MechanicalVector<VecReal>::SPtr(new RawMechanicalVector<VecReal>(index,vecAcc));
                return m_mechanicalVectors[index];
            }
        }

//        m_mechanicalVectors[vid] = typename MechanicalVector<VecReal>::SPtr(new GenericSynchronizer<VecReal>(id,m_vectorOperations));
        m_mechanicalVectors[index] = typename MechanicalVector<VecReal>::SPtr(new SyncMechanicalVector<VecReal>(index,m_vectorOperations));
        return m_mechanicalVectors[index];
    }

//    virtual void print(std::ostream& out, const MechanicalVectorId & v) const override {
//        out << *(getMechanicalVector(v)->read());
//    }

//    virtual std::string getTemplateName() const override {
//        return MatrixType<VecReal,VecInt>::name();
//    }

//    static std::string templateName(const CompressableMatrix<TVecReal,TVecInt>* = NULL) {
//        return MatrixType<VecReal,VecInt>::name();
//    }

protected:
    mutable std::vector<typename MechanicalVector<VecReal>::SPtr> m_mechanicalVectors;
    std::vector<BaseStateVectorOperator::SPtr> m_vectorOperations;
    unsigned m_globalSize;

    BaseStateVectorOperator::SPtr createVectorOperations(core::behavior::BaseMechanicalState * state) {
        auto slaves = state->getSlaves();

        for (unsigned i=0;i<slaves.size();i++) {
            if (BaseStateVectorOperator::SPtr vector = sofa::core::objectmodel::SPtr_dynamic_cast<BaseStateVectorOperator>(slaves[i])) return vector;
        }

        sofa::core::objectmodel::BaseObjectDescription arg;
        arg.setAttribute("type",std::string("StateVectorOperator"));
        arg.setAttribute("template",state->getTemplateName());
        arg.setAttribute("name",state->getName() + "_MV");
        arg.setAttribute("mstate",std::string("@")+state->getPathName());

        sofa::core::objectmodel::BaseObject::SPtr obj = sofa::core::ObjectFactory::getInstance()->createObject(state->getContext(), &arg);
        BaseStateVectorOperator::SPtr vector = sofa::core::objectmodel::SPtr_dynamic_cast<BaseStateVectorOperator>(obj);
        if (vector != NULL) return vector;

        //TODO::
        return sofa::core::objectmodel::New<GenericStateVectorOperator>(state);
    }
};

}
