#pragma once

#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <sofa/msofaplugin/matrix/StateSynchronizer.h>

namespace sofa::msofaplugin::matrix {

class GenericStateVectorOperator : public BaseStateVectorOperator {
public:

    SOFA_CLASS(GenericStateVectorOperator,BaseStateVectorOperator);

    GenericStateVectorOperator(core::behavior::BaseMechanicalState * s) : m_state(s) {
        m_state->addSlave(this);
    }

    virtual core::behavior::BaseMechanicalState * getState() override {
        return m_state;
    }

    virtual unsigned getGlobalDim() const override {
        return m_state->getSize() * m_state->getDerivDimension();
    }

private:
    core::behavior::BaseMechanicalState * m_state;
};

template<class DataTypes>
class StateVectorOperator : public BaseStateVectorOperator, public StateVectorAccessor<typename DataTypes::VecReal> {
public:

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef Data<VecDeriv> DataVecDeriv;

    SOFA_CLASS(SOFA_TEMPLATE(StateVectorOperator,DataTypes),BaseStateVectorOperator);

    /// Pre-construction check method called by ObjectFactory.
    /// Check that DataTypes matches the MechanicalState.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg) {
        std::string link = arg->getAttribute("mstate");
        std::string path(link.begin()+1,link.end());

        sofa::core::behavior::MechanicalState<DataTypes>* state = NULL;
        context->get(state,path);

        if (state == NULL) return false;

        return sofa::core::objectmodel::BaseObject::canCreate(obj, context, arg);
    }

    virtual std::string getTemplateName() const override {
        return templateName(this);
    }

    static std::string templateName(const StateVectorOperator<DataTypes>* = NULL) {
        return DataTypes::Name();
    }

    StateVectorOperator()
    : l_state(initLink("mstate", "MechanicalState used by this vector")) {}

    void parse ( core::objectmodel::BaseObjectDescription* arg ) override {
        core::objectmodel::BaseObject::parse(arg);
        l_state->addSlave(this);
    }

    core::behavior::BaseMechanicalState * getState() {
        return l_state.get();
    }

    unsigned getGlobalDim() const {
        return l_state->getSize() * Deriv::size();
    }

    // Warning this vector will not provide the correct size of the vector and may create memory issues it must only be accessed by the DataSynchronizer
    typename WriteAccessor<VecReal>::SPtr write(unsigned id) {
        core::VecDerivId vid = id;
        return typename WriteAccessor<VecReal>::SPtr(new TWriteAccessor<VecReal,VecDeriv>(l_state->write(vid)->beginEdit()));
    }

    // Warning this vector will not provide the correct size of the vector and may create memory issues it must only be accessed by the DataSynchronizer
    typename ReadAccessor<VecReal>::SPtr read(unsigned id) const {
        core::VecDerivId vid = id;
        return typename ReadAccessor<VecReal>::SPtr(new TReadAccessor<VecReal,VecDeriv>(&l_state->read(vid)->getValue()));
    }

protected:
    sofa::core::objectmodel::SingleLink<StateVectorOperator<DataTypes>,sofa::core::behavior::MechanicalState<DataTypes>,BaseLink::FLAG_STRONGLINK> l_state;

};

}
