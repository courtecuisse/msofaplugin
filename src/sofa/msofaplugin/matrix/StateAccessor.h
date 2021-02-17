#pragma once

#include <sofa/msofaplugin/matrix/MechanicalVector.h>

namespace sofa::msofaplugin::matrix {

template<class DataTypes>
class StateAccessor : public TStateAccessor<typename DataTypes::VecReal> {
public:

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef Data<VecDeriv> DataVecDeriv;

    SOFA_CLASS(SOFA_TEMPLATE(StateAccessor,DataTypes),SOFA_TEMPLATE(TStateAccessor,VecReal));

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

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T*, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = sofa::core::objectmodel::New<T>();
        if (context) context->addObject(obj);
        if (arg) obj->parse(arg);
        return obj;
    }


    virtual std::string getTemplateName() const override {
        return templateName(this);
    }

    static std::string templateName(const StateAccessor<DataTypes>* = NULL) {
        return DataTypes::Name();
    }

    StateAccessor()
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
    sofa::core::objectmodel::SingleLink<StateAccessor<DataTypes>,sofa::core::behavior::MechanicalState<DataTypes>,BaseLink::FLAG_STRONGLINK> l_state;

};

}
