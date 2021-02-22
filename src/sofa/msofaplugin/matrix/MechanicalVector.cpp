#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::msofaplugin::matrix {

class GenericStateAccessor : public BaseStateAccessor {
public:

    SOFA_CLASS(GenericStateAccessor,BaseStateAccessor);

    GenericStateAccessor(core::behavior::BaseMechanicalState * s) : m_state(s) {
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

BaseStateAccessor::SPtr BaseStateAccessor::create(core::behavior::BaseMechanicalState * state) {
    auto slaves = state->getSlaves();

    for (unsigned i=0;i<slaves.size();i++) {
        if (BaseStateAccessor::SPtr vector = sofa::core::objectmodel::SPtr_dynamic_cast<BaseStateAccessor>(slaves[i])) return vector;
    }

    sofa::core::objectmodel::BaseObjectDescription arg;
    arg.setAttribute("type",std::string("StateAccessor"));
    arg.setAttribute("template",state->getTemplateName());
    arg.setAttribute("name",state->getName() + "_MV");
    arg.setAttribute("mstate",std::string("@")+state->getPathName());

    sofa::core::objectmodel::BaseObject::SPtr obj = sofa::core::ObjectFactory::getInstance()->createObject(state->getContext(), &arg);
    BaseStateAccessor::SPtr vector = sofa::core::objectmodel::SPtr_dynamic_cast<BaseStateAccessor>(obj);
    if (vector != NULL) return vector;

    return sofa::core::objectmodel::New<GenericStateAccessor>(state);
}





}
