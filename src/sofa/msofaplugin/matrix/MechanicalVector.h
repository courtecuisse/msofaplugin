#pragma once

#include <sofa/defaulttype/BaseVector.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/behavior/MultiVec.h>

#if 1 // Activate timers
#define DEBUG_PRINT(A) std::cout << A << std::endl;
#define TIMER_START(A) auto A##timer = std::chrono::high_resolution_clock::now()
#define TIMER_END(A) double A = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - A##timer).count()/1000.0
#define TIMER_NEXT(A,B) \
    auto A##B##timer = std::chrono::high_resolution_clock::now(); \
    double A = std::chrono::duration_cast<std::chrono::microseconds>(A##B##timer - A##timer).count()/1000.0; \
    auto B##timer = A##B##timer
#define TIMER_PRINT(A) std::cout << A << std::endl
#define TIMER_DO(A) A
#else
#define DEBUG_PRINT(A)
#define TIMER_START(A)
#define TIMER_END(A)
#define TIMER_NEXT(A,B)
#define TIMER_PRINT(A)
#define TIMER_DO(A)
#endif

namespace sofa::msofaplugin::matrix {

template<class VecReal>
class RealType {
public:
    typedef typename VecReal::value_type Real;
};

template<class T>
class RealType<helper::vector<T>> {
public:
    typedef T Real;
};

class BaseStateVectorOperator : public sofa::core::objectmodel::BaseObject {
public:

    SOFA_ABSTRACT_CLASS(BaseStateVectorOperator,sofa::core::objectmodel::BaseObject);

    virtual core::behavior::BaseMechanicalState * getState() = 0;

    virtual unsigned getGlobalDim() const = 0;

};

template<class VecReal>
class WriteAccessor {
public:
    typedef std::shared_ptr<WriteAccessor<VecReal>> SPtr;
    typedef typename RealType<VecReal>::Real Real;

    virtual Real * data() = 0;
};

template<class VecReal>
class ReadAccessor {
public:
    typedef std::shared_ptr<ReadAccessor<VecReal>> SPtr;
    typedef typename RealType<VecReal>::Real Real;

    virtual const Real * data() = 0;
};

template<class VecReal>
class StateVectorAccessor {
public:
    virtual typename WriteAccessor<VecReal>::SPtr write(unsigned id) = 0;

    virtual typename ReadAccessor<VecReal>::SPtr read(unsigned id) const = 0;

    virtual core::behavior::BaseMechanicalState * getState() = 0;
};

template<class VecReal>
class MechanicalVector {
public:
    typedef std::shared_ptr<MechanicalVector<VecReal> > SPtr;

    MechanicalVector(unsigned id)
    : m_id(id) {}

    virtual void syncToState() const = 0;

    virtual void syncToLocal() const = 0;

    virtual void validateLocal() const = 0;

    virtual void validateState() const = 0;

    virtual void invalidateLocal() const = 0;

    virtual void invalidateState() const = 0;

    virtual void update(const std::vector<BaseStateVectorOperator::SPtr> & vop) = 0;

    virtual typename WriteAccessor<VecReal>::SPtr write() = 0;

    virtual typename ReadAccessor<VecReal>::SPtr read() const = 0;

    unsigned id() const { return m_id; }

protected:
    const unsigned m_id;
};

}
