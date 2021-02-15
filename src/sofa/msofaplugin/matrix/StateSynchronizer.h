#pragma once

#include <sofa/msofaplugin/matrix/MechanicalVector.h>

namespace sofa::msofaplugin::matrix {

template<class VecReal,class VecDeriv>
class TWriteAccessor : public WriteAccessor<VecReal> {
public:
    typedef typename RealType<VecReal>::Real Real;

    TWriteAccessor(VecDeriv * v) : m_vector(v) {}

    Real * data() override { return (Real*) m_vector->data();}

private:
    VecDeriv * m_vector;
};

template<class VecReal,class VecDeriv>
class TReadAccessor : public ReadAccessor<VecReal> {
public:
    typedef typename RealType<VecReal>::Real Real;

    TReadAccessor(const VecDeriv * v) : m_vector(v) {}

    const Real * data() override { return (const Real*) m_vector->data();}

private:
    const VecDeriv * m_vector;
};

class BaseStateSynchronizer {
public:
    BaseStateSynchronizer(unsigned offset,unsigned size)
    : m_offset(offset), m_size(size), m_stateValid(true), m_localValid(false) {}

    typedef std::unique_ptr<BaseStateSynchronizer> Uptr;

    virtual void syncToState() = 0;

    virtual void syncFromState() = 0;

    virtual core::behavior::BaseMechanicalState * getState() = 0;

    inline void validateState() { m_stateValid = true; }

    inline void validateLocal() { m_localValid = true; }

    inline void invalidateState() { m_stateValid = false; }

    inline void invalidateLocal() { m_localValid = false; }

protected:
    const unsigned m_offset,m_size;
    bool m_stateValid, m_localValid;
};

template<class VecReal1, class VecReal2>
class StateSynchronizer : public BaseStateSynchronizer {
public:

    typedef typename RealType<VecReal1>::Real Real1; // vector of the MState
    typedef typename RealType<VecReal2>::Real Real2; // copy vector in vecReal

    StateSynchronizer(StateVectorAccessor<VecReal1> * acc, unsigned id, VecReal2 * lv, unsigned offset,unsigned size)
    : BaseStateSynchronizer(offset,size)
    , m_accessor(acc)
    , m_id(id)
    , m_localVec(lv) {}

    virtual void syncToState() {
        if (m_stateValid) return;

        Real1 * mstateData = m_accessor->write(m_id)->data();
        const Real2 * mlocalData = m_localVec->data();
        for (unsigned i=0;i<m_size;i++) mstateData[i] = mlocalData[m_offset+i];
        m_stateValid = true;
    }

    virtual void syncFromState() {
        if (m_localValid) return;

        const Real1 * mstateData = m_accessor->read(m_id)->data();
        Real2 * mlocalData = (Real2 *) m_localVec->data();
        for (unsigned i=0;i<m_size;i++) mlocalData[m_offset+i] = mstateData[i];
        m_localValid = true;
    }

    inline core::behavior::BaseMechanicalState * getState() {
        return m_accessor->getState();
    }

protected:
    unsigned m_id;
    VecReal2 * m_localVec;
    StateVectorAccessor<VecReal1> * m_accessor;
};

template<class VecReal2>
class GenericSynchronizer : public BaseStateSynchronizer {
public:

    typedef typename RealType<VecReal2>::Real Real2;

    GenericSynchronizer(core::behavior::BaseMechanicalState * state, unsigned id, VecReal2 * lv, unsigned offset, unsigned size)
    : BaseStateSynchronizer(offset, size)
    , m_state(state)
    , m_id(id)
    , m_localVec(lv) {
        m_tmp.resize(size);
    }

    virtual void syncToState() {
        if (m_stateValid) return;

        core::VecDerivId vid = m_id;

        Real2 * mlocalData = (Real2 *) m_localVec->data();
        for (unsigned i=0;i<m_tmp.size();i++) m_tmp[i] = mlocalData[m_offset + i];
        m_state->copyFromBuffer(vid,m_tmp.data(),m_tmp.size());
        m_stateValid = true;
    }

    virtual void syncFromState() {
        if (m_localValid) return;

        core::VecDerivId vid = m_id;

        m_state->copyToBuffer(m_tmp.data(),vid,m_tmp.size());
        Real2 * mlocalData = (Real2 *) m_localVec->data();
        for (unsigned i=0;i<m_tmp.size();i++) mlocalData[m_offset + i] = m_tmp[i];
        m_localValid = true;
    }

    core::behavior::BaseMechanicalState * getState() {
        return m_state;
    }

protected:
    core::behavior::BaseMechanicalState * m_state;
    unsigned m_id;
    VecReal2 * m_localVec;
    helper::vector<SReal> m_tmp;
};

/**
 * Synchro mechanism when mstate is different
 */
template<class VecReal>
class SyncMechanicalVector : public MechanicalVector<VecReal> {
public:

    SyncMechanicalVector(unsigned id, const std::vector<BaseStateVectorOperator::SPtr> & vop)
    : MechanicalVector<VecReal>(id) {
        update(vop);
    }

    void syncToState() const override {
        for (unsigned i=0;i<m_synchrnoizer.size();i++) m_synchrnoizer[i]->syncToState();
    }

    void syncToLocal() const override {
        for (unsigned i=0;i<m_synchrnoizer.size();i++) m_synchrnoizer[i]->syncFromState();
    }

    void invalidateLocal() const override {
        for (unsigned i=0;i<m_synchrnoizer.size();i++) m_synchrnoizer[i]->invalidateLocal();
    }

    void invalidateState() const override {
        for (unsigned i=0;i<m_synchrnoizer.size();i++) m_synchrnoizer[i]->invalidateState();
    }

    void validateLocal() const override {
        for (unsigned i=0;i<m_synchrnoizer.size();i++) m_synchrnoizer[i]->validateLocal();
    }

    void validateState() const override {
        for (unsigned i=0;i<m_synchrnoizer.size();i++) m_synchrnoizer[i]->validateState();
    }

    virtual void update(const std::vector<BaseStateVectorOperator::SPtr> & vop) override {
        unsigned offset = 0;

        m_synchrnoizer.resize(vop.size());

        for (unsigned i=0;i<vop.size();i++) {
            unsigned size = vop[i]->getGlobalDim();
            if (m_synchrnoizer[i] == NULL || (m_synchrnoizer[i]->getState() != vop[i]->getState())) {
                m_synchrnoizer[i] = BaseStateSynchronizer::Uptr(createSynchronizer(vop[i].get(), offset, size));
            }

            offset+=size;
        }

        m_synchrnoizer.resize(vop.size());
        m_localData.fastResize(offset);
    }

    BaseStateSynchronizer * createSynchronizer(BaseStateVectorOperator * v, unsigned offset, unsigned size) {
        if (auto acc = dynamic_cast<StateVectorAccessor<helper::vector<float>> * >(v)) {
            return new StateSynchronizer<helper::vector<float>,VecReal>(acc, this->id(), &m_localData, offset, size);
        } else if (auto acc = dynamic_cast<StateVectorAccessor<helper::vector<double>> * >(v)) {
            return new StateSynchronizer<helper::vector<double>,VecReal>(acc, this->id(), &m_localData, offset, size);
        }

        return new GenericSynchronizer<VecReal>(v->getState(),this->id(), &m_localData,offset,size);
    }

    typename WriteAccessor<VecReal>::SPtr write() override {
//        syncToLocal();
        invalidateState();

        return typename WriteAccessor<VecReal>::SPtr(new TWriteAccessor<VecReal,VecReal>(&m_localData));
    }

    typename ReadAccessor<VecReal>::SPtr read() const override {
//        syncToLocal();
        return typename ReadAccessor<VecReal>::SPtr(new TReadAccessor<VecReal,VecReal>(&m_localData));
    }

protected:
    std::vector<BaseStateSynchronizer::Uptr> m_synchrnoizer;
    VecReal m_localData;
};

template<class VecReal>
class RawMechanicalVector : public MechanicalVector<VecReal> {
public:

    RawMechanicalVector(unsigned id, StateVectorAccessor<VecReal> * acc)
    : MechanicalVector<VecReal>(id)
    , m_accessor(acc) {}

    void syncToState() const override  {}

    void syncToLocal() const override  {}

    void validateLocal() const override  {}

    void validateState() const override  {}

    void invalidateLocal() const override  {}

    void invalidateState() const override  {}

    void update(const std::vector<BaseStateVectorOperator::SPtr> & /*vop*/) override {}

    typename WriteAccessor<VecReal>::SPtr write() override {
        return m_accessor->write(this->id());
    }

    typename ReadAccessor<VecReal>::SPtr read() const override {
        return m_accessor->read(this->id());
    }

protected:
    StateVectorAccessor<VecReal> * m_accessor;
};

}
