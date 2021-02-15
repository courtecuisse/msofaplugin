#pragma once

#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <thread>

namespace sofa::msofaplugin::matrix {

class BuildMVisitor : public simulation::MechanicalVisitor {
public:

    BuildMVisitor(component::linearsolver::DefaultMultiMatrixAccessor * mA)
    : simulation::MechanicalVisitor(core::MechanicalParams::defaultInstance())
    , m_accessor(mA) {
        mparams = *core::MechanicalParams::defaultInstance();
        mparams.setMFactor(1);mparams.setBFactor(0);mparams.setKFactor(0);
    }

    const char* getClassName() const override { return "BuildMVisitor"; }

    Result fwdMass(simulation::Node* /*node*/, core::behavior::BaseMass* ff) override {
        ff->addMToMatrix(&mparams,m_accessor);

        return RESULT_CONTINUE;
    }

    bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* map) override {
        return !map->areMatricesMapped();
    }
private:
    component::linearsolver::DefaultMultiMatrixAccessor * m_accessor;
    core::MechanicalParams mparams;
};

class BuildKVisitor : public simulation::MechanicalVisitor {
public:

    BuildKVisitor(component::linearsolver::DefaultMultiMatrixAccessor * mA)
    : simulation::MechanicalVisitor(core::MechanicalParams::defaultInstance())
    , m_accessor(mA) {
        mparams = *core::MechanicalParams::defaultInstance();
        mparams.setMFactor(0);mparams.setBFactor(0);mparams.setKFactor(1);
    }

    const char* getClassName() const override { return "BuildKVisitor"; }

    Result fwdForceField(simulation::Node* /*node*/, core::behavior::BaseForceField* ff) override {
        if (dynamic_cast<core::behavior::BaseMass*>(ff)) return RESULT_CONTINUE;

        ff->addKToMatrix(&mparams,m_accessor);

        return RESULT_CONTINUE;
    }

    bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* map) override {
        return !map->areMatricesMapped();
    }
private:
    component::linearsolver::DefaultMultiMatrixAccessor * m_accessor;
    core::MechanicalParams mparams;

};

class BuildMBKVisitor : public simulation::MechanicalVisitor {
public:

    BuildMBKVisitor(component::linearsolver::DefaultMultiMatrixAccessor * mA)
    : simulation::MechanicalVisitor(core::MechanicalParams::defaultInstance())
    , m_accessor(mA) {
        mparams = *core::MechanicalParams::defaultInstance();
        mparams.setMFactor(0);mparams.setBFactor(0);mparams.setKFactor(1);
    }

    const char* getClassName() const override { return "BuildKVisitor"; }

    Result fwdForceField(simulation::Node* /*node*/, core::behavior::BaseForceField* ff) override {
        if (dynamic_cast<core::behavior::BaseMass*>(ff)) return RESULT_CONTINUE;

        ff->addKToMatrix(&mparams,m_accessor);

        return RESULT_CONTINUE;
    }

    bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* map) override {
        return !map->areMatricesMapped();
    }
private:
    component::linearsolver::DefaultMultiMatrixAccessor * m_accessor;
    core::MechanicalParams mparams;

};

class ProjectMatrixVisitor : public simulation::MechanicalVisitor {
public:

    ProjectMatrixVisitor(component::linearsolver::DefaultMultiMatrixAccessor * mA)
    : simulation::MechanicalVisitor(core::MechanicalParams::defaultInstance())
    , m_accessor(mA) {}

    const char* getClassName() const override { return "ProjectMatrixVisitor"; }

    Result fwdProjectiveConstraintSet(simulation::Node* /*node*/, core::behavior::BaseProjectiveConstraintSet * c) override {
        c->applyConstraint(core::MechanicalParams::defaultInstance(), m_accessor);

        return RESULT_CONTINUE;
    }

    bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* map) override {
        return !map->areMatricesMapped();
    }
private:
    component::linearsolver::DefaultMultiMatrixAccessor * m_accessor;
};

}
