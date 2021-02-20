#pragma once

#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/msofaplugin/matrix/CompressableMatrix.h>
#include <sofa/msofaplugin/matrix/LocalIncomingSparseMatrix.h>
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
#include <stdlib.h>
#include <Eigen/Sparse>
#include <sofa/helper/AdvancedTimer.h>

namespace sofa::msofaplugin::matrix {

template<class Real>
class EigenCompressedMatrix : public CompressedMatrix<Real> {
public:
    EigenCompressedMatrix(Eigen::SparseMatrix<Real,Eigen::RowMajor> & mat) : m_matrix(mat) {}

    const int * getColptr() const { return m_matrix.outerIndexPtr(); }

    const int * getRowind() const { return m_matrix.innerIndexPtr(); }

    const Real * getValues() const { return m_matrix.valuePtr(); }

    virtual unsigned colSize() const { return m_matrix.cols(); }

    virtual unsigned rowSize() const { return m_matrix.rows(); }

    Eigen::SparseMatrix<Real,Eigen::RowMajor> & m_matrix;
};



template<class Real>
class ProjectIncomingMatrix : public defaulttype::BaseMatrix {
public:
    typedef defaulttype::BaseMatrix::Index Index;

    ProjectIncomingMatrix(helper::vector<int> & clearCols, helper::vector<int> & clearRows)
    : m_clearCols(clearCols)
    , m_clearRows(clearRows)
    , m_needRebuild(true)
    {}

    Index rowSize(void) const override { return m_rowSize; }

    Index colSize(void) const override { return m_colSize;}

    SReal element(Index , Index ) const override { return 0.0; }

    void resize(Index r, Index c) {
        m_rowSize = r;
        m_colSize = c;
    }

    void clear() {
        m_clearCols.clear();
        m_clearRows.clear();
        m_setVal.clear();
        m_needRebuild = true;
    }

    void set(Index i, Index j, double v) override {
        m_setVal.push_back(Eigen::Triplet<Real>(i,j,v));
    }

    inline void add( Index , Index , double ) override {}

    void clearRow(Index i) override { m_clearRows.push_back(i); }

    void clearCol(Index j) override { m_clearCols.push_back(j); }

    void clearRowCol(Index i) override {
        clearCol(i);
        clearRow(i);
    }

    void clearRows(Index imin, Index imax) override {
        for(Index i=imin; i<imax; i++)
            clearRow(i);
    }

    void clearCols(Index imin, Index imax) override {
        for(Index i=imin; i<imax; i++)
            clearCol(i);
    }

    void clearRowsCols(Index imin, Index imax) override {
        for(Index i=imin; i<imax; i++) {
            clearRowCol(i);
        }
    }

    bool needRebuild() {
        return m_needRebuild;
    }

    void buildMatrix(std::vector<BaseStateAccessor::SPtr> & vacc, core::objectmodel::BaseContext * ctx) {
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

        clear();

        component::linearsolver::DefaultMultiMatrixAccessor accessor;
        accessor.setGlobalMatrix(this);
        accessor.setupMatrices();

        for (unsigned i=0;i<vacc.size();i++)
            accessor.addMechanicalState(vacc[i]->getState());

        ProjectMatrixVisitor(&accessor).execute(ctx);

        m_P.resize(m_rowSize,m_colSize);
        m_P.setFromTriplets(m_setVal.begin(),m_setVal.end());

        m_needRebuild = false;
    }

    Eigen::SparseMatrix<Real,Eigen::RowMajor> & getMatrix() {
        return m_P;
    }

private :
    unsigned m_rowSize,m_colSize;
    helper::vector<int> & m_clearCols;
    helper::vector<int> & m_clearRows;
    helper::vector<Eigen::Triplet<Real>> m_setVal;
    Eigen::SparseMatrix<Real,Eigen::RowMajor> m_P;
    bool m_needRebuild;
};

template<class VecReal,class VecInt>
class BaseIncomingSparseMatrix : public CompressableMatrix<VecReal,VecInt> {
public:
    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE2(BaseIncomingSparseMatrix,VecReal,VecInt),SOFA_TEMPLATE2(CompressableMatrix, VecReal,VecInt) );

    typedef typename CompressableMatrix<VecReal,VecInt>::Real Real;

    template<class Real>
    class InternalCompressed : public CompressedMatrix<Real> {
    public:

        InternalCompressed(VecInt & colptr,VecInt & rowInd, VecReal & values, unsigned r, unsigned c)
        : m_colptr(colptr)
        , m_rowptr(rowInd)
        , m_values(values)
        , m_rowSize(r)
        , m_colSize(c) {}

        const int * getRowBegin() const override { return m_colptr.data(); }

        const int * getColsIndex() const override { return m_rowptr.data(); }

        const Real * getColsValue() const override { return m_values.data(); }

        virtual unsigned colSize() const override { return m_colSize; }

        virtual unsigned rowSize() const override { return m_rowSize; }

        VecInt & m_colptr;
        VecInt & m_rowptr;
        VecReal & m_values;
        unsigned m_rowSize,m_colSize;
    };

    virtual void clear() {}

    void buildMatrix() override {
        LocalIncomingSparseMatrix<VecReal,VecInt>::doCreateVisitor(m_matrices,this->getContext());

        if (m_projectionMatrix.needRebuild()) {
            m_projectionMatrix.resize(this->m_globalSize,this->m_globalSize);
            m_projectionMatrix.buildMatrix(this->m_stateAccessor,this->getContext());
        }

        for (unsigned i=0;i<m_matrices.size();i++) {
            if (m_matrices[i]->colSize() == 0 || m_matrices[i]->rowSize() == 0) {
                m_matrices[i]->resize(this->m_globalSize,this->m_globalSize);
                m_matrices[i]->buildMatrix(this->m_stateAccessor,m_clearCols,m_clearRows);
            } else {
                m_matrices[i]->fastReBuild();
            }
        }
    }

    typename CompressedMatrix<Real>::SPtr getCompressedMatrix(const core::MechanicalParams * mparams) {
        M_global = m_projectionMatrix.getMatrix();

        for (unsigned i=0;i<m_matrices.size();i++) {
            Eigen::Map< Eigen::SparseMatrix<Real,Eigen::RowMajor> > map(m_matrices[i]->colSize(),
                                                                        m_matrices[i]->rowSize(),
                                                                        m_matrices[i]->getNnz(),
                                                                        m_matrices[i]->getColptr().data(),
                                                                        m_matrices[i]->getRowind().data(),
                                                                        m_matrices[i]->getValues().data());

            M_global += map*m_matrices[i]->getFactor(mparams);
        }

        return typename CompressedMatrix<Real>::SPtr(new EigenCompressedMatrix<Real>(M_global));
    }

    virtual unsigned colSize() const {
        return m_colSize;
    }

    virtual unsigned rowSize() const {
        return m_rowSize;
    }

    BaseIncomingSparseMatrix()
    : m_projectionMatrix(m_clearCols,m_clearRows) {}

protected:
    unsigned m_rowSize,m_colSize;
    std::vector<typename LocalIncomingSparseMatrix<VecReal,VecInt>::SPtr> m_matrices;
    Eigen::SparseMatrix<Real,Eigen::RowMajor> M_global;
    helper::vector<int> m_clearCols,m_clearRows;
    ProjectIncomingMatrix<Real> m_projectionMatrix;

};

template<class Real>
struct IncomingEigenType;

template<> struct IncomingEigenType<double> {
    typedef Eigen::VectorXd Vector;
};

template<> struct IncomingEigenType<float> {
    typedef Eigen::VectorXf Vector;
};

template<class TReal>
class IncomingSparseMatrix : public BaseIncomingSparseMatrix<helper::vector<TReal>, helper::vector<int>> {
public:


    typedef TReal Real;
    typedef helper::vector<TReal> VecReal;
    typedef helper::vector<int> VecInt;

    SOFA_CLASS(SOFA_TEMPLATE(IncomingSparseMatrix,TReal),SOFA_TEMPLATE2(BaseIncomingSparseMatrix, VecReal, VecInt));

    void mult(core::MechanicalParams * mparams, BaseSystemMatrix::MechanicalVectorId & x, const BaseSystemMatrix::MechanicalVectorId & b, bool acc) override {
        auto X = this->getMechanicalVector(x)->write();
        auto B = this->getMechanicalVector(b)->read();

        TReal * x_ptr = X->data();
        const TReal * b_ptr = B->data();

        Eigen::Map<typename IncomingEigenType<TReal>::Vector> ex(x_ptr,this->m_globalSize);
        Eigen::Map<typename IncomingEigenType<TReal>::Vector> eb((TReal*) b_ptr,this->m_globalSize);

        if (acc) ex += this->m_projectionMatrix.getMatrix() * eb;
        else ex = this->m_projectionMatrix.getMatrix() * eb;

        for (unsigned i=0;i<this->m_matrices.size();i++) {
            Eigen::Map< Eigen::SparseMatrix<Real,Eigen::RowMajor> > M(this->m_matrices[i]->colSize(),
                                                                      this->m_matrices[i]->rowSize(),
                                                                      this->m_matrices[i]->getNnz(),
                                                                      this->m_matrices[i]->getColptr().data(),
                                                                      this->m_matrices[i]->getRowind().data(),
                                                                      this->m_matrices[i]->getValues().data());

            ex += M * eb * this->m_matrices[i]->getFactor(mparams);
        }
    }
};


}
