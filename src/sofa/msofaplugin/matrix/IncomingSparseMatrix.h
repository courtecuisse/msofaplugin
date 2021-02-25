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



template<class VecReal,class VecInt>
class ProjectIncomingMatrix : public defaulttype::BaseMatrix {
public:
    typedef typename CompressableMatrix<VecReal,VecInt>::Real Real;
    typedef defaulttype::BaseMatrix::Index Index;

    ProjectIncomingMatrix(helper::vector<int> & clearCols, helper::vector<int> & clearRows)
    : m_clearCols(clearCols)
    , m_clearRows(clearRows)
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

        m_P_eigen.resize(m_rowSize,m_colSize);
        m_P_eigen.setFromTriplets(m_setVal.begin(),m_setVal.end());

        EigenCompressedMatrix<Real> P(m_P_eigen);

        //push into csr format
        m_colptr.clear();
        m_rowind.clear();
        m_values.clear();
        const int * colptr = P.getColptr();
        const int * rowind = P.getRowind();
        const Real * values = P.getValues();
        m_rows = P.rowSize();
        m_cols = P.colSize();
        for(int i=0; i<m_rows+1; i++){
            m_colptr.push_back(colptr[i]);
        }
        m_nnz = m_colptr[m_rows];
        for(unsigned i=0; i<m_nnz; i++){
            m_rowind.push_back(rowind[i]);
            m_values.push_back(values[i]);
        }

//        for(int i=0; i<m_rows; i++){
//            for(int j=m_colptr[i]; j<m_colptr[i+1]; j++){
//                std::cout<<"row = "<<i<<"  col = "<<rowind[j]<<"  val = "<<values[j]<<std::endl;
//            }
//        }
//
//        std::cout<<"nnz = "<<nnz<<std::endl;

    }

    Eigen::SparseMatrix<Real,Eigen::RowMajor> & getEigenCSRMatrix() {
        return m_P_eigen;
    }

    const int rowSize() {
        return m_rows;
    }

    const int colSize() {
        return m_cols;
    }

    const int getNnz() {
        return m_nnz;
    }

    const VecInt & getColptr() {
        return m_colptr;
    }

    const VecInt & getRowind() {
        return m_rowind;
    }

    const VecReal & getValues() {
        return m_values;
    }

private :
    unsigned m_rowSize,m_colSize;
    helper::vector<int> & m_clearCols;
    helper::vector<int> & m_clearRows;
    helper::vector<Eigen::Triplet<Real>> m_setVal;
    Eigen::SparseMatrix<Real,Eigen::RowMajor> m_P_eigen;

    int m_rows;
    int m_cols;
    int m_nnz;
    VecInt m_colptr;
    VecInt m_rowind;
    VecReal m_values;
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

        if (m_projectionMatrix.colSize() != this->m_globalSize ||
            m_projectionMatrix.rowSize() != this->m_globalSize) {

            m_projectionMatrix.clear();
            m_projectionMatrix.resize(this->m_globalSize,this->m_globalSize);
            m_projectionMatrix.buildMatrix(this->m_stateAccessor,this->getContext());
        }

        for (unsigned i=0;i<m_matrices.size();i++) {
            if (m_matrices[i]->colSize() != this->m_globalSize ||
                m_matrices[i]->rowSize() != this->m_globalSize) {

                m_matrices[i]->clear();
                m_matrices[i]->resize(this->m_globalSize,this->m_globalSize);
                m_matrices[i]->buildMatrix(this->m_stateAccessor,m_clearCols,m_clearRows);
            } else {
                m_matrices[i]->fastReBuild();
            }
        }
    }

    typename CompressedMatrix<Real>::SPtr getCompressedMatrix(const core::MechanicalParams * mparams) {
        M_global = m_projectionMatrix.getEigenCSRMatrix();

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


    virtual std::string getTemplateName() const override {
        return templateName(this);
    }

    static std::string templateName(const BaseIncomingSparseMatrix<VecReal,VecInt>* = NULL) {
        return realName(VecReal().data());
    }

protected:
    unsigned m_rowSize,m_colSize;
    std::vector<typename LocalIncomingSparseMatrix<VecReal,VecInt>::SPtr> m_matrices;
    Eigen::SparseMatrix<Real,Eigen::RowMajor> M_global;
    helper::vector<int> m_clearCols,m_clearRows;
    ProjectIncomingMatrix<VecReal, VecInt> m_projectionMatrix;

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

//        double dot = this->dot(x, x);
//        std::cout<<"dot_x = "<<dot<<std::endl;
//        dot = this->dot(b, b);
//        std::cout<<"dot_b = "<<dot<<std::endl;

        TReal * x_ptr = X->data();
        const TReal * b_ptr = B->data();

        Eigen::Map<typename IncomingEigenType<TReal>::Vector> ex(x_ptr,this->m_globalSize);
        Eigen::Map<typename IncomingEigenType<TReal>::Vector> eb((TReal*) b_ptr,this->m_globalSize);

        if (acc) ex += this->m_projectionMatrix.getEigenCSRMatrix() * eb;
        else ex = this->m_projectionMatrix.getEigenCSRMatrix() * eb;

        for (unsigned i=0;i<this->m_matrices.size();i++) {
            Eigen::Map< Eigen::SparseMatrix<Real,Eigen::RowMajor> > M(this->m_matrices[i]->colSize(),
                                                                      this->m_matrices[i]->rowSize(),
                                                                      this->m_matrices[i]->getNnz(),
                                                                      this->m_matrices[i]->getColptr().data(),
                                                                      this->m_matrices[i]->getRowind().data(),
                                                                      this->m_matrices[i]->getValues().data());

            ex += M * eb * this->m_matrices[i]->getFactor(mparams);
        }

        double dot_res = this->dot(x, x);
        std::cout<<"dot_res = "<<dot_res<<std::endl;
    }
};


}
