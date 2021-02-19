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

    const int * getRowBegin() const { return m_matrix.outerIndexPtr(); }

    const int * getColsIndex() const { return m_matrix.innerIndexPtr(); }

    const Real * getColsValue() const { return m_matrix.valuePtr(); }

    virtual unsigned colSize() const { return m_matrix.cols(); }

    virtual unsigned rowSize() const { return m_matrix.rows(); }

    Eigen::SparseMatrix<Real,Eigen::RowMajor> & m_matrix;
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

    virtual void clear() {
    }

    void buildMatrix() override {
        LocalMIncomingSparseMatrix<VecReal,VecInt>::doCreateVisitor(m_Matrices,this->getContext());
        LocalKIncomingSparseMatrix<VecReal,VecInt>::doCreateVisitor(m_Katrices,this->getContext());

        for (unsigned i=0;i<m_Matrices.size();i++) {
            m_Matrices[i]->resize(this->m_globalSize,this->m_globalSize);
            m_Matrices[i]->buildMatrix(this->m_stateAccessor);
        }

        for (unsigned i=0;i<m_Katrices.size();i++) {
            m_Katrices[i]->resize(this->m_globalSize,this->m_globalSize);
            m_Katrices[i]->buildMatrix(this->m_stateAccessor);
        }

//        sofa::helper::AdvancedTimer::stepNext ("collect", "project");

//        ProjectMatrixVisitor(&accessor).execute(this->getContext());

//        sofa::helper::AdvancedTimer::stepNext ("project", "pattern");

//        m_matrix.rebuildPatternAccess();
//        sofa::helper::AdvancedTimer::stepNext ("pattern", "compress");
//        for (unsigned i=0;i<m_compressedValues.size();i++) m_toCompress.push_back(i);
//        compress();

//        sofa::helper::AdvancedTimer::stepEnd("compress");*/


    }

    typename CompressedMatrix<Real>::SPtr getCompressedMatrix(const core::MechanicalParams * mparams) {
        M_global.resize(0,0);
        M_global.resize(this->m_globalSize,this->m_globalSize);

        for (unsigned i=0;i<m_Matrices.size();i++) {
            Eigen::Map< Eigen::SparseMatrix<Real,Eigen::RowMajor> > map(m_Matrices[i]->colSize(),
                                                                        m_Matrices[i]->rowSize(),
                                                                        m_Matrices[i]->getNnz(),
                                                                        m_Matrices[i]->getColptr().data(),
                                                                        m_Matrices[i]->getRowind().data(),
                                                                        m_Matrices[i]->getValues().data());

            M_global += map*mparams->mFactor();
        }

        for (unsigned i=0;i<m_Katrices.size();i++) {
            Eigen::Map< Eigen::SparseMatrix<Real,Eigen::RowMajor> > map(m_Katrices[i]->colSize(),
                                                                        m_Katrices[i]->rowSize(),
                                                                        m_Katrices[i]->getNnz(),
                                                                        m_Katrices[i]->getColptr().data(),
                                                                        m_Katrices[i]->getRowind().data(),
                                                                        m_Katrices[i]->getValues().data());

            M_global += map*mparams->kFactor();
        }

        return typename CompressedMatrix<Real>::SPtr(new EigenCompressedMatrix<Real>(M_global));
    }

    virtual unsigned colSize() const {
        return m_colSize;
    }

    virtual unsigned rowSize() const {
        return m_rowSize;
    }

protected:
    unsigned m_rowSize,m_colSize;
    std::vector<typename LocalMIncomingSparseMatrix<VecReal,VecInt>::SPtr> m_Matrices;
    std::vector<typename LocalKIncomingSparseMatrix<VecReal,VecInt>::SPtr> m_Katrices;
    Eigen::SparseMatrix<Real,Eigen::RowMajor> M_global;
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

        if (! acc) {
            for (unsigned i=0;i<this->m_globalSize;i++) x_ptr[i] = 0;
        }

        Eigen::Map<typename IncomingEigenType<TReal>::Vector> ex(x_ptr,this->m_globalSize);
        Eigen::Map<typename IncomingEigenType<TReal>::Vector> eb((TReal*) b_ptr,this->m_globalSize);

        for (unsigned i=0;i<this->m_Matrices.size();i++) {
            Eigen::Map< Eigen::SparseMatrix<Real,Eigen::RowMajor> > M(this->m_Matrices[i]->colSize(),
                                                                      this->m_Matrices[i]->rowSize(),
                                                                      this->m_Matrices[i]->getNnz(),
                                                                      this->m_Matrices[i]->getColptr().data(),
                                                                      this->m_Matrices[i]->getRowind().data(),
                                                                      this->m_Matrices[i]->getValues().data());

            ex += M * eb * mparams->mFactor();
        }

        for (unsigned i=0;i<this->m_Katrices.size();i++) {
            Eigen::Map< Eigen::SparseMatrix<Real,Eigen::RowMajor> > K(this->m_Katrices[i]->colSize(),
                                                                      this->m_Katrices[i]->rowSize(),
                                                                      this->m_Katrices[i]->getNnz(),
                                                                      this->m_Katrices[i]->getColptr().data(),
                                                                      this->m_Katrices[i]->getRowind().data(),
                                                                      this->m_Katrices[i]->getValues().data());

            ex += K * eb * mparams->kFactor();
        }
    }
};


}
