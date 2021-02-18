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

#include <sofa/helper/AdvancedTimer.h>

namespace sofa::msofaplugin::matrix {

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
        LocalMIncomingSparseMatrix<VecReal,VecInt>::doCreateVisitor(m_matrices,this->getContext());
        LocalKIncomingSparseMatrix<VecReal,VecInt>::doCreateVisitor(k_matrices,this->getContext());

        for (unsigned i=0;i<m_matrices.size();i++) {
            m_matrices[i]->resize(this->m_globalSize,this->m_globalSize);
            m_matrices[i]->buildMatrix(this->m_stateAccessor);
        }

        for (unsigned i=0;i<k_matrices.size();i++) {
            m_matrices[i]->resize(this->m_globalSize,this->m_globalSize);
            k_matrices[i]->buildMatrix(this->m_stateAccessor);
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
        /*int id = addParams(mparams);
        return m_compressedValues[id];*/
    }

    virtual unsigned colSize() const {
        return m_colSize;
    }

    virtual unsigned rowSize() const {
        return m_rowSize;
    }

protected:
    unsigned m_rowSize,m_colSize;
    std::vector<typename LocalMIncomingSparseMatrix<VecReal,VecInt>::SPtr> m_matrices;
    std::vector<typename LocalKIncomingSparseMatrix<VecReal,VecInt>::SPtr> k_matrices;
};

template<class TReal>
class IncomingSparseMatrix : public BaseIncomingSparseMatrix<helper::vector<TReal>, helper::vector<int>> {
public:


    typedef TReal Real;
    typedef helper::vector<TReal> VecReal;
    typedef helper::vector<int> VecInt;

    SOFA_CLASS(SOFA_TEMPLATE(IncomingSparseMatrix,TReal),SOFA_TEMPLATE2(BaseIncomingSparseMatrix, VecReal, VecInt));

    Data<int> d_thread;

    IncomingSparseMatrix()
    : d_thread(initData(&d_thread, 4,"thread","Number of threads")) {}
/*
    inline void compress() const override {
        if (this->m_toCompress.empty()) return;

        int nnz = this->m_matrix.getColPtr()[this->m_matrix.rowSize()];

        for (unsigned k=0;k<this->m_toCompress.size();k++) {
            auto & ic = this->m_compressedValues[this->m_toCompress[k]];
            ic->m_values.clear();
            ic->m_values.resize(nnz);
        }

        const TReal * sparsevalues = this->m_matrix.getSparseValues().data();
        const int * valptr = this->m_matrix.getValPtr().data();
        const int * sortid = this->m_matrix.getsortId().data();

        auto thread_worker = [&] (int start, int end) {
            for (unsigned i=start;i<end;i++) {
                for (unsigned k=0;k<this->m_toCompress.size();k++) {
                    auto & ic = this->m_compressedValues[this->m_toCompress[k]];
                    for (unsigned j=valptr[i];j<valptr[i+1];j++) {
                        const int id = sortid[j];
                        const double fact = (id < this->m_kid) ? ic->kfact : (id < this->m_mid) ? ic->mfact : 1.0;
                        ic->m_values[i] += sparsevalues[id] * fact;
                    }
                }
            }
        };

        int start = 0;
        int NBTHREAD = d_thread.getValue();
        int NBLOCS=(nnz+NBTHREAD)/NBTHREAD;
        std::vector<std::thread> threads;
        threads.resize(NBTHREAD);
        for (int t=0;t<NBTHREAD;t++) {
            int end = std::min(start+NBLOCS,nnz);
            threads[t]=std::thread(thread_worker, start, end);
            start=end;
        }
        for (int t=0;t<NBTHREAD;t++) threads[t].join();
        this->m_toCompress.clear();
    }
*/
    void mult(core::MechanicalParams * mparams, BaseSystemMatrix::MechanicalVectorId & x, const BaseSystemMatrix::MechanicalVectorId & b, bool acc) override {

        std::cout << "MULT " << std::endl;
        /*
//        TIMER_START(Mult);
        unsigned id = this->addParams(mparams);

        auto X = this->getMechanicalVector(x);
        auto B = this->getMechanicalVector(b);

        TReal * x_ptr = X->write()->data();
        const TReal * b_ptr = B->read()->data();

        if (! acc) memset(x_ptr,0,this->m_matrix.rowSize()*sizeof(TReal));

        const int * colptr = this->m_matrix.getColPtr().data();
        const int * rowind = this->m_matrix.getRowInd().data();
        const Real * values = this->m_compressedValues[id]->m_values.data();

        auto thread_worker = [&] (int start, int end) {
            for (unsigned j=start;j<end;j++) {
                for (unsigned i=colptr[j];i<colptr[j+1];i++) {
                    x_ptr[j] += b_ptr[rowind[i]] * values[i];
                }
            }
        };

        const int size = this->m_matrix.rowSize();
        int NBTHREAD = d_thread.getValue();
        int NBLOCS=(size+NBTHREAD-1)/NBTHREAD;
        int start = 0;
        std::vector<std::thread> threads(NBTHREAD);
        for (int t=0;t<NBTHREAD;t++) {
            int end = std::min(start+NBLOCS,size);
            threads[t]=std::thread(thread_worker, start, end);
            start=end;
        }
        for (int t=0;t<NBTHREAD;t++) threads[t].join();


//        TIMER_END(Mult);
//        TIMER_PRINT("mult CG " << Mult << " ms");
*/
    }
};

//template<> std::string IncomingSparseMatrix<float>::Name() { return "float"; }
//template<> std::string IncomingSparseMatrix<double>::Name() { return "double"; }

}
