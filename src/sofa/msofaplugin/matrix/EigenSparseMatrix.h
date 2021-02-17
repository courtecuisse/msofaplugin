#pragma once

#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/msofaplugin/matrix/CompressableMatrix.h>
#include <sofa/msofaplugin/matrix/BuildMatrixVisitors.h>
#include <Eigen/Sparse>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/MechanicalMatrixVisitor.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sofa/helper/AdvancedTimer.h>

namespace sofa::msofaplugin::matrix {

template<class Real>
struct EigenType;

template<> struct EigenType<double> {
    typedef Eigen::VectorXd Vector;
};

template<> struct EigenType<float> {
    typedef Eigen::VectorXf Vector;
};

template<class TReal>
class EigenBaseMatrix : public defaulttype::BaseMatrix {
public:
    typedef TReal Real;
    typedef defaulttype::BaseMatrix::Index Index;

    Index rowSize(void) const override {
        return m_rowSize;
    }

    Index colSize(void) const override {
        return m_colSize;
    }

    SReal element(Index i, Index j) const override {
        compress();
        return m_compressed.coeffRef(i,j);
    }

    void resize(Index nbRow, Index nbCol) {
        m_rowSize = nbRow;
        m_colSize = nbCol;
    }

    void clear() {
        m_rowSize = 0;
        m_colSize = 0;
        m_incoming.clear();
        m_compressed.resize(0,0);
    }

    inline const Eigen::SparseMatrix<Real,Eigen::RowMajor> & compress() const {
        if (m_incoming.empty()) {
            if (m_compressed.size() == 0) m_compressed.resize(m_rowSize,m_colSize);
            return m_compressed;
        }

        if (m_compressed.size() == 0) {
            m_compressed.resize(m_rowSize,m_colSize);
            m_compressed.setFromTriplets(m_incoming.begin(),m_incoming.end());
        } else {
            Eigen::SparseMatrix<Real,Eigen::RowMajor> tmp;
            tmp.resize(m_rowSize,m_colSize);
            tmp.setFromTriplets(m_incoming.begin(),m_incoming.end());
            m_compressed += tmp;
        }
        m_incoming.clear();

        return m_compressed;
    }

    std::vector<Eigen::Triplet<Real> > & getIncoming() {
        return m_incoming;
    }

protected :
    unsigned m_colSize,m_rowSize;
    mutable Eigen::SparseMatrix<Real,Eigen::RowMajor> m_compressed;
    mutable std::vector<Eigen::Triplet<Real> > m_incoming;
};

template<class Real>
class IncomingEigenMatrix : public EigenBaseMatrix<Real> {
public:
    typedef defaulttype::BaseMatrix::Index Index;

    void set(Index i, Index j, double v) override {}

    void add( Index i, Index j, double v ) override {
        this->m_incoming.push_back(Eigen::Triplet<Real>(i,j,v));
    }

    inline void add(Index i, Index j, const defaulttype::Mat3x3d & M) override {
        this->m_incoming.push_back(Eigen::Triplet<Real>(i+0,j+0,M[0][0]));
        this->m_incoming.push_back(Eigen::Triplet<Real>(i+0,j+1,M[0][1]));
        this->m_incoming.push_back(Eigen::Triplet<Real>(i+0,j+2,M[0][2]));

        this->m_incoming.push_back(Eigen::Triplet<Real>(i+1,j+0,M[1][0]));
        this->m_incoming.push_back(Eigen::Triplet<Real>(i+1,j+1,M[1][1]));
        this->m_incoming.push_back(Eigen::Triplet<Real>(i+1,j+2,M[1][2]));

        this->m_incoming.push_back(Eigen::Triplet<Real>(i+2,j+0,M[2][0]));
        this->m_incoming.push_back(Eigen::Triplet<Real>(i+2,j+1,M[2][1]));
        this->m_incoming.push_back(Eigen::Triplet<Real>(i+2,j+2,M[2][2]));
    }

    void clearRow(Index i) override {
        this->compress();
        for (typename Eigen::SparseMatrix<Real,Eigen::RowMajor>::InnerIterator it(this->m_compressed,i); it; ++it)
            it.valueRef() = 0;
    }

    void clearCol(Index j) override {
        this->compress();
        for (typename Eigen::SparseMatrix<Real,Eigen::RowMajor>::InnerIterator it(this->m_compressed,j); it; ++it)
            this->m_compressed.coeffRef(it.col(),j) = 0;
    }
};


template<class Real>
class ProjectEigenMatrix : public EigenBaseMatrix<Real> {
public:
    typedef defaulttype::BaseMatrix::Index Index;

    void set(Index i, Index j, double v) override {
        this->m_incoming.push_back(Eigen::Triplet<Real>(i,j,v));
    }

    void add( Index i, Index j, double v ) override {}

    void clearRow(Index i) override {}

    void clearCol(Index j) override {}
};

template<class Real>
class EigenCompressedMatrix : public CompressedMatrix<Real> {
public:
    EigenCompressedMatrix(Eigen::SparseMatrix<Real,Eigen::RowMajor> * mat) : m_matrix(mat) {}

    const int * getRowBegin() const {
        return m_matrix->outerIndexPtr();
    }

    const int * getColsIndex() const {
        return m_matrix->innerIndexPtr();
    }

    const Real * getColsValue() const {
        return m_matrix->valuePtr();
    }

    virtual unsigned colSize() const {
        return m_matrix->cols();
    }

    virtual unsigned rowSize() const {
        return m_matrix->rows();
    }

    Eigen::SparseMatrix<Real,Eigen::RowMajor> * m_matrix;
};

template<class TReal>
class EigenSparseMatrix : public CompressableMatrix<helper::vector<TReal>,helper::vector<int> > {
public:
    typedef TReal Real;
    typedef helper::vector<Real> VecReal;
    typedef helper::vector<int> VecInt;

    SOFA_CLASS(SOFA_TEMPLATE(EigenSparseMatrix,TReal),SOFA_TEMPLATE2(CompressableMatrix, VecReal,VecInt));

    void clear() {
        m_M.clear();
        m_K.clear();
        m_P.clear();
    }

    void buildMatrix() override {
        sofa::helper::AdvancedTimer::stepBegin("prepare");

        component::linearsolver::DefaultMultiMatrixAccessor accessorM;
        component::linearsolver::DefaultMultiMatrixAccessor accessorK;
        component::linearsolver::DefaultMultiMatrixAccessor accessorP;

        accessorM.setGlobalMatrix(&m_M);
        accessorK.setGlobalMatrix(&m_K);
        accessorP.setGlobalMatrix(&m_P);

        accessorM.setupMatrices();
        accessorK.setupMatrices();
        accessorP.setupMatrices();

        for (unsigned i=0;i<this->m_vectorOperations.size();i++) {
            accessorM.addMechanicalState(this->m_vectorOperations[i]->getState());
            accessorK.addMechanicalState(this->m_vectorOperations[i]->getState());
            accessorP.addMechanicalState(this->m_vectorOperations[i]->getState());
        }
        unsigned size = this->m_globalSize;

        m_M.resize(size,size);
        m_K.resize(size,size);
        m_P.resize(size,size);

        sofa::helper::AdvancedTimer::stepNext ("prepare", "build");

        BuildMVisitor(&accessorM).execute(this->getContext());
        BuildKVisitor(&accessorK).execute(this->getContext());

        sofa::helper::AdvancedTimer::stepNext ("build", "project+compress");

        ProjectMatrixVisitor(&accessorM).execute(this->getContext());
        ProjectMatrixVisitor(&accessorK).execute(this->getContext());
        ProjectMatrixVisitor(&accessorP).execute(this->getContext());

        sofa::helper::AdvancedTimer::stepEnd("project+compress");

//        TIMER_PRINT("Prepare=" << prepare << " ms    Build=" << build << " ms project=" << project << " ms  total=" << (prepare+build+project) << " ms");;
    }

    void mult(core::MechanicalParams * mp, BaseSystemMatrix::MechanicalVectorId & x, const BaseSystemMatrix::MechanicalVectorId & b, bool acc) override {
//        TIMER_START(mult);

        mult(m_M.compress(),mp->mFactor(), x, b, acc);
        mult(m_K.compress(),mp->kFactor(), x, b);
        mult(m_P.compress(),1.0, x, b);

//        TIMER_END(mult);

//        TIMER_PRINT("mult =" << mult << " ms");

    }

//    virtual void mult(const core::MultiVecDerivId & x, const core::MultiVecDerivId & b,bool acc) override {
////        TIMER_START(multf);

//        mult(m_compressed,1.0,x,b,acc);

////        TIMER_END(multf);
////        TIMER_PRINT("multf=" << multf << " ms");
//    }

    virtual std::shared_ptr<CompressedMatrix<TReal>> getCompressedMatrix(const core::MechanicalParams * mparams) {
        m_compressed = m_M.compress()*mparams->mFactor() +
                       m_K.compress()*mparams->kFactor() +
                       m_P.compress();

//#if 1
//    std::ofstream myfile;
//    myfile.open ("example.txt");

//    myfile << "E=zeros("<<rowSize()<< ","<<colSize()<<");"<<std::endl;
//    for (unsigned j=0;j<rowSize();j++) {
//        for (unsigned i=getRowBegin()[j];i<getRowBegin()[j+1];i++) {
//            myfile << "E("<<j+1<< ","<<getColsIndex()[i]+1<<")="<<getColsValue()[i]<<";"<<std::endl;
//        }
//    }

//    myfile.close();
//    exit(1);
//#endif

        return std::shared_ptr<CompressedMatrix<TReal>>(new EigenCompressedMatrix<TReal>(&m_compressed));
    }

    virtual unsigned colSize() const {
        return m_M.cols();
    }

    virtual unsigned rowSize() const {
        return m_M.rows();
    }

//    static std::string Name();

//    virtual std::string getTemplateName() const override {
//        return Name();
//    }

protected:
    ProjectEigenMatrix<TReal> m_P;
    IncomingEigenMatrix<TReal> m_M,m_K;
    Eigen::SparseMatrix<TReal,Eigen::RowMajor> m_compressed;

    inline void mult(const Eigen::SparseMatrix<TReal,Eigen::RowMajor> & M, double fact, BaseSystemMatrix::MechanicalVectorId & x, const BaseSystemMatrix::MechanicalVectorId & b, bool acc = true) {
        auto X = this->getMechanicalVector(x)->write();
        auto B = this->getMechanicalVector(b)->read();

        TReal * x_ptr = X->data();
        const TReal * b_ptr = B->data();

        Eigen::Map<typename EigenType<TReal>::Vector> ex(x_ptr,rowSize());
        Eigen::Map<typename EigenType<TReal>::Vector> eb((TReal*) b_ptr,colSize());

        if (acc) ex += M * eb * fact ;
        else ex = M * eb * fact;
    }


};

//template<> std::string EigenSparseMatrix<float>::Name() { return "float"; }
//template<> std::string EigenSparseMatrix<double>::Name() { return "double"; }

}
