#pragma once

#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/msofaplugin/matrix/CompressableMatrix.h>
#include <sofa/msofaplugin/matrix/BuildMatrixVisitors.h>
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
class IncomingBaseMatrix : public defaulttype::BaseMatrix {
public:
    typedef typename RealType<VecReal>::Real Real;
    typedef defaulttype::BaseMatrix::Index Index;

    struct MatrixCoord {
        MatrixCoord() {}
        MatrixCoord(unsigned r,unsigned c) : row(r), col(c) {}
        unsigned row;
        unsigned col;
    };

    IncomingBaseMatrix() {
        m_keepStruct = false;
        m_compressed = false;
    }

    Index rowSize(void) const override {
        return m_rowSize;
    }

    Index colSize(void) const override {
        return m_colSize;
    }

    Index nnz(void) const {
        return m_nnz;
    }

    SReal element(Index i, Index j) const override { return 0.0; }

    void resize(Index nbRow, Index nbCol) {
        m_rowSize = nbRow;
        m_colSize = nbCol;
    }

    void clear() {
        m_rowSize = 0;
        m_colSize = 0;
        m_writeId = 0;
        m_compressed = false;
        m_clearCols.clear();
        m_clearRows.clear();
        m_setValId.clear();
        m_sparseValuesPtr = m_sparseValuesVec.data();
    }

    void set(Index i, Index j, double v) override {
        m_setValId.push_back(m_writeId);
        add(i,j,v);
    }

    inline void add( Index i, Index j, double v ) override {
        m_keepStruct = (m_keepStruct) &&
                       (m_sparseIndices.size() > m_writeId) &&
                       (m_sparseIndices[m_writeId].row == i) && //not necessary?
                       (m_sparseIndices[m_writeId].col == j);

        if (! m_keepStruct) {
            m_sparseValuesVec.resize(m_writeId+1);
            m_sparseValuesPtr = m_sparseValuesVec.data();

            m_sparseIndices.resize(m_writeId+1);
            m_sparseIndices[m_writeId] = MatrixCoord(i,j);
        }

        m_sparseValuesPtr[m_writeId] = v;
        m_writeId++;
    }

    inline void add(Index i, Index j, const defaulttype::Mat3x3d & M) override {
        // only check for for topleft coord
        m_keepStruct = m_keepStruct &&
                       m_sparseIndices.size() > (m_writeId+9) &&
                       m_sparseIndices[m_writeId].row == i && //not necessary?
                       m_sparseIndices[m_writeId].col == j;

        if (! m_keepStruct) {
            m_sparseValuesVec.fastResize(m_writeId+9);
            m_sparseValuesPtr = m_sparseValuesVec.data();

            m_sparseIndices.resize(m_writeId+9);
            m_sparseIndices[m_writeId+0] = MatrixCoord(i+0,j+0);
            m_sparseIndices[m_writeId+1] = MatrixCoord(i+0,j+1);
            m_sparseIndices[m_writeId+2] = MatrixCoord(i+0,j+2);

            m_sparseIndices[m_writeId+3] = MatrixCoord(i+1,j+0);
            m_sparseIndices[m_writeId+4] = MatrixCoord(i+1,j+1);
            m_sparseIndices[m_writeId+5] = MatrixCoord(i+1,j+2);

            m_sparseIndices[m_writeId+6] = MatrixCoord(i+2,j+0);
            m_sparseIndices[m_writeId+7] = MatrixCoord(i+2,j+1);
            m_sparseIndices[m_writeId+8] = MatrixCoord(i+2,j+2);
        }

        m_sparseValuesPtr[m_writeId++] = M[0][0];
        m_sparseValuesPtr[m_writeId++] = M[0][1];
        m_sparseValuesPtr[m_writeId++] = M[0][2];

        m_sparseValuesPtr[m_writeId++] = M[1][0];
        m_sparseValuesPtr[m_writeId++] = M[1][1];
        m_sparseValuesPtr[m_writeId++] = M[1][2];

        m_sparseValuesPtr[m_writeId++] = M[2][0];
        m_sparseValuesPtr[m_writeId++] = M[2][1];
        m_sparseValuesPtr[m_writeId++] = M[2][2];
    }

    void clearRow(Index i) override {
        m_clearRows.push_back(i);
    }

    void clearCol(Index j) override {
        m_clearCols.push_back(j);
    }

    void clearRowCol(Index i) override {
        clearCol(i); // first delete the comun as it uses row values
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

    void inline rebuildPatternAccess() {
        if (m_keepStruct) return;

        std::vector<int> countvec;

        std::vector<int> tmp_colptr;//CRS format
        std::vector<int> tmp_rowind;

        std::vector<int> tmp_tran_colptr;//Transposed uncollapsed CRS format
        std::vector<int> tmp_tran_rowind;
        std::vector<int> tmp_tran_sortid;//index of the incoming values in the transposed format

        unsigned tmp_nnz = m_writeId;

        tmp_tran_colptr.resize(m_colSize+1);
        tmp_tran_rowind.resize(tmp_nnz);
        tmp_tran_sortid.resize(tmp_nnz);

        tmp_colptr.resize(m_rowSize+1);
        tmp_rowind.resize(tmp_nnz);
        m_sortid.resize(tmp_nnz);


        //First we build the uncollapsed transposed matrix
        countvec.resize(m_colSize);
        for (unsigned i=0;i<tmp_nnz;i++) {
            unsigned col = m_sparseIndices[i].col;
            countvec[col]++;
        }

        //clear the columns set the count values to 0
        for (unsigned i=0;i<m_clearCols.size();i++) countvec[m_clearCols[i]] = 0;

        tmp_tran_colptr[0] = 0;
        for (int j=0;j<m_colSize;j++) tmp_tran_colptr[j+1] = tmp_tran_colptr[j] + countvec[j];

        countvec.clear();
        countvec.resize(m_colSize); // reset to zero

        for (unsigned i=0;i<tmp_nnz;i++) {
            unsigned col = m_sparseIndices[i].col;

            unsigned write_id = tmp_tran_colptr[col] + countvec[col];
            if (write_id>=tmp_tran_colptr[col+1]) continue; //skip clear columns

            unsigned row = m_sparseIndices[i].row;

            tmp_tran_rowind[write_id] = row;
            tmp_tran_sortid[write_id] = i;

            countvec[col]++;
        }

        //Compute the transpose ==> provides implicit sorting

        countvec.clear();
        countvec.resize(m_rowSize); // reset to zero

        //Compute the transpose count number of values per lines
        for (unsigned j=0;j<m_colSize;j++) {
            for (unsigned i=tmp_tran_colptr[j];i<tmp_tran_colptr[j+1];i++) {
                int row = tmp_tran_rowind[i];
                countvec[row]++;
            }
        }

        //clear the columns set the count values to 0
        for (unsigned i=0;i<m_clearRows.size();i++) countvec[m_clearRows[i]] = 0;

        for (unsigned i=0;i<m_setValId.size();i++) {
            int id=m_setValId[i];
            int row = m_sparseIndices[id].row;
            countvec[row] = 1; // only one value can be set per line
        }

        tmp_colptr[0] = 0;
        for (int j=0;j<m_rowSize;j++) tmp_colptr[j+1] = tmp_colptr[j] + countvec[j];

        countvec.clear();
        countvec.resize(m_rowSize); // reset to zero

        for (unsigned j=0;j<m_colSize;j++) {
            for (unsigned i=tmp_tran_colptr[j];i<tmp_tran_colptr[j+1];i++) {
                unsigned row = tmp_tran_rowind[i];

                unsigned write_id = tmp_colptr[row] + countvec[row];
                if (write_id>=tmp_colptr[row+1]) continue; //skip clear rows

                tmp_rowind[write_id] = j; // we are building the column j
                m_sortid[write_id] = tmp_tran_sortid[i]; // build the sortidvector

                countvec[row]++;
            }
        }

        for (unsigned i=0;i<m_setValId.size();i++) {
            int id=m_setValId[i];
            int row = m_sparseIndices[id].row;
            tmp_rowind[tmp_colptr[row]] = m_sparseIndices[id].col;
            m_sortid[tmp_colptr[row]] = id;
        }

        m_valptr.clear();
        m_rowind.clear();
        m_colptr.clear();
        m_colptr.push_back(0);

        //Compute the compressed row sparse format and collapse values
        for (unsigned j=0;j<m_rowSize;j++) {
            int newcol=-1;
            for (unsigned i=tmp_colptr[j];i<tmp_colptr[j+1];i++) {
                // add new values only if needed
                if (tmp_rowind[i] != newcol) {
                    newcol = tmp_rowind[i];
                    m_rowind.push_back(newcol);
                    m_valptr.push_back(i);
                }
            }

            m_colptr.push_back(m_rowind.size());
        }

        m_valptr.push_back(tmp_colptr[m_rowSize]);

        m_nnz = m_colptr[m_rowSize];

        m_keepStruct = true;
    }

    const VecInt & getColPtr() const {
        return m_colptr;
    }

    const VecInt & getRowInd() const {
        return m_rowind;
    }

    const VecInt & getValPtr() const {
        return m_valptr;
    }

    const VecInt & getsortId() const {
        return m_sortid;
    }

    const VecReal & getSparseValues() const {
        return m_sparseValuesVec;
    }

    unsigned getWriteId() {
        return m_writeId;
    }

    int getSize(){
        return m_sparseIndices.size();
    }

private :
    unsigned m_colSize,m_rowSize,m_nnz;

    Real * m_sparseValuesPtr; // this ptr is used for gpu data It allows to not check the avalaibility of data every tim a new incoming value is pused in the vector
    mutable VecReal m_sparseValuesVec; // incoming new values
    mutable helper::vector<MatrixCoord> m_sparseIndices; // incoming vector of coordinated used to check the consistency of the structure
    mutable bool m_keepStruct;
    mutable bool m_compressed;
    mutable int m_writeId;

    mutable VecInt m_colptr;//CRS format
    mutable VecInt m_rowind;

    mutable VecInt m_valptr;// start and end point of each sum values in the CSR m_values vector
    mutable VecInt m_sortid;// corresponding indices in the m_sparseIndices vector

    mutable helper::vector<int> m_clearCols,m_clearRows,m_setValId;
};

template<class Real, class VecReal, class VecInt>
class InternalCompressed : public CompressedMatrix<Real> {
public:

    InternalCompressed(const IncomingBaseMatrix<VecReal,VecInt> * matrix, double m, double k)
    : m_matrix(matrix)
    , mfact(m), kfact(k), used(true) {}

    const int * getRowBegin() const override {
        return m_matrix->getColPtr().data();
    }

    const int * getColsIndex() const override {
        return m_matrix->getRowInd().data();
    }

    const Real * getColsValue() const override {
        return m_values.data();
    }

    virtual unsigned colSize() const override {
        return m_matrix->colSize();
    }

    virtual unsigned rowSize() const override {
        return m_matrix->rowSize();
    }

    const IncomingBaseMatrix<VecReal,VecInt> * m_matrix;
    VecReal m_values;
    double mfact;
    double kfact;
    bool used;
};

template<class VecReal,class VecInt>
class BaseIncomingSparseMatrix : public CompressableMatrix<VecReal,VecInt> {
public:

    typedef typename CompressableMatrix<VecReal,VecInt>::Real Real;
    typedef std::shared_ptr<InternalCompressed<Real, VecReal,VecInt>> FullCompressedMatrix;

    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE2(BaseIncomingSparseMatrix,VecReal,VecInt),SOFA_TEMPLATE2(CompressableMatrix, VecReal,VecInt) );

    virtual void clear() {
        m_matrix.clear();

        //Clear unused parameters
        for (int i=m_compressedValues.size()-1;i>=0;i--) {
            if (!m_compressedValues[i]->used) m_compressedValues.erase(m_compressedValues.begin()+i);
        }
    }

    void buildMatrix() override {
        sofa::helper::AdvancedTimer::stepBegin("prepare");
        TIMER_START(prepare);

        component::linearsolver::DefaultMultiMatrixAccessor accessor;
        accessor.setGlobalMatrix(&m_matrix);
        accessor.setupMatrices();

        for (unsigned i=0;i<this->m_vectorOperations.size();i++)
            accessor.addMechanicalState(this->m_vectorOperations[i]->getState());

        unsigned size = this->m_globalSize;
        m_matrix.resize(size,size);

        sofa::helper::AdvancedTimer::stepNext ("prepare", "collect");
        TIMER_NEXT(prepare,build);

        BuildKVisitor(&accessor).execute(this->getContext());
        m_kid = m_matrix.getWriteId();
        BuildMVisitor(&accessor).execute(this->getContext());
        m_mid = m_matrix.getWriteId();

        sofa::helper::AdvancedTimer::stepNext ("collect", "project");
        TIMER_NEXT(build,project);

        ProjectMatrixVisitor(&accessor).execute(this->getContext());

        sofa::helper::AdvancedTimer::stepNext ("project", "pattern");
        TIMER_NEXT(project,pattern);

        m_matrix.rebuildPatternAccess();
        sofa::helper::AdvancedTimer::stepNext ("pattern", "compress");
        for (unsigned i=0;i<m_compressedValues.size();i++) m_toCompress.push_back(i);
        compress();

        sofa::helper::AdvancedTimer::stepEnd("compress");
        TIMER_END(pattern);

//        TIMER_PRINT("Prepare=" << prepare << " ms    Build=" << build << " ms project=" << project << " ms pattern=" << pattern << " ms  total=" << (prepare+build+project+pattern) << " ms");;
    }

    std::shared_ptr<CompressedMatrix<Real>> getCompressedMatrix(const core::MechanicalParams * mparams) {
        int id = addParams(mparams);
        return m_compressedValues[id];
    }

//    static std::string Name();

//    virtual std::string getTemplateName() const override {
//        return Name();
//    }

    virtual unsigned colSize() const {
        return m_matrix.cols();
    }

    virtual unsigned rowSize() const {
        return m_matrix.rows();
    }

protected:
    IncomingBaseMatrix<VecReal,VecInt> m_matrix;
    unsigned m_mid,m_kid; // last index of mass and stiffness contributions
    mutable std::vector<FullCompressedMatrix> m_compressedValues; // memorize parameter in order to group compression step
    mutable std::vector<unsigned> m_toCompress;    

    unsigned addParams(const core::MechanicalParams * mparams) {
        for (unsigned i=0;i<m_compressedValues.size();i++) {
            if (m_compressedValues[i]->mfact == mparams->mFactor() &&
                m_compressedValues[i]->kfact == mparams->kFactor()) {
                m_compressedValues[i]->used = true;
                return i;
            }
        }

        unsigned id = m_compressedValues.size();
        m_compressedValues.push_back(FullCompressedMatrix(new InternalCompressed<Real, VecReal,VecInt>(&m_matrix, mparams->mFactor(), mparams->kFactor())));
        m_toCompress.push_back(id);
        compress();
        return id;
    }

    virtual void compress() const = 0;

};

template<class TReal>
class IncomingSparseMatrix : public BaseIncomingSparseMatrix<helper::vector<TReal>, helper::vector<int>> {
public:

    typedef TReal Real;
    typedef helper::vector<TReal> VecReal;
    typedef helper::vector<int> VecInt;

    SOFA_CLASS(SOFA_TEMPLATE(IncomingSparseMatrix,TReal),SOFA_TEMPLATE2(BaseIncomingSparseMatrix, VecReal, VecInt));

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
        int NBTHREAD = atoi(getenv("MSOFA_THREADS"));
        if (NBTHREAD<1) NBTHREAD=1;
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

    void mult(core::MechanicalParams * mparams, BaseSystemMatrix::MechanicalVectorId & x, const BaseSystemMatrix::MechanicalVectorId & b, bool acc) override {
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
        int NBTHREAD = atoi(getenv("MSOFA_THREADS"));
        if (NBTHREAD<1) NBTHREAD=1;
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
    }
};

//template<> std::string IncomingSparseMatrix<float>::Name() { return "float"; }
//template<> std::string IncomingSparseMatrix<double>::Name() { return "double"; }

}
