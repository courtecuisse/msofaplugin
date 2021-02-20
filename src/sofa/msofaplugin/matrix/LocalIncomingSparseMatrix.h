#pragma once

#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/msofaplugin/matrix/CompressableMatrix.h>
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
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <Eigen/Sparse>

namespace sofa::msofaplugin::matrix {

template<class VecReal,class VecInt>
class BuildingIncomingMatrix : public defaulttype::BaseMatrix {
public:
    typedef typename RealType<VecReal>::Real Real;
    typedef defaulttype::BaseMatrix::Index Index;

    struct MatrixCoord {
        MatrixCoord() {}
        MatrixCoord(unsigned r,unsigned c) : row(r), col(c) {}
        unsigned row;
        unsigned col;
    };

    BuildingIncomingMatrix(unsigned & r,unsigned & c,VecInt & colptr,VecInt & rowind,VecReal & values,VecInt & mapping)
    : m_rowSize(r)
    , m_colSize(c)
    , m_colptr(colptr)
    , m_rowind(rowind)
    , m_values(values)
    , m_mapping(mapping) {}

    Index rowSize(void) const override { return m_rowSize; }

    Index colSize(void) const override { return m_colSize; }

    SReal element(Index /*i*/, Index /*j*/) const override { return 0.0; }

    void resize(Index nbRow, Index nbCol) {
        m_rowSize = nbRow;
        m_colSize = nbCol;
    }

    void clear() {
        m_sparseIndices.clear();
        m_sparseValuesVec.clear();
        m_clearCols.clear();
        m_clearRows.clear();
    }

    void set(Index , Index , double ) override {}

    inline void add( Index i, Index j, double v ) override {
        m_sparseIndices.push_back(MatrixCoord(i,j));
        m_sparseValuesVec.push_back(v);
    }

    inline void add(Index i, Index j, const defaulttype::Mat3x3d & M) override {
        m_sparseIndices.push_back(MatrixCoord(i+0,j+0));
        m_sparseIndices.push_back(MatrixCoord(i+0,j+1));
        m_sparseIndices.push_back(MatrixCoord(i+0,j+2));

        m_sparseIndices.push_back(MatrixCoord(i+1,j+0));
        m_sparseIndices.push_back(MatrixCoord(i+1,j+1));
        m_sparseIndices.push_back(MatrixCoord(i+1,j+2));

        m_sparseIndices.push_back(MatrixCoord(i+2,j+0));
        m_sparseIndices.push_back(MatrixCoord(i+2,j+1));
        m_sparseIndices.push_back(MatrixCoord(i+2,j+2));

        m_sparseValuesVec.push_back(M[0][0]);
        m_sparseValuesVec.push_back(M[0][1]);
        m_sparseValuesVec.push_back(M[0][2]);

        m_sparseValuesVec.push_back(M[1][0]);
        m_sparseValuesVec.push_back(M[1][1]);
        m_sparseValuesVec.push_back(M[1][2]);

        m_sparseValuesVec.push_back(M[2][0]);
        m_sparseValuesVec.push_back(M[2][1]);
        m_sparseValuesVec.push_back(M[2][2]);
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
        std::vector<int> countvec;

        std::vector<int> tmp_colptr;//CRS format
        std::vector<int> tmp_rowind;

        std::vector<int> tmp_tran_colptr;//Transposed uncollapsed CRS format
        std::vector<int> tmp_tran_rowind;
        std::vector<int> tmp_tran_sortid;//index of the incoming values in the transposed format

        std::vector<int> valptr;// start and end point of each sum values in the CSR m_values vector
        std::vector<int> sortid;// corresponding indices in the m_sparseIndices vector


        unsigned tmp_nnz = m_sparseValuesVec.size();

        tmp_tran_colptr.resize(m_colSize+1);
        tmp_tran_rowind.resize(tmp_nnz);
        tmp_tran_sortid.resize(tmp_nnz);

        tmp_colptr.resize(m_rowSize+1);
        tmp_rowind.resize(tmp_nnz);
        sortid.resize(tmp_nnz);

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
                sortid[write_id] = tmp_tran_sortid[i]; // build the sortidvector

                countvec[row]++;
            }
        }

        m_rowind.clear();
        m_colptr.clear();
        m_colptr.push_back(0);

        //Compute the compressed row sparse format and collapse values
        for (unsigned j=0;j<m_rowSize;j++) {
            int newcol=-1;
            for (int i=tmp_colptr[j];i<tmp_colptr[j+1];i++) {
                // add new values only if needed
                if (tmp_rowind[i] != newcol) {
                    newcol = tmp_rowind[i];
                    m_rowind.push_back(newcol);
                    valptr.push_back(i);
                }
            }

            m_colptr.push_back(m_rowind.size());
        }

        valptr.push_back(tmp_colptr[m_rowSize]);

        //Build the mapping
        m_mapping.clear();
        //by default all the values will be sumed to the last value in the m_sparseValuesVec, we initialize the vector with the last index of the vector (size+1)
        m_mapping.resize(m_sparseValuesVec.size()+1,m_sparseValuesVec.size());
        for (unsigned i=0;i<m_rowind.size();i++) {
            for (unsigned j=valptr[i];j<valptr[i+1];j++) {
                m_mapping[sortid[j]] = i;
            }
        }

        //Build the values
        m_values.clear();
        m_values.resize(m_rowind.size());
        for (unsigned i=0;i<m_sparseValuesVec.size();i++) {
            m_values[m_mapping[i]] += m_sparseValuesVec[i];
        }
    }

private :
    VecReal m_sparseValuesVec; // incoming new values
    helper::vector<MatrixCoord> m_sparseIndices; // incoming vector of coordinated used to check the consistency of the structure

    unsigned & m_colSize;
    unsigned & m_rowSize;

    VecInt & m_colptr;//CRS format
    VecInt & m_rowind;
    VecReal & m_values;

    VecInt & m_mapping;

    helper::vector<int> m_clearCols,m_clearRows;
};

template<class VecReal,class VecInt>
class MappedIncomingMatrix : public defaulttype::BaseMatrix {
public:
    typedef typename RealType<VecReal>::Real Real;
    typedef defaulttype::BaseMatrix::Index Index;

    MappedIncomingMatrix(unsigned & r,unsigned & c,VecReal & values,VecInt & mapping)
    : m_rowSize(r)
    , m_colSize(c)
    , m_values(values)
    , m_mapping(mapping) {}

    Index rowSize(void) const override { return m_rowSize; }

    Index colSize(void) const override { return m_colSize; }

    SReal element(Index /*i*/, Index /*j*/) const override { return 0.0; }

    void resize(Index nbRow, Index nbCol) {
        m_rowSize = nbRow;
        m_colSize = nbCol;
    }

    void clear() {
        m_writeId=0;
        unsigned nnz = m_values.size();
        m_values.clear();
        m_values.resize(nnz);
    }

    void set(Index , Index , double ) override {}

    inline void add( Index , Index , double v ) override {
        m_values[m_mapping[m_writeId++]] += v;
    }

    inline void add(Index , Index , const defaulttype::Mat3x3d & M) override {
        m_values[m_mapping[m_writeId++]] += M[0][0];
        m_values[m_mapping[m_writeId++]] += M[0][1];
        m_values[m_mapping[m_writeId++]] += M[0][2];

        m_values[m_mapping[m_writeId++]] += M[1][0];
        m_values[m_mapping[m_writeId++]] += M[1][1];
        m_values[m_mapping[m_writeId++]] += M[1][2];

        m_values[m_mapping[m_writeId++]] += M[2][0];
        m_values[m_mapping[m_writeId++]] += M[2][1];
        m_values[m_mapping[m_writeId++]] += M[2][2];
    }

    void clearRow(Index ) override {}

    void clearCol(Index ) override {}

    void clearRowCol(Index ) override {}

    void clearRows(Index , Index ) override {}

    void clearCols(Index , Index ) override {}

    void clearRowsCols(Index , Index ) override {}

private :
    unsigned & m_colSize;
    unsigned & m_rowSize;
    unsigned m_writeId;

    VecReal & m_values;
    VecInt & m_mapping;
};

template<class VecReal, class VecInt>
class BaseLocalIncomingSparseMatrix : public sofa::core::objectmodel::BaseObject {
public:

    SOFA_CLASS(SOFA_TEMPLATE2(BaseLocalIncomingSparseMatrix,VecReal,VecInt),sofa::core::objectmodel::BaseObject);

    void resize(unsigned r,unsigned c) {
        m_rowSize = r;
        m_colSize = c;
    }

    void buildMatrix(std::vector<BaseStateAccessor::SPtr> & vacc) {
        if (m_mapping.empty()) {
            m_buildingMatrix.clear();

            component::linearsolver::DefaultMultiMatrixAccessor accessor;
            accessor.setGlobalMatrix(&m_buildingMatrix);
            accessor.setupMatrices();

            for (unsigned i=0;i<vacc.size();i++)
                accessor.addMechanicalState(vacc[i]->getState());

            doBuild(&accessor);

            m_buildingMatrix.rebuildPatternAccess();
        } else {
            m_mappedMatrix.clear();

            component::linearsolver::DefaultMultiMatrixAccessor accessor;
            accessor.setGlobalMatrix(&m_mappedMatrix);
            accessor.setupMatrices();

            for (unsigned i=0;i<vacc.size();i++)
                accessor.addMechanicalState(vacc[i]->getState());

            doBuild(&accessor);
        }
    }

    unsigned rowSize() { return m_rowSize; }

    unsigned colSize() { return m_colSize; }

    unsigned getNnz() { return m_rowind.size(); }

    VecInt & getColptr() { return m_colptr; }

    VecInt & getRowind() { return m_rowind; }

    VecReal & getValues() { return m_values; }

protected:
    BaseLocalIncomingSparseMatrix()
    : m_buildingMatrix(m_rowSize,m_colSize,m_colptr,m_rowind,m_values,m_mapping)
    , m_mappedMatrix(m_rowSize,m_colSize,m_values,m_mapping) {}

    virtual void doBuild(component::linearsolver::DefaultMultiMatrixAccessor * accessor) = 0;

    core::MechanicalParams m_params;
    unsigned m_rowSize,m_colSize;
    BuildingIncomingMatrix<VecReal,VecInt> m_buildingMatrix;
    MappedIncomingMatrix<VecReal,VecInt> m_mappedMatrix;

    VecInt m_colptr;//CRS format
    VecInt m_rowind;
    VecReal m_values;

    VecInt m_mapping;
};

template<class VecReal, class VecInt>
class LocalKIncomingSparseMatrix : public BaseLocalIncomingSparseMatrix<VecReal,VecInt> {
public:

    SOFA_CLASS(SOFA_TEMPLATE2(LocalKIncomingSparseMatrix,VecReal,VecInt),SOFA_TEMPLATE2(BaseLocalIncomingSparseMatrix,VecReal,VecInt));

    static void doCreateVisitor(std::vector<typename LocalKIncomingSparseMatrix::SPtr> & locals, core::objectmodel::BaseContext * ctx) {
        class BuildKLocalVisitor : public simulation::MechanicalVisitor {
        public:

            BuildKLocalVisitor(std::vector<typename LocalKIncomingSparseMatrix::SPtr> & locals)
            : simulation::MechanicalVisitor(core::MechanicalParams::defaultInstance())
            , m_locals(locals) {}

            const char* getClassName() const override { return "BuildLocalKVisitor"; }

            Result fwdForceField(simulation::Node* /*node*/, core::behavior::BaseForceField* ff) override {
                if (dynamic_cast<core::behavior::BaseMass*>(ff)) return RESULT_CONTINUE;

                auto slaves = ff->getSlaves();

                for (unsigned i=0;i<slaves.size();i++) {
                    if (typename LocalKIncomingSparseMatrix::SPtr local = sofa::core::objectmodel::SPtr_dynamic_cast<LocalKIncomingSparseMatrix>(slaves[i])) {
                        m_locals.push_back(local);
                        return RESULT_CONTINUE;
                    }
                }

                sofa::core::objectmodel::BaseObjectDescription arg;
                arg.setAttribute("type",std::string("LocalKIncomingSparseMatrix"));
                arg.setAttribute("template",typeid(VecReal).name());
                arg.setAttribute("name",ff->getName() + "_LIM");

                sofa::core::objectmodel::BaseObject::SPtr obj = sofa::core::ObjectFactory::getInstance()->createObject(ff->getContext(), &arg);
                LocalKIncomingSparseMatrix::SPtr local = sofa::core::objectmodel::SPtr_dynamic_cast<LocalKIncomingSparseMatrix>(obj);
                if (local== NULL) return RESULT_CONTINUE;
                local->setFF(ff);
                ff->addSlave(local);

                m_locals.push_back(local);

                return RESULT_CONTINUE;
            }

            bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* map) override {
                return !map->areMatricesMapped();
            }

        private:
            std::vector<typename LocalKIncomingSparseMatrix::SPtr> & m_locals;
        };

        locals.clear();
        BuildKLocalVisitor(locals).execute(ctx);
    }

    void setFF(core::behavior::BaseForceField* ff) {
        m_ff = ff;
    }

    virtual std::string getTemplateName() const override {
        return templateName(this);
    }

    static std::string templateName(const LocalKIncomingSparseMatrix* = NULL) {
        return typeid(VecReal).name();
    }

    void doBuild(component::linearsolver::DefaultMultiMatrixAccessor * accessor) override {
        m_ff->addKToMatrix(&this->m_params,accessor);
    }

protected:
    LocalKIncomingSparseMatrix()
    : BaseLocalIncomingSparseMatrix<VecReal,VecInt>()
    , m_ff(NULL){
        this->m_params.setMFactor(0.0);
        this->m_params.setBFactor(0.0);
        this->m_params.setKFactor(1.0);
    }

    core::behavior::BaseForceField * m_ff;
};

template<class VecReal, class VecInt>
class LocalMIncomingSparseMatrix : public BaseLocalIncomingSparseMatrix<VecReal,VecInt> {
public:

    SOFA_CLASS(SOFA_TEMPLATE2(LocalMIncomingSparseMatrix,VecReal,VecInt),SOFA_TEMPLATE2(BaseLocalIncomingSparseMatrix,VecReal,VecInt));

    static void doCreateVisitor(std::vector<typename LocalMIncomingSparseMatrix::SPtr> & locals, core::objectmodel::BaseContext * ctx) {
        class BuildMLocalVisitor : public simulation::MechanicalVisitor {
        public:

            BuildMLocalVisitor(std::vector<typename LocalMIncomingSparseMatrix::SPtr> & locals)
            : simulation::MechanicalVisitor(core::MechanicalParams::defaultInstance())
            , m_locals(locals) {}

            const char* getClassName() const override { return "BuildLocalMVisitor"; }

            Result fwdMass(simulation::Node* /*node*/, core::behavior::BaseMass* ff) override {
                auto slaves = ff->getSlaves();

                for (unsigned i=0;i<slaves.size();i++) {
                    if (typename LocalMIncomingSparseMatrix::SPtr local = sofa::core::objectmodel::SPtr_dynamic_cast<LocalMIncomingSparseMatrix>(slaves[i])) {
                        m_locals.push_back(local);
                        return RESULT_CONTINUE;
                    }
                }

//                typename LocalMIncomingSparseMatrix::SPtr local = core::objectmodel::New<LocalMIncomingSparseMatrix>(ff);

                sofa::core::objectmodel::BaseObjectDescription arg;
                arg.setAttribute("type",std::string("LocalMIncomingSparseMatrix"));
                arg.setAttribute("template",typeid(VecReal).name());
                arg.setAttribute("name",ff->getName() + "_LIM");

                sofa::core::objectmodel::BaseObject::SPtr obj = sofa::core::ObjectFactory::getInstance()->createObject(ff->getContext(), &arg);
                LocalMIncomingSparseMatrix::SPtr local = sofa::core::objectmodel::SPtr_dynamic_cast<LocalMIncomingSparseMatrix>(obj);
                if (local== NULL) return RESULT_CONTINUE;
                local->setFF(ff);
                ff->addSlave(local);

                m_locals.push_back(local);

                return RESULT_CONTINUE;
            }

            bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* map) override {
                return !map->areMatricesMapped();
            }

        private:
            std::vector<typename LocalMIncomingSparseMatrix::SPtr> & m_locals;
        };

        locals.clear();
        BuildMLocalVisitor(locals).execute(ctx);
    }

    void setFF(core::behavior::BaseMass * ff) {
        m_ff = ff;
    }

    virtual std::string getTemplateName() const override {
        return templateName(this);
    }

    static std::string templateName(const LocalMIncomingSparseMatrix* = NULL) {
        return typeid(VecReal).name();
    }

    void doBuild(component::linearsolver::DefaultMultiMatrixAccessor * accessor) override {
        m_ff->addMToMatrix(&this->m_params,accessor);
    }

protected:
    LocalMIncomingSparseMatrix()
    : BaseLocalIncomingSparseMatrix<VecReal,VecInt>()
    , m_ff(NULL){
        this->m_params.setMFactor(1.0);
        this->m_params.setBFactor(0.0);
        this->m_params.setKFactor(0.0);
    }

    core::behavior::BaseMass * m_ff;
};







}
