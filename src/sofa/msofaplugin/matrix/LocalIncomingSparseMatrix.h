#pragma once

#include <sofa/msofaplugin/matrix/BaseSystemMatrix.h>
#include <sofa/msofaplugin/matrix/CompressableMatrix.h>
#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/BaseMass.h>
#include <sofa/linearalgebra/BaseMatrix.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/msofaplugin/matrix/MechanicalVector.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>
#include <sofa/simulation/MechanicalOperations.h>
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
class BuildingIncomingMatrix : public sofa::linearalgebra::BaseMatrix {
public:
    typedef typename RealType<VecReal>::Real Real;
    typedef sofa::linearalgebra::BaseMatrix::Index Index;

    struct MatrixCoord {
        MatrixCoord() {}
        MatrixCoord(unsigned r,unsigned c) : row(r), col(c) {}
        unsigned row;
        unsigned col;
    };

    BuildingIncomingMatrix(unsigned & r,unsigned & c,VecInt & colptr,VecInt & rowind,VecReal & values,std::vector<int> & mapping)
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
        m_colptr.clear();
        m_rowind.clear();
        m_values.clear();
        m_mapping.clear();
        m_rowSize = 0;
        m_colSize = 0;
    }

    void set(Index , Index , double ) override {}

    inline void add( Index i, Index j, double v ) override {
        m_sparseIndices.push_back(MatrixCoord(i,j));
        m_sparseValuesVec.push_back(v);
    }

    inline void add(Index i, Index j, const type::Mat3x3d & M) override {
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

    void clearRow(Index ) override {}

    void clearCol(Index ) override {}

    void clearRowCol(Index ) override {}

    void clearRows(Index , Index ) override {}

    void clearCols(Index , Index ) override {}

    void clearRowsCols(Index , Index ) override {}

    void inline rebuildPatternAccess(sofa::type::vector<int> & clearCols, sofa::type::vector<int> & clearRows) {
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
        for (unsigned i=0;i<clearCols.size();i++) countvec[clearCols[i]] = 0;

        tmp_tran_colptr[0] = 0;
        for (unsigned j=0;j<m_colSize;j++) tmp_tran_colptr[j+1] = tmp_tran_colptr[j] + countvec[j];

        countvec.clear();
        countvec.resize(m_colSize); // reset to zero

        for (unsigned i=0;i<tmp_nnz;i++) {
            unsigned col = m_sparseIndices[i].col;

            int write_id = tmp_tran_colptr[col] + countvec[col];
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
            for (int i=tmp_tran_colptr[j];i<tmp_tran_colptr[j+1];i++) {
                int row = tmp_tran_rowind[i];
                countvec[row]++;
            }
        }

        //clear the columns set the count values to 0
        for (unsigned i=0;i<clearRows.size();i++) countvec[clearRows[i]] = 0;

        tmp_colptr[0] = 0;
        for (unsigned j=0;j<m_rowSize;j++) tmp_colptr[j+1] = tmp_colptr[j] + countvec[j];

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
        m_mapping.resize(m_sparseValuesVec.size(),m_rowind.size());
        for (unsigned i=0;i<m_rowind.size();i++) {
            for (unsigned j=valptr[i];j<valptr[i+1];j++) {
                m_mapping[sortid[j]] = i;
            }
        }

        //Build the values
        m_values.clear();
        m_values.resize(m_rowind.size()+1); // +1 to store contribution of projected points
        for (unsigned i=0;i<m_sparseValuesVec.size();i++) {
            m_values[m_mapping[i]] += m_sparseValuesVec[i];
        }
    }

private :
    std::vector<Real> m_sparseValuesVec; // incoming new values
    sofa::type::vector<MatrixCoord> m_sparseIndices; // incoming vector of coordinated used to check the consistency of the structure

    unsigned & m_rowSize;
    unsigned & m_colSize;

    VecInt & m_colptr;//CRS format
    VecInt & m_rowind;
    VecReal & m_values;

    std::vector<int> & m_mapping;
};

template<class VecReal,class VecInt>
class MappedIncomingMatrix : public sofa::linearalgebra::BaseMatrix {
public:
    typedef typename RealType<VecReal>::Real Real;
    typedef sofa::linearalgebra::BaseMatrix::Index Index;

    MappedIncomingMatrix(unsigned & r,unsigned & c,VecReal & values,const std::vector<int> & mapping)
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
        m_valuesPtr = m_values.data();
        memset(m_valuesPtr,0,m_values.size()*sizeof(Real));
    }

    void set(Index , Index , double ) override {}

    inline void add( Index , Index , double v ) override {
        m_valuesPtr[m_mapping[m_writeId++]] += v;
    }

    inline void add(Index , Index , const type::Mat3x3d & M) override {
        m_valuesPtr[m_mapping[m_writeId++]] += M[0][0];
        m_valuesPtr[m_mapping[m_writeId++]] += M[0][1];
        m_valuesPtr[m_mapping[m_writeId++]] += M[0][2];

        m_valuesPtr[m_mapping[m_writeId++]] += M[1][0];
        m_valuesPtr[m_mapping[m_writeId++]] += M[1][1];
        m_valuesPtr[m_mapping[m_writeId++]] += M[1][2];

        m_valuesPtr[m_mapping[m_writeId++]] += M[2][0];
        m_valuesPtr[m_mapping[m_writeId++]] += M[2][1];
        m_valuesPtr[m_mapping[m_writeId++]] += M[2][2];
    }

    void clearRow(Index ) override {}

    void clearCol(Index ) override {}

    void clearRowCol(Index ) override {}

    void clearRows(Index , Index ) override {}

    void clearCols(Index , Index ) override {}

    void clearRowsCols(Index , Index ) override {}

private :
    unsigned m_writeId;
    unsigned & m_rowSize;
    unsigned & m_colSize;

    VecReal & m_values;
    const std::vector<int> & m_mapping;
    Real * m_valuesPtr; // this pointer may be accessed on GPU
};

template<class Real>
static std::string realName(const Real * );

template<>
std::string realName<float>(const float * ) { return "float";  }

template<>
std::string realName<double>(const double * ) { return "double";  }


template<class VecReal, class VecInt>
class LocalIncomingSparseMatrix : public sofa::core::objectmodel::BaseObject {
public:

    class BaseInternalBuilder {
    public:
        typedef std::unique_ptr<BaseInternalBuilder> UPtr;

        virtual void build(component::linearsolver::DefaultMultiMatrixAccessor * accessor) = 0;

        virtual double getFactor(const core::MechanicalParams * params) = 0;
    };

    class KBuilder : public BaseInternalBuilder {
    public:
        KBuilder(core::behavior::BaseForceField* f) : m_ff(f) {}

        void build(component::linearsolver::DefaultMultiMatrixAccessor *accessor) override {
            core::MechanicalParams mparams;
            mparams.setMFactor(0.0);mparams.setBFactor(0.0);mparams.setKFactor(1.0);
            m_ff->addKToMatrix(&mparams,accessor);
        };

        double getFactor(const core::MechanicalParams *params) override { return params->kFactor(); };
    protected:
        core::behavior::BaseForceField* m_ff;
    };

    class MBuilder : public BaseInternalBuilder {
    public:
        MBuilder(core::behavior::BaseMass * m) : m_ff(m) {}

        void build(component::linearsolver::DefaultMultiMatrixAccessor *accessor) override {
            core::MechanicalParams mparams;
            mparams.setMFactor(1.0);mparams.setBFactor(0.0);mparams.setKFactor(0.0);
            m_ff->addMToMatrix(&mparams,accessor);
        };

        double getFactor(const core::MechanicalParams *params) override { return params->mFactor(); };
    protected:
        core::behavior::BaseMass * m_ff;
    };

    SOFA_CLASS(SOFA_TEMPLATE2(LocalIncomingSparseMatrix,VecReal,VecInt),sofa::core::objectmodel::BaseObject);

    void clear() {
        m_buildingMatrix.clear();
    }

    void resize(unsigned r,unsigned c) {
        m_rowSize = r;
        m_colSize = c;
    }

    inline void buildMatrix(std::vector<BaseStateAccessor::SPtr> & vacc,sofa::type::vector<int> & clearCols,sofa::type::vector<int> & clearRows) {
        m_accessor.clear();
        m_accessor.setGlobalMatrix(&m_buildingMatrix);
        m_accessor.setupMatrices();

        for (unsigned i=0;i<vacc.size();i++)
            m_accessor.addMechanicalState(vacc[i]->getState());

        m_builder->build(&m_accessor);

        m_buildingMatrix.rebuildPatternAccess(clearCols,clearRows);
    }

    inline void fastReBuild() {
        m_mappedMatrix.clear();

        //prepare the next step using the compressed matrix
        m_accessor.setGlobalMatrix(&m_mappedMatrix);

        m_builder->build(&m_accessor);
    }

    unsigned rowSize() { return m_rowSize; }

    unsigned colSize() { return m_colSize; }

    unsigned getNnz() { return m_rowind.size(); }

    VecInt & getColptr() { return m_colptr; }

    VecInt & getRowind() { return m_rowind; }

    VecReal & getValues() { return m_values; }

    static void doCreateVisitor(std::vector<typename LocalIncomingSparseMatrix::SPtr> & locals, core::objectmodel::BaseContext * ctx) {
        class BuildLocalMatrixVisitor : public simulation::MechanicalVisitor {
        public:

            BuildLocalMatrixVisitor(std::vector<typename LocalIncomingSparseMatrix::SPtr> & locals)
            : simulation::MechanicalVisitor(core::MechanicalParams::defaultInstance())
            , m_locals(locals) {}

            const char* getClassName() const override { return "BuildLocalKVisitor"; }

            LocalIncomingSparseMatrix::SPtr getLocal(core::objectmodel::BaseObject * obj) {
                auto slaves = obj->getSlaves();

                for (unsigned i=0;i<slaves.size();i++) {
                    if (typename LocalIncomingSparseMatrix::SPtr local = sofa::core::objectmodel::SPtr_dynamic_cast<LocalIncomingSparseMatrix>(slaves[i])) {
                        return local;
                    }
                }

                sofa::core::objectmodel::BaseObjectDescription arg;
                arg.setAttribute("type",std::string("LocalIncomingSparseMatrix"));
                arg.setAttribute("template",realName(VecReal().data()));
                arg.setAttribute("name",obj->getName() + "_LIM");

                sofa::core::objectmodel::BaseObject::SPtr newobj = sofa::core::ObjectFactory::getInstance()->createObject(obj->getContext(), &arg);
                obj->addSlave(newobj);
                return sofa::core::objectmodel::SPtr_dynamic_cast<LocalIncomingSparseMatrix>(newobj);
            }

            Result fwdMass(simulation::Node* /*node*/, core::behavior::BaseMass* ff) override {
                LocalIncomingSparseMatrix::SPtr local = getLocal(ff);
                if (local== NULL) return RESULT_CONTINUE;
                local->setBuilder(new MBuilder(ff));
                m_locals.push_back(local);
                return RESULT_CONTINUE;
            }

            Result fwdForceField(simulation::Node* /*node*/, core::behavior::BaseForceField* ff) override {
                if (dynamic_cast<core::behavior::BaseMass*>(ff)) return RESULT_CONTINUE;

                LocalIncomingSparseMatrix::SPtr local = getLocal(ff);
                if (local== NULL) return RESULT_CONTINUE;
                local->setBuilder(new KBuilder(ff));
                m_locals.push_back(local);
                return RESULT_CONTINUE;
            }

            bool stopAtMechanicalMapping(simulation::Node* /*node*/, core::BaseMapping* map) override {
                return !map->areMatricesMapped();
            }

        private:
            std::vector<typename LocalIncomingSparseMatrix::SPtr> & m_locals;
        };

        locals.clear();
        BuildLocalMatrixVisitor(locals).execute(ctx);
    }

    double getFactor(const core::MechanicalParams * params) { return m_builder->getFactor(params); }

    void setBuilder(BaseInternalBuilder * b) { m_builder = typename BaseInternalBuilder::UPtr(b); }

    virtual std::string getTemplateName() const override {
        return templateName(this);
    }

    static std::string templateName(const LocalIncomingSparseMatrix<VecReal,VecInt>* = NULL) {
        return realName(VecReal().data());
    }

protected:
    LocalIncomingSparseMatrix()
    : m_buildingMatrix(m_rowSize,m_colSize,m_colptr,m_rowind,m_values,m_mapping)
    , m_mappedMatrix(m_rowSize,m_colSize,m_values,m_mapping)
    , m_rowSize(0)
    , m_colSize(0) {}

    BuildingIncomingMatrix<VecReal,VecInt> m_buildingMatrix;
    MappedIncomingMatrix<VecReal,VecInt> m_mappedMatrix;
    unsigned m_rowSize,m_colSize;

    VecInt m_colptr;//CRS format
    VecInt m_rowind;
    VecReal m_values;

    std::vector<int> m_mapping;

    component::linearsolver::DefaultMultiMatrixAccessor m_accessor;
    typename BaseInternalBuilder::UPtr m_builder;
};



}
