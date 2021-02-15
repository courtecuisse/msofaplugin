#include <sofa/core/ObjectFactory.h>
#include <sofa/msofaplugin/matrix/EigenSparseMatrix.h>

namespace sofa::msofaplugin::matrix {

int EigenSparseMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< EigenSparseMatrix< double > >()
.add< EigenSparseMatrix< float > >()
;

}
