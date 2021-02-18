#include <sofa/core/ObjectFactory.h>
#include <sofa/msofaplugin/matrix/IncomingSparseMatrix.h>

namespace sofa::msofaplugin::matrix {

int IncomingBaseMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< IncomingSparseMatrix< double > >()
.add< IncomingSparseMatrix< float > >()
;

}
