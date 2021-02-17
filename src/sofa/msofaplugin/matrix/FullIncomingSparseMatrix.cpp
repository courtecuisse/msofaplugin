#include <sofa/core/ObjectFactory.h>
#include <sofa/msofaplugin/matrix/FullIncomingSparseMatrix.h>

namespace sofa::msofaplugin::matrix {

int IncomingBaseMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< FullIncomingSparseMatrix< double > >()
.add< FullIncomingSparseMatrix< float > >()
;

}
