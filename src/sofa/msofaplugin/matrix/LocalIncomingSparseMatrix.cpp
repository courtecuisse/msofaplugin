#include <sofa/core/ObjectFactory.h>
#include <sofa/msofaplugin/matrix/LocalIncomingSparseMatrix.h>

namespace sofa::msofaplugin::matrix {

int LocalMIncomingBaseMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< LocalIncomingSparseMatrix< helper::vector<double>, helper::vector<int> > >()
.add< LocalIncomingSparseMatrix< helper::vector<float>, helper::vector<int> > >()
;

}
