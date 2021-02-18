#include <sofa/core/ObjectFactory.h>
#include <sofa/msofaplugin/matrix/LocalIncomingSparseMatrix.h>

namespace sofa::msofaplugin::matrix {

int LocalMIncomingBaseMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< LocalMIncomingSparseMatrix< helper::vector<double>, helper::vector<int> > >()
.add< LocalMIncomingSparseMatrix< helper::vector<double>, helper::vector<int> > >()
;

int LocalKIncomingBaseMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< LocalKIncomingSparseMatrix< helper::vector<double>, helper::vector<int> > >()
.add< LocalKIncomingSparseMatrix< helper::vector<double>, helper::vector<int> > >()
;

}
