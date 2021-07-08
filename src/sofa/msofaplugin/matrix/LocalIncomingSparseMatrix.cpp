#include <sofa/core/ObjectFactory.h>
#include <sofa/msofaplugin/matrix/LocalIncomingSparseMatrix.h>

namespace sofa::msofaplugin::matrix {

int LocalMIncomingBaseMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< LocalIncomingSparseMatrix< sofa::type::vector<double>, sofa::type::vector<int> > >()
.add< LocalIncomingSparseMatrix< sofa::type::vector<float>, sofa::type::vector<int> > >()
;

}
