#include <sofa/core/ObjectFactory.h>
#include <sofa/msofaplugin/matrix/MUnbuiltMatrix.h>

namespace sofa::msofaplugin::matrix {

int UnbuiltMatrixClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< MUnbuiltMatrix >()
;

}
