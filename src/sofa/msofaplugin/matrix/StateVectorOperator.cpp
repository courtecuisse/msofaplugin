#include <sofa/msofaplugin/matrix/StateVectorOperator.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::msofaplugin::matrix {

int StateVectorOperatorClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< StateVectorOperator< defaulttype::Vec1Types > >()
.add< StateVectorOperator< defaulttype::Vec2Types > >()
.add< StateVectorOperator< defaulttype::Vec3Types > >()
.add< StateVectorOperator< defaulttype::Vec6Types > >()
.add< StateVectorOperator< defaulttype::Rigid3Types > >()
.add< StateVectorOperator< defaulttype::Rigid2Types > >();

}
