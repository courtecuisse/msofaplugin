#include <sofa/msofaplugin/matrix/StateAccessor.h>
#include <sofa/type/Vec.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::msofaplugin::matrix {

int StateVectorOperatorClass = core::RegisterObject("Direct Linear Solver using a Sparse LDL^T factorization.")
.add< StateAccessor< defaulttype::Vec1Types > >()
.add< StateAccessor< defaulttype::Vec2Types > >()
.add< StateAccessor< defaulttype::Vec3Types > >()
.add< StateAccessor< defaulttype::Vec6Types > >()
.add< StateAccessor< defaulttype::Rigid3Types > >()
.add< StateAccessor< defaulttype::Rigid2Types > >();

}
