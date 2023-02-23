#include <sofa/msofaplugin/solver/MShewchukPCGLinearSolver.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/linearsolver/iterative/GraphScatteredTypes.h>

namespace sofa::msofaplugin::solver {

int MPCGSolverClass = core::RegisterObject("MShewchukPCGLinearSolver")
.add<MShewchukPCGLinearSolver>();

}
