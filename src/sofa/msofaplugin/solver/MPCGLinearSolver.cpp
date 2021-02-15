#include <sofa/msofaplugin/solver/MPCGLinearSolver.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseLinearSolver/GraphScatteredTypes.h>

namespace sofa::msofaplugin::solver {

int MPCGSolverClass = core::RegisterObject("MPCGLinearSolver")
.add<MPCGLinearSolver>();

}
