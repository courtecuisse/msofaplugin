#pragma once

#include <sofa/msofaplugin/solver/MBaseLinearSolver.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/helper/map.h>

namespace sofa::msofaplugin::solver {

class MPCGLinearSolver : public MBaseLinearSolver {
public:

    SOFA_CLASS(MPCGLinearSolver, MBaseLinearSolver);

    Data<std::map < std::string, sofa::helper::vector<double> > > d_graph; ///< Graph of residuals at each iteration
    Data<int> d_iteration;
    Data<double> d_tolerance;
    Data<double> d_threshold;
    Data<bool> d_usePrecond; ///< Graph of residuals at each iteration
    core::objectmodel::SingleLink<MPCGLinearSolver,BaseSystemMatrix,BaseLink::FLAG_STRONGLINK|BaseLink::FLAG_STOREPATH> l_matrix;
    core::objectmodel::SingleLink<MPCGLinearSolver,MBaseLinearSolver,BaseLink::FLAG_STRONGLINK|BaseLink::FLAG_STOREPATH> l_preconditioner;

    MPCGLinearSolver()
    : d_graph(initData(&d_graph, "graph", "Error graph evolution"))
    , d_iteration(initData(&d_iteration, 100,"iterations","Number of iterations"))
    , d_tolerance(initData(&d_tolerance, 0.00001,"tolerance","Tolerance of the CG"))
    , d_threshold(initData(&d_threshold, 0.00001,"threshold","Threshold of the CG"))
    , d_usePrecond(initData(&d_usePrecond, true, "usePrecond", "Activate or desactivate the preconditioner"))
    , l_matrix(initLink("matrix", "Link to the matrix"))
    , l_preconditioner(initLink("preconditioner", "Link to the perconditioner (optional)"))
    {
        l_matrix.setPath("@.");
    }

    virtual void buildSystemMatrix() override {
        this->l_matrix->buildSystemMatrix(); // unbuilt:0 ms    eigen : 14 ms
        if (d_usePrecond.getValue() && l_preconditioner != NULL) l_preconditioner->setSystemMBKMatrix(this->getMParams());
    }

    //compute x=inv(A) b
    virtual void solve(BaseSystemMatrix::MechanicalVectorId & x, BaseSystemMatrix::MechanicalVectorId & b) override {
        sofa::helper::AdvancedTimer::stepBegin("SolveSystem");

        BaseSystemMatrix::MechanicalVectorId q = getMatrix()->createMVecId();
        BaseSystemMatrix::MechanicalVectorId d = getMatrix()->createMVecId();
        BaseSystemMatrix::MechanicalVectorId r = getMatrix()->createMVecId();

        this->l_matrix->mult(&this->m_params, q, x); // 0.28    0.6
//        auto A = this->l_matrix->getCompressedMatrix(&this->m_params); // 1.5 ms
//        A->mult(q,x); //0.24 ms

        this->l_matrix->eq(r,b,q,-1); // r = b - A*x

        if (! d_usePrecond.getValue() || l_preconditioner == NULL) this->l_matrix->eq(d,r); // d=r
        else {
            helper::AdvancedTimer::stepBegin("PCG - apply precond");
            l_preconditioner->solve(d,r);
            helper::AdvancedTimer::stepEnd("PCG - apply precond");
        }


        double d_new = this->l_matrix->dot(r,d); // d_new=r*d

        double d_0 = d_new;

        double tolerance;
        if(d_tolerance.getValue() < 0.0)
            tolerance = -1.0;
        else
            tolerance = d_tolerance.getValue()*d_tolerance.getValue()*d_0;

        sofa::helper::vector<double> & graph = (*d_graph.beginEdit())["tolerance"];
        graph.clear();

        int cg_iter = 0;
        while (cg_iter < d_iteration.getValue() && d_new > tolerance)
        {
//            helper::AdvancedTimer::stepBegin("CG ite runtime");
            graph.push_back(d_new);

//            sofa::helper::AdvancedTimer::stepBegin("SpMV");
            // q = Ad;
            this->l_matrix->mult(&this->m_params,q,d); // 0.28 ms
//            A->mult(q,d); // 0.24 ms
//            sofa::helper::AdvancedTimer::stepEnd("SpMV");

//            sofa::helper::AdvancedTimer::stepBegin("dot");
            double den = this->l_matrix->dot(d,q);
//            sofa::helper::AdvancedTimer::stepEnd("dot");


            if( fabs(den)<d_threshold.getValue()) break;
//            sofa::helper::AdvancedTimer::stepBegin("peq");

            double alpha = d_new / den;

            this->l_matrix->peq(x,d,alpha);

            this->l_matrix->peq(r,q,-alpha);

            double d_old = d_new;
//            sofa::helper::AdvancedTimer::stepEnd("peq");

            if (! d_usePrecond.getValue() || l_preconditioner == NULL) {
//                sofa::helper::AdvancedTimer::stepBegin("dot");
                d_new = this->l_matrix->dot(r,r); // r*r
//                sofa::helper::AdvancedTimer::stepEnd("dot");

                double beta = d_new / d_old;

                this->l_matrix->eq(d,r,d,beta);

            } else {
                helper::AdvancedTimer::stepBegin("PCG - apply precond");
                l_preconditioner->solve(q,r);
                helper::AdvancedTimer::stepEnd("PCG - apply precond");

                d_new = this->l_matrix->dot(r,q); // d_new=r * q

                double beta = d_new / d_old;

                this->l_matrix->eq(d,q,d,beta);
            }

            ++cg_iter;
//            helper::AdvancedTimer::stepEnd("CG ite runtime");
        }
        d_graph.endEdit();

        sofa::helper::AdvancedTimer::stepEnd("SolveSystem");

        sofa::helper::AdvancedTimer::valSet("PCG rho", d_new);
        sofa::helper::AdvancedTimer::valSet("PCG iterations", cg_iter);
    }    

    BaseSystemMatrix * getMatrix() override {
        return this->l_matrix.get();
    }

};

}
