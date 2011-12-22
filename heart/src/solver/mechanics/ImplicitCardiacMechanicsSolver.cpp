/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "ImplicitCardiacMechanicsSolver.hpp"
#include "Kerchoffs2003ContractionModel.hpp"
#include "NhsModelWithBackwardSolver.hpp"
#include "NonPhysiologicalContractionModel.hpp"

template<class ELASTICITY_SOLVER,unsigned DIM>
ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::ImplicitCardiacMechanicsSolver(
                                  ContractionModelName contractionModel,
                                  QuadraticMesh<DIM>& rQuadMesh,
                                  SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                  std::string outputDirectory)
    : AbstractCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>(rQuadMesh,
                                                            rProblemDefinition,
                                                            outputDirectory)
{
    InitialiseContractionModels(contractionModel);

    // initialise stores
    mStretchesLastTimeStep.resize(this->mTotalQuadPoints, 1.0);
}


template<class ELASTICITY_SOLVER,unsigned DIM>
void ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::InitialiseContractionModels(ContractionModelName contractionModel)
{
    switch(contractionModel)
    {
        case NONPHYSIOL1:
        case NONPHYSIOL2:
        case NONPHYSIOL3:
        {
            unsigned option = (contractionModel==NONPHYSIOL1 ? 1 : (contractionModel==NONPHYSIOL2? 2 : 3));
            for(unsigned i=0; i<this->mTotalQuadPoints; i++)
            {
                this->mContractionModelSystems.push_back(new NonPhysiologicalContractionModel(option));
            }
            break;
        }
        case NHS:
        {
            for(unsigned i=0; i<this->mTotalQuadPoints; i++)
            {
                this->mContractionModelSystems.push_back(new NhsModelWithBackwardSolver);
            }
            break;
        }
        case KERCHOFFS2003: //stretch dependent
        {
            for(unsigned i=0; i<this->mTotalQuadPoints; i++)
            {
                this->mContractionModelSystems.push_back(new Kerchoffs2003ContractionModel());
            }
            break;
        }
        default:
        {
            #define COVERAGE_IGNORE // currently all available contraction models are acceptable for implicit
            EXCEPTION("Unknown or disallowed contraction model");
            #undef COVERAGE_IGNORE
        }
    }
}

template<class ELASTICITY_SOLVER,unsigned DIM>
ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::~ImplicitCardiacMechanicsSolver()
{
    for(unsigned i=0; i<this->mContractionModelSystems.size(); i++)
    {
        delete this->mContractionModelSystems[i];
    }
}


template<class ELASTICITY_SOLVER,unsigned DIM>
std::vector<double>& ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::rGetFibreStretches()
{
    return this->mStretches;
}


template<class ELASTICITY_SOLVER,unsigned DIM>
void ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::Solve(double time, double nextTime, double odeTimestep)
{
    // set the times, which are used in AssembleOnElement
    assert(time < nextTime);
    this->mCurrentTime = time;
    this->mNextTime = nextTime;
    this->mOdeTimestep = odeTimestep;

    // solve
    ELASTICITY_SOLVER::Solve();

    // assemble residual again (to solve the cell models implicitly again
    // using the correct value of the deformation x (in case this wasn't the
    // last thing that was done
    if(this->GetNumNewtonIterations() > 0)
    {
        this->AssembleSystem(true,false);
    }

    // now update state variables, and set lambda at last timestep. Note
    // lambda was set in AssembleOnElement
    for(unsigned i=0; i<this->mContractionModelSystems.size(); i++)
    {
         this->mContractionModelSystems[i]->UpdateStateVariables();
         mStretchesLastTimeStep[i] = this->mStretches[i];
    }
}



template<class ELASTICITY_SOLVER,unsigned DIM>
void ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::GetActiveTensionAndTensionDerivs(double currentFibreStretch,
                                                                                             unsigned currentQuadPointGlobalIndex,
                                                                                             bool assembleJacobian,
                                                                                             double& rActiveTension,
                                                                                             double& rDerivActiveTensionWrtLambda,
                                                                                             double& rDerivActiveTensionWrtDLambdaDt)
{
    // save this fibre stretch
    this->mStretches[currentQuadPointGlobalIndex] = currentFibreStretch;

    // compute dlam/dt
    double dlam_dt = (currentFibreStretch-mStretchesLastTimeStep[currentQuadPointGlobalIndex])/(this->mNextTime-this->mCurrentTime);

    AbstractContractionModel* p_contraction_model = this->mContractionModelSystems[currentQuadPointGlobalIndex];

    // Set this stretch and stretch rate
    p_contraction_model->SetStretchAndStretchRate(currentFibreStretch, dlam_dt);

    // Call RunDoNotUpdate() on the contraction model to solve it using this stretch, and get the active tension
    try
    {
        p_contraction_model->RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
        rActiveTension = p_contraction_model->GetNextActiveTension();
    }
    catch (Exception& e)
    {
        #define COVERAGE_IGNORE
        // if this failed during assembling the Jacobian this is a fatal error.
        if(assembleJacobian)
        {
            // probably shouldn't be able to get here
            EXCEPTION("Failure in solving contraction models using current stretches for assembling Jacobian");
        }
        // if this failed during assembling the residual, the stretches are too large, so we just
        // set the active tension to infinity so that the residual will be infinite
        rActiveTension = DBL_MAX;
        std::cout << "WARNING: could not solve contraction model with this stretch and stretch rate. "
                  << "Setting active tension to infinity (DBL_MAX) so that the residual(-norm) is also infinite\n" << std::flush;
        assert(0); // just to see if we ever get here, can be removed..
        return;
        #undef COVERAGE_IGNORE
    }

    // if assembling the Jacobian, numerically evaluate dTa/dlam & dTa/d(lamdot)
    if(assembleJacobian)
    {
        // get active tension for (lam+h,dlamdt)
        double h1 = std::max(1e-6, currentFibreStretch/100);
        p_contraction_model->SetStretchAndStretchRate(currentFibreStretch+h1, dlam_dt);
        p_contraction_model->RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
        double active_tension_at_lam_plus_h = p_contraction_model->GetNextActiveTension();

        // get active tension for (lam,dlamdt+h)
        double h2 = std::max(1e-6, dlam_dt/100);
        p_contraction_model->SetStretchAndStretchRate(currentFibreStretch, dlam_dt+h2);
        p_contraction_model->RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
        double active_tension_at_dlamdt_plus_h = p_contraction_model->GetNextActiveTension();

        rDerivActiveTensionWrtLambda = (active_tension_at_lam_plus_h - rActiveTension)/h1;
        rDerivActiveTensionWrtDLambdaDt = (active_tension_at_dlamdt_plus_h - rActiveTension)/h2;
    }

    // Re-set the stretch and stretch rate and recompute the active tension so that
    // if this guess turns out to the solution, we can just update the state variables
    //  -- not needed as AssembleSystem(true,false) [ie assemble residual] is
    //     called in ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::Solve() above
    //     after the solve and before the update.
    //  -- The SetActiveTensionInitialGuess() would make this very fast
    //     (compared to AssembleSystem(true,false) above), but the NHS class uses the last
    //     active tension as the initial guess anyway..
    //p_contraction_model->SetStretchAndStretchRate(currentFibreStretch, dlam_dt);
    //p_contraction_model->SetActiveTensionInitialGuess(rActiveTension);   // this line would now need a cast
    //p_contraction_model->RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
    //assert( fabs(p_contraction_model->GetNextActiveTension()-rActiveTension)<1e-8);
}



template class ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2>;
template class ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<3>,3>;
template class ImplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<2>,2>;
template class ImplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<3>,3>;


