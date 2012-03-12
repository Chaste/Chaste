/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
            #define COVERAGE_IGNORE
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


