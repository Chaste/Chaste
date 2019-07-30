/*

Copyright (c) 2005-2019, University of Oxford.
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

template<class ELASTICITY_SOLVER,unsigned DIM>
ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::ImplicitCardiacMechanicsSolver(
                                  QuadraticMesh<DIM>& rQuadMesh,
                                  ElectroMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                  std::string outputDirectory)
    : AbstractCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>(rQuadMesh,
                                                            rProblemDefinition,
                                                            outputDirectory)
{

}

template<class ELASTICITY_SOLVER,unsigned DIM>
ImplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::~ImplicitCardiacMechanicsSolver()
{
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
    if (this->GetNumNewtonIterations() > 0)
    {
        this->AssembleSystem(true,false);
    }

    // now update state variables, and set lambda at last timestep. Note
    // stretches were set in AssembleOnElement
    for (std::map<unsigned,DataAtQuadraturePoint>::iterator iter = this->mQuadPointToDataAtQuadPointMap.begin();
         iter != this->mQuadPointToDataAtQuadPointMap.end();
         iter++)
    {
        AbstractContractionModel* p_contraction_model = iter->second.ContractionModel;
        iter->second.StretchLastTimeStep = iter->second.Stretch;
        p_contraction_model->UpdateStateVariables();
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
    // The iterator should be pointing to the right place (note: it is incremented at the end of this method)
    // This iterator is used so that we don't have to search the map
    assert(this->mMapIterator->first==currentQuadPointGlobalIndex);

    DataAtQuadraturePoint& r_data_at_quad_point = this->mMapIterator->second;

    // save this fibre stretch
    r_data_at_quad_point.Stretch = currentFibreStretch;

    // compute dlam/dt
    double dlam_dt = (currentFibreStretch-r_data_at_quad_point.StretchLastTimeStep)/(this->mNextTime-this->mCurrentTime);

    // Set this stretch and stretch rate on contraction model
    AbstractContractionModel* p_contraction_model = r_data_at_quad_point.ContractionModel;
    p_contraction_model->SetStretchAndStretchRate(currentFibreStretch, dlam_dt);

    // Call RunDoNotUpdate() on the contraction model to solve it using this stretch, and get the active tension
    try
    {
        p_contraction_model->RunDoNotUpdate(this->mCurrentTime,this->mNextTime,this->mOdeTimestep);
        rActiveTension = p_contraction_model->GetNextActiveTension();
    }
    // LCOV_EXCL_START
    catch (Exception&)
    {
        // if this failed during assembling the Jacobian this is a fatal error.
        if (assembleJacobian)
        {
            // probably shouldn't be able to get here
            EXCEPTION("Failure in solving contraction models using current stretches for assembling Jacobian");
        }
        // if this failed during assembling the residual, the stretches are too large, so we just
        // set the active tension to infinity so that the residual will be infinite
        rActiveTension = DBL_MAX;
        std::cout << "WARNING: could not solve contraction model with this stretch and stretch rate. "
                  << "Setting active tension to infinity (DBL_MAX) so that the residual(-norm) is also infinite\n" << std::flush;
        NEVER_REACHED;
        return;
    }
    // LCOV_EXCL_STOP

    // if assembling the Jacobian, numerically evaluate dTa/dlam & dTa/d(lamdot)
    if (assembleJacobian)
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

    // Increment the iterator
    this->mMapIterator++;
    if (this->mMapIterator==this->mQuadPointToDataAtQuadPointMap.end())
    {
        this->mMapIterator = this->mQuadPointToDataAtQuadPointMap.begin();
    }
}

// Explicit instantiation
template class ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2>;
template class ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<3>,3>;
template class ImplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<2>,2>;
template class ImplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<3>,3>;
