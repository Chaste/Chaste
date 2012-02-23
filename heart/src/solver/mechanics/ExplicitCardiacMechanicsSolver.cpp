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

#include "ExplicitCardiacMechanicsSolver.hpp"
#include "Nash2004ContractionModel.hpp"
#include "Kerchoffs2003ContractionModel.hpp"
#include "NonPhysiologicalContractionModel.hpp"


template<class ELASTICITY_SOLVER,unsigned DIM>
ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::ExplicitCardiacMechanicsSolver(ContractionModelName contractionModel,
                                                                                      QuadraticMesh<DIM>& rQuadMesh,
                                                                                      SolidMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                                                                      std::string outputDirectory)
    : AbstractCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>(rQuadMesh,
                                                            rProblemDefinition,
                                                            outputDirectory)
{
    InitialiseContractionModels(contractionModel);
}



template<class ELASTICITY_SOLVER,unsigned DIM>
void ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::InitialiseContractionModels(ContractionModelName contractionModel)
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
        case NASH2004: //stretch dependent, will this work with explicit??
        {
            for(unsigned i=0; i<this->mTotalQuadPoints; i++)
            {
                this->mContractionModelSystems.push_back(new Nash2004ContractionModel);
            }
            break;
        }
        case KERCHOFFS2003: //stretch dependent, will this work with explicit? Answer: can be unstable
        {
            for(unsigned i=0; i<this->mTotalQuadPoints; i++)
            {
                this->mContractionModelSystems.push_back(new Kerchoffs2003ContractionModel);
            }
            break;
        }
        default:
        {
            EXCEPTION("Unknown or stretch-rate-dependent contraction model");
        }
    }

    assert(!(this->mContractionModelSystems[0]->IsStretchRateDependent()));
}

template<class ELASTICITY_SOLVER,unsigned DIM>
ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::~ExplicitCardiacMechanicsSolver()
{
    for(unsigned i=0; i<this->mContractionModelSystems.size(); i++)
    {
        delete this->mContractionModelSystems[i];
    }
}

template<class ELASTICITY_SOLVER,unsigned DIM>
void ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::GetActiveTensionAndTensionDerivs(double currentFibreStretch,
                                                                                             unsigned currentQuadPointGlobalIndex,
                                                                                             bool assembleJacobian,
                                                                                             double& rActiveTension,
                                                                                             double& rDerivActiveTensionWrtLambda,
                                                                                             double& rDerivActiveTensionWrtDLambdaDt)
{
    // the active tensions have already been computed for each contraction model, so can
    // return it straightaway..
    rActiveTension = this->mContractionModelSystems[currentQuadPointGlobalIndex]->GetActiveTension();

    // these are unset
    rDerivActiveTensionWrtLambda = 0.0;
    rDerivActiveTensionWrtDLambdaDt = 0.0;

    // store the value of given for this quad point, so that it can be used when computing
    // the active tension at the next timestep
    this->mStretches[currentQuadPointGlobalIndex] = currentFibreStretch;
}

template<class ELASTICITY_SOLVER,unsigned DIM>
void ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::Solve(double time, double nextTime, double odeTimestep)
{
    assert(time < nextTime);
    this->mCurrentTime = time;
    this->mNextTime = nextTime;
    this->mOdeTimestep = odeTimestep;

    // assemble the residual again so that mStretches is set (in GetActiveTensionAndTensionDerivs)
    // using the current deformation.
    this->AssembleSystem(true,false);

    // integrate contraction models
    for(unsigned i=0; i<this->mContractionModelSystems.size(); i++)
    {
        this->mContractionModelSystems[i]->SetStretchAndStretchRate(this->mStretches[i], 0.0 /*dlam_dt*/);
        this->mContractionModelSystems[i]->RunAndUpdate(time, nextTime, odeTimestep);
    }

    // solve
    ELASTICITY_SOLVER::Solve();
}



template class ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2>;
template class ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<3>,3>;
template class ExplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<2>,2>;
template class ExplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<3>,3>;
