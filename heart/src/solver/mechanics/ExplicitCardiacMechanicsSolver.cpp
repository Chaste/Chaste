/*

Copyright (C) University of Oxford, 2005-2012

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
