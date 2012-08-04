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
#include "FakeBathContractionModel.hpp"

template<class ELASTICITY_SOLVER,unsigned DIM>
ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::ExplicitCardiacMechanicsSolver(ContractionModelName contractionModelName,
                                                                                      QuadraticMesh<DIM>& rQuadMesh,
                                                                                      ElectroMechanicsProblemDefinition<DIM>& rProblemDefinition,
                                                                                      std::string outputDirectory)
    : AbstractCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>(rQuadMesh,
                                                            contractionModelName,
                                                            rProblemDefinition,
                                                            outputDirectory)
{

}



template<class ELASTICITY_SOLVER,unsigned DIM>
void ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::InitialiseContractionModels(ContractionModelName contractionModelName)
{
    for(std::map<unsigned,DataAtQuadraturePoint>::iterator iter = this->mQuadPointToDataAtQuadPointMap.begin();
        iter != this->mQuadPointToDataAtQuadPointMap.end();
        iter++)
    {
        AbstractContractionModel* p_contraction_model;

        if (iter->second.Active == true)
        {
            //tissue node
            switch(contractionModelName)
            {
                case NONPHYSIOL1:
                case NONPHYSIOL2:
                case NONPHYSIOL3:
                {
                    unsigned option = (contractionModelName==NONPHYSIOL1 ? 1 : (contractionModelName==NONPHYSIOL2? 2 : 3));
                    p_contraction_model = new NonPhysiologicalContractionModel(option);
                    break;
                }
                case NASH2004: //stretch dependent, will this work with explicit??
                {
                    p_contraction_model = new Nash2004ContractionModel;
                    break;
                }
                case KERCHOFFS2003: //stretch dependent, will this work with explicit? Answer: can be unstable
                {
                    p_contraction_model = new Kerchoffs2003ContractionModel;
                    break;
                }
                default:
                {
                    EXCEPTION("Unknown or stretch-rate-dependent contraction model");
                }
            }

            iter->second.ContractionModel = p_contraction_model;
        }
        else
        {
            //bath
            p_contraction_model = new FakeBathContractionModel;
            iter->second.ContractionModel = p_contraction_model;
        }
    }
}

template<class ELASTICITY_SOLVER,unsigned DIM>
ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::~ExplicitCardiacMechanicsSolver()
{
}

template<class ELASTICITY_SOLVER,unsigned DIM>
void ExplicitCardiacMechanicsSolver<ELASTICITY_SOLVER,DIM>::GetActiveTensionAndTensionDerivs(double currentFibreStretch,
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

    // the active tensions have already been computed for each contraction model, so can
    // return it straightaway..
    rActiveTension = r_data_at_quad_point.ContractionModel->GetActiveTension();

    // these are unset
    rDerivActiveTensionWrtLambda = 0.0;
    rDerivActiveTensionWrtDLambdaDt = 0.0;

    // store the value of given for this quad point, so that it can be used when computing
    // the active tension at the next timestep
    r_data_at_quad_point.Stretch = currentFibreStretch;

    // increment the iterator
    this->mMapIterator++;
    if(this->mMapIterator==this->mQuadPointToDataAtQuadPointMap.end())
    {
        this->mMapIterator = this->mQuadPointToDataAtQuadPointMap.begin();
    }

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
    for(std::map<unsigned,DataAtQuadraturePoint>::iterator iter = this->mQuadPointToDataAtQuadPointMap.begin();
        iter != this->mQuadPointToDataAtQuadPointMap.end();
        iter++)
    {
        AbstractContractionModel* p_contraction_model = iter->second.ContractionModel;
        double stretch = iter->second.Stretch;
        p_contraction_model->SetStretchAndStretchRate(stretch, 0.0 /*dlam_dt*/);
        p_contraction_model->RunAndUpdate(time, nextTime, odeTimestep);
    }

    // solve
    ELASTICITY_SOLVER::Solve();
}



template class ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<2>,2>;
template class ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<3>,3>;
template class ExplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<2>,2>;
template class ExplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<3>,3>;
