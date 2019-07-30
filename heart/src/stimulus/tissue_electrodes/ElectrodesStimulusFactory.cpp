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

#include "ElectrodesStimulusFactory.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "IsNan.hpp"
#include "HeartConfig.hpp"
#include "GaussianQuadratureRule.hpp"
#include "RegularStimulus.hpp"

template<unsigned DIM>
ElectrodesStimulusFactory<DIM>::ElectrodesStimulusFactory(std::vector<std::pair<AbstractChasteRegion<DIM>*, AbstractChasteRegion<DIM>*> >& rElectrodePairs,
                                                          std::vector<double>& rStimulusMagnitudes,
                                                          std::vector<double>& rDurations,
                                                          std::vector<double>& rPeriods,
                                                          std::vector<double>& rStarts,
                                                          std::vector<double>& rEnds)
    : mrElectrodePairs(rElectrodePairs),
      mrMagnitudes(rStimulusMagnitudes),
      mrDurations(rDurations),
      mrPeriods(rPeriods),
      mrStarts(rStarts),
      mrEnds(rEnds),
      mGroundSecondElectrode(false)
{
    if ((rElectrodePairs.size() != rStimulusMagnitudes.size()) ||
        (rElectrodePairs.size() != rDurations.size()) ||
        (rElectrodePairs.size() != rPeriods.size()) ||
        (rElectrodePairs.size() != rStarts.size()) ||
        (rElectrodePairs.size() != rEnds.size()))
    {
        EXCEPTION ("Vector of electrode pairs and vector of stimulation paremeters must have the same size");
    }

    mMagnitudesElectrode1 = mrMagnitudes;
    mMagnitudesElectrode2 = mrMagnitudes;
}

template<unsigned DIM>
ElectrodesStimulusFactory<DIM>::~ElectrodesStimulusFactory()
{
}

template<unsigned DIM>
void ElectrodesStimulusFactory<DIM>::CheckForElectrodesIntersection()
{
    std::vector<unsigned> nodes_in_all_electrodes;
    for (unsigned global_node_index = 0; global_node_index < this->mpMesh->GetNumNodes(); global_node_index++)
    {
        if (this->mpMesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(global_node_index))
        {
            for (unsigned pair_index = 0; pair_index <mrElectrodePairs.size(); pair_index++)
            {
                if (mrElectrodePairs[pair_index].first->DoesContain( this->mpMesh->GetNode(global_node_index)->GetPoint()))
                {
                    nodes_in_all_electrodes.push_back( this->mpMesh->GetNode(global_node_index)->GetIndex() );
                }
                if (mrElectrodePairs[pair_index].second->DoesContain( this->mpMesh->GetNode(global_node_index)->GetPoint()))
                {
                    nodes_in_all_electrodes.push_back( this->mpMesh->GetNode(global_node_index)->GetIndex() );
                }
            }
        }
    }
    PetscTools::Barrier();
    for (unsigned node_index = 0; node_index < nodes_in_all_electrodes.size(); node_index++)
    {
        unsigned number_of_hits = 0;
        for (unsigned node_to_check = 0; node_to_check < nodes_in_all_electrodes.size(); node_to_check++)
        {
            if (nodes_in_all_electrodes[node_index] == nodes_in_all_electrodes[node_to_check] )
            {
                number_of_hits++;
            }
        }
        if (number_of_hits>1)
        {
            EXCEPTION("Two or more electrodes intersect with each other");
        }
    }
}

template<unsigned DIM>
void ElectrodesStimulusFactory<DIM>::GroundSecondElectrode(bool grounded)
{
    mGroundSecondElectrode = grounded;
}

template<unsigned DIM>
void ElectrodesStimulusFactory<DIM>::SetCompatibleExtracellularStimulus()
{
    assert(this->mpMesh!=NULL);
    try
    {
        CheckForElectrodesIntersection();
    }
    catch (Exception &e)
    {
        PetscTools::ReplicateException(true); //Electrodes intersect
        throw e;
    }
    PetscTools::ReplicateException(false);

    for (unsigned pair_index = 0; pair_index < mrElectrodePairs.size(); pair_index++)
    {

        if (!mGroundSecondElectrode)//no grounding of second electrode
        {
            //compute the two fluxes for this pair
            double flux_electrode_1  = ComputeElectrodeTotalFlux(mrElectrodePairs[pair_index].first, mMagnitudesElectrode1[pair_index]);
            double flux_electrode_2  = ComputeElectrodeTotalFlux(mrElectrodePairs[pair_index].second, mMagnitudesElectrode2[pair_index]);

            //Scale the magnitude of the second electrode for this pair
            mMagnitudesElectrode2[pair_index] = -mMagnitudesElectrode1[pair_index]*flux_electrode_1/flux_electrode_2;

            // Paranoia.
            assert( flux_electrode_2 != 0.0);
            assert( ! std::isnan(mMagnitudesElectrode2[pair_index]));

        }
        else//second electrode is grounded
        {
            this->mGroundedRegions.push_back( mrElectrodePairs[pair_index].second );
        }
    }
}

template<unsigned DIM>
boost::shared_ptr<AbstractStimulusFunction> ElectrodesStimulusFactory<DIM>::CreateStimulusForNode(Node<DIM>* pNode)
{
    boost::shared_ptr<RegularStimulus> p_stimulus;
    for (unsigned pair_index = 0; pair_index < mrElectrodePairs.size(); pair_index++)
    {
        if (mrElectrodePairs[pair_index].first->DoesContain(pNode->GetPoint()) )
        {
            p_stimulus.reset( new RegularStimulus(mMagnitudesElectrode1[pair_index], mrDurations[pair_index], mrPeriods[pair_index], mrStarts[pair_index], mrEnds[pair_index]));

        }
        else if (mrElectrodePairs[pair_index].second->DoesContain(pNode->GetPoint()) )
        {
            p_stimulus.reset ( new RegularStimulus(mMagnitudesElectrode2[pair_index], mrDurations[pair_index], mrPeriods[pair_index], mrStarts[pair_index], mrEnds[pair_index]));

        }
        else//no stimulus here
        {
            p_stimulus.reset ( new RegularStimulus(0.0, mrDurations[pair_index], mrPeriods[pair_index], mrStarts[pair_index], mrEnds[pair_index]));
        }
    }
    return p_stimulus;
}

template<unsigned DIM>
double ElectrodesStimulusFactory<DIM>::ComputeElectrodeTotalFlux(AbstractChasteRegion<DIM>* pRegion, double stimulusMagnitude)
{
    GaussianQuadratureRule<DIM>* pQuadRule = new GaussianQuadratureRule<DIM>(2);

    //the basis functions
    c_vector<double, DIM+1> phi;

    double total_electrode_flux = 0.0;
    double ret;

    for (typename AbstractTetrahedralMesh<DIM,DIM>::NodeIterator node_iter=this->mpMesh->GetNodeIteratorBegin();
         node_iter != this->mpMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        if (pRegion->DoesContain((*node_iter).GetPoint()))
        {
            unsigned node_index = node_iter->GetIndex();
            assert(node_index < this->mpMesh->GetNumNodes());
            double contribution_of_this_node = 0.0;
            // Loop over the elements where this node is contained
            for (std::set<unsigned>::iterator iter = this->mpMesh->GetNode(node_index)->rGetContainingElementIndices().begin();
                 iter != this->mpMesh->GetNode(node_index)->rGetContainingElementIndices().end();
                 ++iter)
            {
                Element<DIM,DIM>* p_element = this->mpMesh->GetElement(*iter);

                /*Determine jacobian for this element*/
                c_matrix<double, DIM, DIM> jacobian;
                c_matrix<double, DIM, DIM> inverse_jacobian;//unused here, but needed for the function below
                double jacobian_determinant;
                this->mpMesh->GetInverseJacobianForElement(p_element->GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

                double contribution_of_this_element = 0.0;//...to this node
                 // loop over Gauss points
                for (unsigned quad_index=0; quad_index < pQuadRule->GetNumQuadPoints(); quad_index++)
                {
                    const ChastePoint<DIM>& quad_point = pQuadRule->rGetQuadPoint(quad_index);
                    BasisFunction::ComputeBasisFunctions(quad_point, phi);

                    double interpolated_stimulus = 0.0;
                    //loop over nodes in this element
                    for (unsigned node_index_in_element = 0; node_index_in_element < p_element->GetNumNodes(); node_index_in_element++)
                    {
                        //const Node<DIM>* p_node = p_element->GetNode(node_index_in_element);
                        assert(p_element->GetNumNodes() == DIM+1);
                        interpolated_stimulus += stimulusMagnitude*phi(node_index_in_element);
                        contribution_of_this_element += interpolated_stimulus * phi(node_index_in_element) * jacobian_determinant * pQuadRule->GetWeight(quad_index);
                    }

                }/*end of loop over gauss points*/
                contribution_of_this_node += contribution_of_this_element;

            }/*end of loop over elements where the node is contained*/
            total_electrode_flux += contribution_of_this_node;
        }/* end of if that checks if node is in the electrode*/
    }/* end of loop over nodes in the mesh*/
#ifndef NDEBUG
    int mpi_ret = MPI_Allreduce(&total_electrode_flux, &ret, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    assert(mpi_ret == MPI_SUCCESS);
#else
    MPI_Allreduce(&total_electrode_flux, &ret, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
#endif

    //clear up memory
    delete pQuadRule;

    assert(ret < DBL_MAX);
    return ret;
}

// Explicit instantiation
template class ElectrodesStimulusFactory<1>;
template class ElectrodesStimulusFactory<2>;
template class ElectrodesStimulusFactory<3>;
