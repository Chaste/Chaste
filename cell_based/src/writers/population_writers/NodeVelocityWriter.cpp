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

#include "NodeVelocityWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NodeVelocityWriter<ELEMENT_DIM, SPACE_DIM>::NodeVelocityWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("nodevelocities.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodeVelocityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // We should never encounter deleted nodes when calling this method
        assert(!node_iter->IsDeleted());

        unsigned node_index = node_iter->GetIndex();

        // Check that results should be written for this node
        bool is_real_node = !(pCellPopulation->IsGhostNode(node_index));

        // We should never encounter nodes associated with dead cells when calling this method
        if (is_real_node)
        {
            assert(!pCellPopulation->GetCellUsingLocationIndex(node_index)->IsDead());
        }

        if (is_real_node)
        {
            // Write this node's index to file
            *this->mpOutStream << node_index  << " ";

            // Write this node's position to file
            const c_vector<double, SPACE_DIM>& position = node_iter->rGetLocation();
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                *this->mpOutStream << position[i] << " ";
            }

            // Write this node's velocity to file
            double time_step = SimulationTime::Instance()->GetTimeStep(); ///\todo correct time step? (#2404)
            double damping_constant = pCellPopulation->GetDampingConstant(node_index);
            c_vector<double, SPACE_DIM> velocity = time_step * node_iter->rGetAppliedForce() / damping_constant;
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                *this->mpOutStream << velocity[i] << " ";
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodeVelocityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("NodeVelocityWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodeVelocityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // We should never encounter deleted nodes when calling this method
        assert(!node_iter->IsDeleted());

        unsigned node_index = node_iter->GetIndex();

        // Check that results should be written for this node
        bool is_real_node = !(pCellPopulation->IsGhostNode(node_index));

        // We should never encounter nodes associated with dead cells when calling this method
        if (is_real_node)
        {
            assert(!pCellPopulation->GetCellUsingLocationIndex(node_index)->IsDead());
        }

        if (is_real_node)
        {
            // Write this node's index to file
            *this->mpOutStream << node_index  << " ";

            // Write this node's position to file
            const c_vector<double, SPACE_DIM>& position = node_iter->rGetLocation();
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                *this->mpOutStream << position[i] << " ";
            }

            // Write this node's velocity to file
            double time_step = SimulationTime::Instance()->GetTimeStep(); ///\todo correct time step? (#2404)
            double damping_constant = pCellPopulation->GetDampingConstant(node_index);
            c_vector<double, SPACE_DIM> velocity = time_step * node_iter->rGetAppliedForce() / damping_constant;
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                *this->mpOutStream << velocity[i] << " ";
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodeVelocityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("NodeVelocityWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodeVelocityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    for (typename AbstractMesh<SPACE_DIM, SPACE_DIM>::NodeIterator node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
         node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // We should never encounter deleted nodes when calling this method
        assert(!node_iter->IsDeleted());

        // Write this node's index to file
        unsigned node_index = node_iter->GetIndex();
        *this->mpOutStream << node_index  << " ";

        // Write this node's position to file
        const c_vector<double, SPACE_DIM>& position = node_iter->rGetLocation();
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << position[i] << " ";
        }

        // Write this node's velocity to file
        double time_step = SimulationTime::Instance()->GetTimeStep(); ///\todo correct time step? (#2404)
        double damping_constant = pCellPopulation->GetDampingConstant(node_index);
        c_vector<double, SPACE_DIM> velocity = time_step * node_iter->rGetAppliedForce() / damping_constant;
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << velocity[i] << " ";
        }
    }
}

// Explicit instantiation
template class NodeVelocityWriter<1,1>;
template class NodeVelocityWriter<1,2>;
template class NodeVelocityWriter<2,2>;
template class NodeVelocityWriter<1,3>;
template class NodeVelocityWriter<2,3>;
template class NodeVelocityWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NodeVelocityWriter)
