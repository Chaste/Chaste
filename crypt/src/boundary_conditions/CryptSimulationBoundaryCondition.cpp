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

#include "CryptSimulationBoundaryCondition.hpp"
#include "WntConcentration.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"

template<unsigned DIM>
CryptSimulationBoundaryCondition<DIM>::CryptSimulationBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mUseJiggledBottomCells(false)
{
}

template<unsigned DIM>
void CryptSimulationBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    // We only allow jiggling of bottom cells in 2D
    if (DIM == 1)
    {
        mUseJiggledBottomCells = false;
    }

    // Check whether a WntConcentration singleton has been set up
    bool is_wnt_included = WntConcentration<DIM>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<DIM>::Destroy();
    }

    // We iterate differently depending on whether we are using a centre- or vertex-based model
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(this->mpCellPopulation))
    {
        // Iterate over all nodes associated with real cells to update their positions
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            // Get index of node associated with cell
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

            // Get pointer to this node
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            if (!is_wnt_included)
            {
                /*
                 * If WntConcentration is not set up then stem cells must be pinned,
                 * so we reset the location of each stem cell.
                 */
                if (cell_iter->GetCellProliferativeType()->template IsType<StemCellProliferativeType>())
                {
                    // Get old node location
                    c_vector<double, DIM> old_node_location = rOldLocations.find(p_node)->second;

                    // Return node to old location
                    p_node->rGetModifiableLocation() = old_node_location;
                }
            }

            // Any cell that has moved below the bottom of the crypt must be moved back up
            if (p_node->rGetLocation()[DIM-1] < 0.0)
            {
                p_node->rGetModifiableLocation()[DIM-1] = 0.0;

                if (mUseJiggledBottomCells)
                {
                   /*
                    * Here we give the cell a push upwards so that it doesn't
                    * get stuck on the bottom of the crypt (as per #422).
                    *
                    * Note that all stem cells may get moved to the same height, so
                    * we use a random perturbation to help ensure we are not simply
                    * faced with the same problem at a different height!
                    */
                    p_node->rGetModifiableLocation()[DIM-1] = 0.05*RandomNumberGenerator::Instance()->ranf();
                }
            }
            assert(p_node->rGetLocation()[DIM-1] >= 0.0);
        }
    }
    else
    {
        // Iterate over all nodes to update their positions
        for (unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
        {
            // Get pointer to this node
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            if (!is_wnt_included)
            {
                /*
                 * If WntConcentration is not set up then stem cells must be pinned,
                 * so we reset the location of each node whose height was close to zero.
                 */
                double node_height = rOldLocations.find(p_node)->second[DIM-1];
                if (node_height < DBL_EPSILON)
                {
                    // Return node to its old height, but allow it to slide left or right
                    p_node->rGetModifiableLocation()[DIM-1] = node_height;
                }
            }

            // Any node that has moved below the bottom of the crypt must be moved back up
            if (p_node->rGetLocation()[DIM-1] < 0.0)
            {
                p_node->rGetModifiableLocation()[DIM-1] = 0.0;

                if (mUseJiggledBottomCells)
                {
                   /*
                    * Here we give the node a push upwards so that it doesn't
                    * get stuck on the bottom of the crypt.
                    */
                    p_node->rGetModifiableLocation()[DIM-1] = 0.05*RandomNumberGenerator::Instance()->ranf();
                }
            }
            assert(p_node->rGetLocation()[DIM-1] >= 0.0);
        }
    }
}

template<unsigned DIM>
bool CryptSimulationBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool boundary_condition_satisfied = true;

    /*
     * Here we verify that the boundary condition is still satisfied by simply
     * checking that no cells lies below the y=0 boundary.
     */
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Get pointer to this node
        Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

        // If this node lies below the y=0 boundary, break and return false
        if (p_node->rGetLocation()[DIM-1] < 0.0)
        {
            boundary_condition_satisfied = false;
            break;
        }
    }

    return boundary_condition_satisfied;
}

template<unsigned DIM>
void CryptSimulationBoundaryCondition<DIM>::SetUseJiggledBottomCells(bool useJiggledBottomCells)
{
    mUseJiggledBottomCells = useJiggledBottomCells;
}

template<unsigned DIM>
bool CryptSimulationBoundaryCondition<DIM>::GetUseJiggledBottomCells()
{
    return mUseJiggledBottomCells;
}

template<unsigned DIM>
void CryptSimulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<UseJiggledBottomCells>" << mUseJiggledBottomCells << "</UseJiggledBottomCells>\n";
    ///\todo Can we abstract these XML out methods and do automatic indentation?
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class CryptSimulationBoundaryCondition<1>;
template class CryptSimulationBoundaryCondition<2>;
template class CryptSimulationBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptSimulationBoundaryCondition)
