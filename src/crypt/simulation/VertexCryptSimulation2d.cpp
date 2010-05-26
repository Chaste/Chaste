/*

Copyright (C) University of Oxford, 2005-2010

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

#include "VertexCryptSimulation2d.hpp"
#include "WntConcentration.hpp"

VertexCryptSimulation2d::VertexCryptSimulation2d(AbstractTissue<2>& rTissue,
                  std::vector<AbstractForce<2>*> forceCollection,
                  bool deleteTissueAndForceCollection,
                  bool initialiseCells)
    : TissueSimulation<2>(rTissue,
                          forceCollection,
                          deleteTissueAndForceCollection,
                          initialiseCells),
      mUseJiggledBottomCells(false)
{
    mpStaticCastTissue = static_cast<VertexBasedTissue<2>*>(&mrTissue);
}


void VertexCryptSimulation2d::UseJiggledBottomCells()
{
    mUseJiggledBottomCells = true;
}


void VertexCryptSimulation2d::WriteVisualizerSetupFile()
{
    *mpVizSetupFile << "MeshWidth\t" << mpStaticCastTissue->rGetMesh().GetWidth(0u) << "\n";
}


c_vector<double, 2> VertexCryptSimulation2d::CalculateCellDivisionVector(TissueCell& rParentCell)
{
    c_vector<double, 2> axis_of_division = zero_vector<double>(2);

    // We don't need to prescribe how 'stem' cells divide if Wnt is present
    bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<2>::Destroy();
        if (rParentCell.GetCellCycleModel()->GetCellProliferativeType() == STEM)
        {
            axis_of_division(0) = 1.0;
            axis_of_division(1) = 0.0;
        }
    }
    return axis_of_division;
}


void VertexCryptSimulation2d::ApplyTissueBoundaryConditions(const std::vector< c_vector<double, 2> >& rOldLocations)
{
    bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
    if (!is_wnt_included)
    {
        WntConcentration<2>::Destroy();
    }

    // Iterate over all nodes and update their positions according to the boundary conditions
    for (unsigned node_index=0; node_index<this->mrTissue.GetNumNodes(); node_index++)
    {
        Node<2>* p_node = this->mrTissue.GetNode(node_index);
        c_vector<double, 2> old_location = rOldLocations[node_index];

        if (!is_wnt_included)
        {
            /**
             * If WntConcentration is not set up then the stem cells must be pinned
             * to y=0, so any node whose old hieght was close to zero is moved back
             * to zero.
             */
            if (rOldLocations[node_index][1] < DBL_EPSILON)
            {
                // Return node to its old height, but allow it to slide left or right
                p_node->rGetModifiableLocation()[1] = rOldLocations[node_index][1];
            }
        }

        // Any node that has moved below the bottom of the crypt must be moved back up
        if (p_node->rGetLocation()[1] < 0.0)
        {
            p_node->rGetModifiableLocation()[1] = 0.0;
            if (this->mUseJiggledBottomCells)
            {
                // Give the node a push upwards so that it doesn't get stuck on the bottom of the crypt.
                p_node->rGetModifiableLocation()[1] = 0.05*mpRandomGenerator->ranf();
            }
        }
        assert(p_node->rGetLocation()[1] >= 0.0);
    }
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VertexCryptSimulation2d)
