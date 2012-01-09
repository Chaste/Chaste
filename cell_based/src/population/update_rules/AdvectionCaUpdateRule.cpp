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

#include "AdvectionCaUpdateRule.hpp"
#include "Exception.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
AdvectionCaUpdateRule<DIM>::AdvectionCaUpdateRule(unsigned advectionDirection, double advectionSpeed)
    : AbstractCaUpdateRule<DIM>(),
      mAdvectionDirection(advectionDirection),
      mAdvectionSpeed(advectionSpeed)
{
}

template<unsigned DIM>
AdvectionCaUpdateRule<DIM>::AdvectionCaUpdateRule()
    : AbstractCaUpdateRule<DIM>()
{
}

template<unsigned DIM>
AdvectionCaUpdateRule<DIM>::~AdvectionCaUpdateRule()
{
}

/**
 * For now this is accepting an unsigned flow input. Ultimately we want
 * to have this accepting a std::vector<unsigned> so a flow field can be
 * read in.
 */
template<unsigned DIM>
unsigned AdvectionCaUpdateRule<DIM>::GetNewLocationOfCell(unsigned currentLocationIndex,
                                                        CaBasedCellPopulation<DIM>& rCellPopulation,
                                                        double dt)
{
    // This method only works in 2D at present
    assert(DIM == 2);

    // Make sure we have a cell at this node
    if (rCellPopulation.IsEmptySite(currentLocationIndex))
    {
        EXCEPTION("There is no cell at the current location.");
    }

    double probability_of_moving = dt*mAdvectionSpeed;

    unsigned new_location_index = currentLocationIndex;

    if (RandomNumberGenerator::Instance()->ranf() < probability_of_moving)
    {
        unsigned flow_induced_new_index = currentLocationIndex;
        double width = rCellPopulation.rGetMesh().GetWidth(0);
        unsigned nodes_across = (unsigned)width + 1;
        double height = rCellPopulation.rGetMesh().GetWidth(1);
        unsigned nodes_up = (unsigned)height + 1;

        // Work out whether this node lies on any edge of the mesh
        bool on_south_edge = (currentLocationIndex < nodes_across);
        bool on_north_edge = (currentLocationIndex > nodes_up*(nodes_across - 1)-1);
        bool on_west_edge = (currentLocationIndex%nodes_across == 0);
        bool on_east_edge = (currentLocationIndex%nodes_across == nodes_across - 1);

        switch (mAdvectionDirection)
        {
            case 0:
            {
                if (!on_north_edge)
                {
                    flow_induced_new_index += nodes_across;
                }
                break;
            }
            case 1:
            {
                if (!on_north_edge && !on_west_edge)
                {
                    flow_induced_new_index += nodes_across - 1;
                }
                break;
            }
            case 2:
            {
                if (!on_west_edge)
                {
                    flow_induced_new_index -= 1;
                }
                break;
            }
            case 3:
            {
                if (!on_south_edge && !on_west_edge)
                {
                    flow_induced_new_index -= nodes_across + 1;
                }
                break;
            }
            case 4:
            {
                if (!on_south_edge)
                {
                    flow_induced_new_index -= nodes_across;
                }
                break;
            }
            case 5:
            {
                if (!on_south_edge && !on_east_edge)
                {
                    flow_induced_new_index -= nodes_across - 1;
                }
                break;
            }
            case 6:
            {
                if (!on_east_edge)
                {
                    flow_induced_new_index += 1;
                }
                break;
            }
            case 7:
            {
                if (!on_north_edge && !on_east_edge)
                {
                    flow_induced_new_index += nodes_across + 1;
                }
                break;
            }
            default:
                NEVER_REACHED;
        }

        if (rCellPopulation.IsEmptySite(flow_induced_new_index))
        {
            new_location_index = flow_induced_new_index;
        }
    }
    return new_location_index;
}

template<unsigned DIM>
unsigned AdvectionCaUpdateRule<DIM>::GetAdvectionDirection()
{
    return mAdvectionDirection;
}

template<unsigned DIM>
double AdvectionCaUpdateRule<DIM>::GetAdvectionSpeed()
{
    return mAdvectionSpeed;
}

template<unsigned DIM>
void AdvectionCaUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AdvectionDirection>" << mAdvectionDirection << "</AdvectionDirection>\n";
    *rParamsFile << "\t\t\t<AdvectionSpeed>" << mAdvectionSpeed << "</AdvectionSpeed>\n";

    // Call method on direct parent class
    AbstractCaUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AdvectionCaUpdateRule<1>;
template class AdvectionCaUpdateRule<2>;
template class AdvectionCaUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AdvectionCaUpdateRule)
