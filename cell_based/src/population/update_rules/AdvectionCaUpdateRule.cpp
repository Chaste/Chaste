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
