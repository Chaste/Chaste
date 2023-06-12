/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "ExtrinsicPullModifier.hpp"

template<unsigned DIM>
ExtrinsicPullModifier<DIM>::ExtrinsicPullModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mApplyExtrinsicPullToAllNodes(true),
      mSpeed(1.0)
{
}

template<unsigned DIM>
ExtrinsicPullModifier<DIM>::~ExtrinsicPullModifier()
{
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double epsilon = 0.8;

    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    ChasteCuboid<DIM> bounds = rCellPopulation.rGetMesh().CalculateBoundingBox();
    double x_min = bounds.rGetLowerCorner()[0];
    double x_max = bounds.rGetUpperCorner()[0];

    if (mApplyExtrinsicPullToAllNodes)
    {
        // Pull on all nodes, with a constant strain rate
        double width = x_max - x_min;
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            double speed = mSpeed * (p_node->rGetLocation()[0] - x_min) / width;

            if (p_node->rGetLocation()[0] > x_min + epsilon)
            {
                    p_node->rGetModifiableLocation()[0] += speed*dt;
            }
        }
    }
    else
    {
        // Pull on the right-most nodes only, with a constant speed
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            if (fabs(p_node->rGetLocation()[0] - x_max) < 0.1)
            {
                p_node->rGetModifiableLocation()[0] += mSpeed*dt;
            }
        }
    }
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::SetApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes)
{
    mApplyExtrinsicPullToAllNodes = applyExtrinsicPullToAllNodes;
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::SetSpeed(double speed)
{
    mSpeed = speed;
}

template<unsigned DIM>
bool ExtrinsicPullModifier<DIM>::GetApplyExtrinsicPullToAllNodes()
{
    return mApplyExtrinsicPullToAllNodes;
}

template<unsigned DIM>
double ExtrinsicPullModifier<DIM>::GetSpeed()
{
    return mSpeed;
}

template<unsigned DIM>
void ExtrinsicPullModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ApplyExtrinsicPullToAllNodes>" << mApplyExtrinsicPullToAllNodes << "</ApplyExtrinsicPullToAllNodes>\n";
    *rParamsFile << "\t\t\t<Speed>" << mSpeed << "</Speed>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ExtrinsicPullModifier<1>;
template class ExtrinsicPullModifier<2>;
template class ExtrinsicPullModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtrinsicPullModifier)