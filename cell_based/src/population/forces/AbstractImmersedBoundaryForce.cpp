/*

Copyright (c) 2005-2025, University of Oxford.
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

#include <random>

#include "AbstractImmersedBoundaryForce.hpp"

template <unsigned DIM>
AbstractImmersedBoundaryForce<DIM>::AbstractImmersedBoundaryForce()
    : mAdditiveNormalNoise(false),
      mNormalNoiseMean(1.0),
      mNormalNoiseStdDev(0.0)
{
}

template<unsigned DIM>
AbstractImmersedBoundaryForce<DIM>::~AbstractImmersedBoundaryForce()
{
}

template<unsigned DIM>
void AbstractImmersedBoundaryForce<DIM>::AddNormalNoiseToNodes(
    ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    // Add random forces to elements
    for (unsigned elem_idx = 0; elem_idx < rCellPopulation.GetNumElements(); ++elem_idx)
    {
        unsigned num_nodes_this_elem = rCellPopulation.GetElement(elem_idx)->GetNumNodes();

        // Generate a vector with 0, 1, 2, ..., num_nodes-1
        std::vector<unsigned> random_node_order;
        random_node_order.reserve(num_nodes_this_elem);
        for (unsigned node_idx = 0; node_idx < num_nodes_this_elem; ++node_idx)
        {
            random_node_order.push_back(node_idx);
        }

        // Shuffle the vector
        std::shuffle(random_node_order.begin(), random_node_order.end(), std::mt19937(std::random_device()()));

        // Only go to half way (rounded down - forget about the very last node for now)
        unsigned half_way = num_nodes_this_elem / 2;
        for (unsigned i = 0; i < half_way; ++i)
        {
            unsigned rand_node_a = random_node_order[2 * i];
            unsigned rand_node_b = random_node_order[2 * i + 1];

            c_vector<double, DIM> random_force;
            for (unsigned dim = 0; dim < DIM; ++dim)
            {
                random_force[dim] = p_gen->NormalRandomDeviate(mNormalNoiseMean, mNormalNoiseStdDev);
            }

            Node<DIM>* p_node_a = rCellPopulation.GetElement(elem_idx)->GetNode(rand_node_a);
            Node<DIM>* p_node_b = rCellPopulation.GetElement(elem_idx)->GetNode(rand_node_b);

            double avg_magnitude = 0.5 * (norm_2(p_node_a->rGetAppliedForce()) +
                                          norm_2(p_node_b->rGetAppliedForce()));

            random_force *= avg_magnitude;
            p_node_a->AddAppliedForceContribution(random_force);

            random_force *= -1.0;
            p_node_b->AddAppliedForceContribution(random_force);
        }
    }

    // Do the same for laminas
    for (unsigned lam_idx = 0; lam_idx < rCellPopulation.GetNumLaminas(); ++lam_idx)
    {
        unsigned num_nodes_this_elem = rCellPopulation.GetLamina(lam_idx)->GetNumNodes();

        // Generate a vector with 0, 1, 2, ..., num_nodes-1
        std::vector<unsigned> random_node_order;
        random_node_order.reserve(num_nodes_this_elem);
        for (unsigned node_idx = 0; node_idx < num_nodes_this_elem; ++node_idx)
        {
            random_node_order.push_back(node_idx);
        }

        // Shuffle the vector
        std::shuffle(random_node_order.begin(), random_node_order.end(), std::mt19937(std::random_device()()));

        // Only go to half way (rounded down - forget about the very last node for now)
        unsigned half_way = num_nodes_this_elem / 2;
        for (unsigned i = 0; i < half_way; ++i)
        {
            unsigned rand_node_a = random_node_order[2 * i];
            unsigned rand_node_b = random_node_order[2 * i + 1];

            c_vector<double, DIM> random_force;
            for (unsigned dim = 0; dim < DIM; ++dim)
            {
                random_force[dim] = p_gen->NormalRandomDeviate(mNormalNoiseMean, mNormalNoiseStdDev);
            }

            Node<DIM>* p_node_a = rCellPopulation.GetLamina(lam_idx)->GetNode(rand_node_a);
            Node<DIM>* p_node_b = rCellPopulation.GetLamina(lam_idx)->GetNode(rand_node_b);

            double avg_magnitude = 0.5 * (norm_2(p_node_a->rGetAppliedForce()) +
                                          norm_2(p_node_b->rGetAppliedForce()));

            random_force *= avg_magnitude;
            p_node_a->AddAppliedForceContribution(random_force);

            random_force *= -1.0;
            p_node_b->AddAppliedForceContribution(random_force);
        }
    }
}

template <unsigned DIM>
bool AbstractImmersedBoundaryForce<DIM>::GetAdditiveNormalNoise() const
{
    return mAdditiveNormalNoise;
}

template <unsigned DIM>
void AbstractImmersedBoundaryForce<DIM>::SetAdditiveNormalNoise(
    bool additiveNormalNoise)
{
    mAdditiveNormalNoise = additiveNormalNoise;
}

template <unsigned DIM>
double AbstractImmersedBoundaryForce<DIM>::GetNormalNoiseMean() const
{
    return mNormalNoiseMean;
}

template <unsigned DIM>
void AbstractImmersedBoundaryForce<DIM>::SetNormalNoiseMean(
    double normalNoiseMean)
{
    mNormalNoiseMean = normalNoiseMean;
}

template <unsigned DIM>
double AbstractImmersedBoundaryForce<DIM>::GetNormalNoiseStdDev() const
{
    return mNormalNoiseStdDev;
}

template <unsigned DIM>
void AbstractImmersedBoundaryForce<DIM>::SetNormalNoiseStdDev(
    double normalNoiseStdDev)
{
    mNormalNoiseStdDev = normalNoiseStdDev;
}

template<unsigned DIM>
void AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(
    out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AdditiveNormalNoise>" << mAdditiveNormalNoise << "</AdditiveNormalNoise>\n";
    *rParamsFile << "\t\t\t<NormalNoiseMean>" << mNormalNoiseMean << "</NormalNoiseMean>\n";
    *rParamsFile << "\t\t\t<NormalNoiseStdDev>" << mNormalNoiseStdDev << "</NormalNoiseStdDev>\n";
}

// Explicit instantiation
template class AbstractImmersedBoundaryForce<1>;
template class AbstractImmersedBoundaryForce<2>;
template class AbstractImmersedBoundaryForce<3>;
