/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "AbstractSemRandomForce.hpp"

#include "Exception.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
AbstractSemRandomForce<DIM>::AbstractSemRandomForce()
    : AbstractForce<DIM>(),
      mDiffusionConstant(0.0),
      mCoolingStartTime(0.0),
      mCoolingEndTime(0.0)
{
}

template<unsigned DIM>
double AbstractSemRandomForce<DIM>::GetDiffusionConstant() const
{
    return mDiffusionConstant;
}

template<unsigned DIM>
void AbstractSemRandomForce<DIM>::SetDiffusionConstant(double diffusionConstant)
{
    assert(diffusionConstant >= 0.0);
    mDiffusionConstant = diffusionConstant;
}

template<unsigned DIM>
void AbstractSemRandomForce<DIM>::SetCoolingWindow(double startTime, double endTime)
{
    if (endTime < startTime)
    {
        EXCEPTION("AbstractSemRandomForce: the cooling window must not end before it starts");
    }
    mCoolingStartTime = startTime;
    mCoolingEndTime = endTime;
}

template<unsigned DIM>
double AbstractSemRandomForce<DIM>::GetCurrentDiffusionConstant() const
{
    // A window of zero width means no cooling was requested, so the noise never changes
    if (mCoolingEndTime <= mCoolingStartTime)
    {
        return mDiffusionConstant;
    }

    const double time = SimulationTime::Instance()->GetTime();

    if (time <= mCoolingStartTime)
    {
        return mDiffusionConstant;
    }
    if (time >= mCoolingEndTime)
    {
        return 0.0;
    }

    const double fraction_remaining
        = (mCoolingEndTime - time) / (mCoolingEndTime - mCoolingStartTime);
    return mDiffusionConstant * fraction_remaining;
}

template<unsigned DIM>
void AbstractSemRandomForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    auto p_sem_population = dynamic_cast<SemBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (p_sem_population == nullptr)
    {
        EXCEPTION("AbstractSemRandomForce is to be used with a SemBasedCellPopulation only");
    }

    std::vector<Node<DIM>*> nodes;
    nodes.reserve(p_sem_population->GetNumNodes());
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = p_sem_population->rGetMesh().GetNodeIteratorBegin();
         node_iter != p_sem_population->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        nodes.push_back(&(*node_iter));
    }

    const std::vector<c_vector<double, DIM> > unit_noise = GenerateUnitNoise(nodes);
    assert(unit_noise.size() == nodes.size());

    const double dt = SimulationTime::Instance()->GetTimeStep();
    const double force_scale = std::sqrt(2.0 * GetCurrentDiffusionConstant() / dt);

    for (unsigned node_index = 0; node_index < nodes.size(); ++node_index)
    {
        const double damping_constant = p_sem_population->GetDampingConstant(nodes[node_index]->GetIndex());
        nodes[node_index]->AddAppliedForceContribution(damping_constant * force_scale * unit_noise[node_index]);
    }
}

template<unsigned DIM>
void AbstractSemRandomForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DiffusionConstant>" << mDiffusionConstant << "</DiffusionConstant>\n";
    *rParamsFile << "\t\t\t<CoolingStartTime>" << mCoolingStartTime << "</CoolingStartTime>\n";
    *rParamsFile << "\t\t\t<CoolingEndTime>" << mCoolingEndTime << "</CoolingEndTime>\n";

    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractSemRandomForce<1>;
template class AbstractSemRandomForce<2>;
template class AbstractSemRandomForce<3>;
