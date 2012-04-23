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

#include "DiffusionForce.hpp"

template<unsigned DIM>
DiffusionForce<DIM>::DiffusionForce()
    : AbstractForce<DIM>(),
      mDiffusionConstant(0.01),
      mMechanicsCutOffLength(10.0)
{
}

template<unsigned DIM>
DiffusionForce<DIM>::~DiffusionForce()
{
}

template<unsigned DIM>
void DiffusionForce<DIM>::SetCutOffLength(double cutOffLength)
{
    assert(cutOffLength > 0.0);
    mMechanicsCutOffLength=cutOffLength;
}


template<unsigned DIM>
double DiffusionForce<DIM>::GetCutOffLength()
{
    return mMechanicsCutOffLength;
}

template<unsigned DIM>
void DiffusionForce<DIM>::SetDiffusionConstant(double diffusionConstant)
{
    mDiffusionConstant = diffusionConstant;
}

template<unsigned DIM>
double DiffusionForce<DIM>::GetDiffusionConstant()
{
    return mDiffusionConstant;
}

template<unsigned DIM>
void DiffusionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                          AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(&rCellPopulation)==NULL)
    {
        EXCEPTION("You should use the diffusion force with an OffLatticePopulation.");
    }

    AbstractOffLatticeCellPopulation<DIM>* pStaticCastCellPopulation = static_cast<AbstractOffLatticeCellPopulation<DIM>*>(&rCellPopulation);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        int node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        c_vector<double, DIM> force_contribution;
        for (unsigned i=0; i<DIM; i++)
        {
            double nu = pStaticCastCellPopulation->GetDampingConstant(node_index);
            double dt = SimulationTime::Instance()->GetTimeStep();
            double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
            force_contribution[i] = nu*sqrt(2.0*mDiffusionConstant*dt)/dt*xi;
        }
        rForces[node_index] += force_contribution;
    }
}

template<unsigned DIM>
void DiffusionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DiffusionConstant>" << mDiffusionConstant << "</DiffusionConstant> \n";
    *rParamsFile << "\t\t\t<MechanicsCutOffLength>" << mMechanicsCutOffLength << "</MechanicsCutOffLength> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DiffusionForce<1>;
template class DiffusionForce<2>;
template class DiffusionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionForce)
