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
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
DiffusionForce<DIM>::DiffusionForce()
    : AbstractForce<DIM>(),
      mDiffusionConstant(0.01),
      mAbsoluteTemperature(296.0), // default to room temperature
      mViscosity(3.204e-6), // default to viscosity of water at room temperature in (using microns and hours)
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
    assert(diffusionConstant > 0.0);
    mDiffusionConstant = diffusionConstant;
}

template<unsigned DIM>
double DiffusionForce<DIM>::GetDiffusionConstant()
{
    return mDiffusionConstant;
}

template<unsigned DIM>
void DiffusionForce<DIM>::SetAbsoluteTemperature(double newValue)
{
    assert(newValue > 0.0);
    mAbsoluteTemperature = newValue;
}

template<unsigned DIM>
double DiffusionForce<DIM>::GetAbsoluteTemperature()
{
    return mAbsoluteTemperature;
}

template<unsigned DIM>
void DiffusionForce<DIM>::SetViscosity(double newValue)
{
    assert(newValue > 0.0);
    mViscosity = newValue;
}

template<unsigned DIM>
double DiffusionForce<DIM>::GetViscosity()
{
    return mViscosity;
}

template<unsigned DIM>
void DiffusionForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                          AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    // Loop over the cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the radius, damping constant and node index associated with this cell
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        double nu = dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(&rCellPopulation)->GetDampingConstant(node_index);
        double cell_radius = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh().GetCellRadius(node_index);

        /*
         * Compute the diffusion coefficient D as D = k*T/(6*pi*eta*r), where
         *
         * k = Boltzmann's constant,
         * T = absolute temperature,
         * eta = dynamic viscosity,
         * r = cell radius.
         */
        double boltzmann_constant = 1.3806488e-23;
        double diffusion_const_scaling = boltzmann_constant*mAbsoluteTemperature/(6.0*mViscosity*M_PI);
        double diffusion_constant = diffusion_const_scaling/cell_radius;

        c_vector<double, DIM> force_contribution;
        for (unsigned i=0; i<DIM; i++)
        {
            /*
             * The force on this cell is scaled with the timestep such that when it is
             * used in the discretised equation of motion for the cell, we obtain the
             * correct formula
             *
             * x_new = x_old + sqrt(2*D*dt)*W
             *
             * where W is a standard normal random variable.
             */
            double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

            force_contribution[i] = (nu*sqrt(2.0*diffusion_constant*dt)/dt)*xi;
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
