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

#include "DiffusionForce.hpp"
#include "NodeBasedCellPopulation.hpp"

//Static constant is instantiated here.
template<unsigned DIM>
const double DiffusionForce<DIM>::msBoltzmannConstant = 4.97033568e-7;

template<unsigned DIM>
DiffusionForce<DIM>::DiffusionForce()
    : AbstractForce<DIM>(),
      mAbsoluteTemperature(296.0), // default to room temperature
      mViscosity(3.204e-6) // default to viscosity of water at room temperature in (using 10 microns and hours)
{
}

template<unsigned DIM>
DiffusionForce<DIM>::~DiffusionForce()
{
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
double DiffusionForce<DIM>::GetDiffusionScalingConstant()
{
    return msBoltzmannConstant*mAbsoluteTemperature/(6.0*mViscosity*M_PI);
}

template<unsigned DIM>
void DiffusionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // Get the index, radius and damping constant of this node
        unsigned node_index = node_iter->GetIndex();
        double node_radius = node_iter->GetRadius();

        // If the node radius is zero, then it has not been set...
        if (node_radius == 0.0)
        {
            // ...so throw an exception to avoid dividing by zero when we compute diffusion_constant below
            EXCEPTION("SetRadius() must be called on each Node before calling DiffusionForce::AddForceContribution() to avoid a division by zero error");
        }

        double nu = dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(&rCellPopulation)->GetDampingConstant(node_index);

        /*
         * Compute the diffusion coefficient D as D = k*T/(6*pi*eta*r), where
         *
         * k = Boltzmann's constant,
         * T = absolute temperature,
         * eta = dynamic viscosity,
         * r = cell radius.
         */
        double diffusion_const_scaling = GetDiffusionScalingConstant();
        double diffusion_constant = diffusion_const_scaling/node_radius;

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
        node_iter->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void DiffusionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AbsoluteTemperature>" << mAbsoluteTemperature << "</AbsoluteTemperature> \n";
    *rParamsFile << "\t\t\t<Viscosity>" << mViscosity << "</Viscosity> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class DiffusionForce<1>;
template class DiffusionForce<2>;
template class DiffusionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionForce)
