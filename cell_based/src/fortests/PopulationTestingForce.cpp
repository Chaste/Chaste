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

#include "PopulationTestingForce.hpp"

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
PopulationTestingForce<ELEMENT_DIM, SPACE_DIM>::PopulationTestingForce(bool hasPositionDependence)
    :AbstractForce<ELEMENT_DIM, SPACE_DIM>(),
     mWithPositionDependence(hasPositionDependence)
{
}

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationTestingForce<ELEMENT_DIM, SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    for (unsigned i=0; i<rCellPopulation.GetNumNodes(); i++)
    {
        c_vector<double, SPACE_DIM> force;

        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            if (mWithPositionDependence)
            {
              force[j] = (j+1)*i*0.01*rCellPopulation.GetNode(i)->rGetLocation()[j];
            }
            else
            {
                force[j] = (j+1)*i*0.01;
            }
        }
        rCellPopulation.GetNode(i)->ClearAppliedForce();
        rCellPopulation.GetNode(i)->AddAppliedForceContribution(force);
    }
}

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PopulationTestingForce<ELEMENT_DIM, SPACE_DIM>::GetExpectedOneStepLocationFE(unsigned nodeIndex,
                                                                                       double damping,
                                                                                       c_vector<double, SPACE_DIM>& oldLocation,
                                                                                       double dt)
{
    c_vector<double, SPACE_DIM> result;
    for (unsigned j = 0; j < SPACE_DIM; j++)
    {
        if (mWithPositionDependence)
        {
            result[j] = oldLocation[j] + dt * (j+1)*0.01*nodeIndex * oldLocation[j] / damping;
        }
        else
        {
            result[j] = oldLocation[j] + dt * (j+1)*0.01*nodeIndex/damping;
        }
    }
    return result;
}

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PopulationTestingForce<ELEMENT_DIM, SPACE_DIM>::GetExpectedOneStepLocationRK4(unsigned nodeIndex,
                                                                                      double damping,
                                                                                      c_vector<double, SPACE_DIM>& oldLocation,
                                                                                      double dt)
{
    c_vector<double, SPACE_DIM> result;
    for (unsigned j = 0; j < SPACE_DIM; j++)
    {
        double k1 = (j+1)*0.01*nodeIndex * oldLocation[j] / damping;
        double k2 = (j+1)*0.01*nodeIndex * (oldLocation[j] + (dt/2.0)*k1) / damping;
        double k3 = (j+1)*0.01*nodeIndex * (oldLocation[j] + (dt/2.0)*k2) / damping;
        double k4 = (j+1)*0.01*nodeIndex * (oldLocation[j] + dt*k3) / damping;
        result[j] = oldLocation[j] + (1.0/6.0)*dt*(k1 + 2*k2 + 2*k3 + k4);
    }
    return result;
}

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PopulationTestingForce<ELEMENT_DIM, SPACE_DIM>::GetExpectedOneStepLocationAM2(unsigned nodeIndex,
                                                                                    double damping,
                                                                                    c_vector<double, SPACE_DIM>& oldLocation,
                                                                                    double dt)
{
    c_vector<double, SPACE_DIM> result;
    for (unsigned j = 0; j < SPACE_DIM; j++)
    {
        result[j] = oldLocation[j] * (1 + 0.5*dt*0.01*(j+1)*nodeIndex / damping) / (1 - 0.5*dt*0.01*(j+1)*nodeIndex / damping);
    }
    return result;
}

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PopulationTestingForce<ELEMENT_DIM, SPACE_DIM>::GetExpectedOneStepLocationBE(unsigned nodeIndex,
                                                                                      double damping,
                                                                                      c_vector<double, SPACE_DIM>& oldLocation,
                                                                                      double dt)
{
    c_vector<double, SPACE_DIM> result;
    for (unsigned j = 0; j < SPACE_DIM; j++)
    {
        result[j] = oldLocation[j] / (1 - dt*0.01*(j+1)*nodeIndex/damping);
    }
    return result;
}

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationTestingForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // This force is not used in a simulation so this method is never called
    NEVER_REACHED;
}

// Explicit instantiation
template class PopulationTestingForce<1,1>;
template class PopulationTestingForce<1,2>;
template class PopulationTestingForce<2,2>;
template class PopulationTestingForce<1,3>;
template class PopulationTestingForce<2,3>;
template class PopulationTestingForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PopulationTestingForce)
