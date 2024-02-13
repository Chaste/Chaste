/*

Copyright (c) 2005-2024, University of Oxford.
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

#include "SemIntraCellularForce.hpp"

template<unsigned DIM>
SemIntraCellularForce<DIM>::SemIntraCellularForce()
   : GeneralisedLinearSpringForce<DIM>()
{
}

template<unsigned DIM>
void SemIntraCellularForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    std::cout << "AddForceContribution: " << rCellPopulation.GetNumRealCells() << std::endl;
    // Throw an exception message if not using a SemBasedCellPopulation
    if (dynamic_cast<SemBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SemIntraCellularForce is to be used with a SemBasedCellPopulation only");
    }
    
    auto cell_population = dynamic_cast<SemBasedCellPopulation<DIM>*>(&rCellPopulation);
    
    // For each cell
    for (unsigned int i = 0; i < rCellPopulation.GetNumAllCells(); i++) {
      auto p_cell = rCellPopulation.GetCellUsingLocationIndex(i);
      
      // Intra cellular force acts between all nodes in the cell
      auto p_element = cell_population->GetElementCorrespondingToCell(p_cell);
        
      auto nodes = p_element->rGetNodes();
      std::cout << nodes.size() << std::endl;

      for (unsigned index_a = 0; index_a < nodes.size(); index_a++) {
        for (unsigned index_b = index_a + 1; index_b < nodes.size(); index_b++) {
          Node<DIM>* p_node_a = nodes[index_a];
          Node<DIM>* p_node_b = nodes[index_b];
          
          // Get the node locations
          const c_vector<double, DIM>& r_node_a_location = p_node_a->rGetLocation();
          const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

          // Get the unit vector parallel to the line joining the two nodes
          c_vector<double, DIM> unit_difference;

          unit_difference = r_node_b_location - r_node_a_location;

          // Calculate the value of the rest length
          double interaction_radius = 0.2;
          double rest_length = 0.1;

          if (norm_2(unit_difference) < interaction_radius)
          {
              // Calculate the force between nodes
              c_vector<double, DIM> force = unit_difference * 0.1;
              if (norm_2(unit_difference) < rest_length) {
                force = -1.0 * force;
              }
              if (norm_2(unit_difference) < 0.02) {
                force = 20.0 * force;
              }
              c_vector<double, DIM> negative_force = -1.0 * force;

              for (unsigned j=0; j<DIM; j++)
              {
                  assert(!std::isnan(force[j]));
              }
              // Add the force contribution to each node
              p_node_a->AddAppliedForceContribution(force);
              p_node_b->AddAppliedForceContribution(negative_force);
          }
        }
      }
    }
}

template<unsigned DIM>
void SemIntraCellularForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    GeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SemIntraCellularForce<1>;
template class SemIntraCellularForce<2>;
template class SemIntraCellularForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemIntraCellularForce)
