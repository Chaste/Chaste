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

#include "AirwayRemesher.hpp"

AirwayRemesher::AirwayRemesher(TetrahedralMesh<1,3>& rMesh,
                               unsigned rootIndex) : mrMesh(rMesh),
                                                     mOutletNodeIndex(rootIndex),
                                                     mWalker(mrMesh, mOutletNodeIndex),
                                                     mCalculator(mrMesh, mOutletNodeIndex)
{
}

void AirwayRemesher::Remesh(MutableMesh<1,3>& rOutputMesh)
{
    Remesh(rOutputMesh, DBL_MAX);
}

void AirwayRemesher::Remesh(MutableMesh<1,3>& outputMesh, double maximumResistance)
{
  std::vector<AirwayBranch*> branches = mCalculator.GetBranches();
  std::map<Node<3>*, Node<3>*> node_map;

  for (std::vector<AirwayBranch*>::iterator iter = branches.begin(); iter != branches.end(); ++iter)
  {
      double resistance = (*iter)->GetPoiseuilleResistance();

      Node<3>* start_node = (*iter)->GetProximalNode();
      Node<3>* distal_node = (*iter)->GetDistalNode();
      c_vector<double, 3> direction = distal_node->rGetLocation() - start_node->rGetLocation();

      double equivalent_radius = std::pow(norm_2(direction)/resistance, 1.0/4.0);
      unsigned subdivisions = (unsigned)(ceil(resistance/maximumResistance));

      Node<3>* new_start_node;
      if (node_map.count(start_node) == 0)
      {
          new_start_node = new Node<3>(0, start_node->rGetLocation(), start_node->IsBoundaryNode());
          outputMesh.AddNode(new_start_node);
          node_map[start_node] = new_start_node;
      }
      else
      {
          new_start_node = node_map[start_node];
      }

      for (unsigned division = 0; division < (subdivisions - 1); ++division) //Note we handle the final element as a special case
      {
          Node<3>* new_end_node = new Node<3>(0, start_node->rGetLocation() + ((double)division + 1)/subdivisions*direction);

          outputMesh.AddNode(new_end_node);

          std::vector<Node<3>* > nodes;
          nodes.push_back(new_start_node);
          nodes.push_back(new_end_node);

          Element<1,3>* elem = new Element<1,3>(0, nodes);
          elem->SetAttribute(equivalent_radius);
          outputMesh.AddElement(elem);

          new_start_node = new_end_node;
      }

      Node<3>* new_distal_node = new Node<3>(0, distal_node->rGetLocation(), distal_node->IsBoundaryNode());
      outputMesh.AddNode(new_distal_node);
      node_map[distal_node] = new_distal_node;

      std::vector<Node<3>* > nodes;
      nodes.push_back(new_start_node);
      nodes.push_back(new_distal_node);

      Element<1,3>* elem = new Element<1,3>(0, nodes);
      elem->SetAttribute(equivalent_radius);
      outputMesh.AddElement(elem);
  }

  // If average resistance is > maximumResistance, sub divide
}
