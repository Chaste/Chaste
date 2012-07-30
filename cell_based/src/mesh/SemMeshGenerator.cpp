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

#include "SemMeshGenerator.hpp"

template<unsigned DIM>
SemMeshGenerator<DIM>::SemMeshGenerator(unsigned numCellsAcross,
                                            unsigned numCellsUp,
                                            unsigned numCellsDeep,
                                            unsigned numSubCellularElementsPerCellAcross,
                                            unsigned numSubCellularElementsPerCellUp,
                                            unsigned numSubCellularElementsPerCellDeep)
{
    mTotalSubcellularElementsPerCell = numSubCellularElementsPerCellAcross + numSubCellularElementsPerCellUp + numSubCellularElementsPerCellDeep;

    // Check input is of correct type.
    if(numCellsAcross == 0)
    {
        EXCEPTION("An SemMeshGenerator must have numCellsAcross > 0");
    }
    if(numSubCellularElementsPerCellAcross == 0)
    {
        EXCEPTION("An SemMeshGenerator must have more than 1 subcellular element per cell. Set numSubCellularElementsPerCell > 0");
    }

    // Make sure the dimensions match the number of cells up / deep.
    if(DIM == 1)
    {
        EXCEPTION("Trying to create a 1D mesh, SemMesh only defined for 2 and 3D");
    }
    if(numCellsDeep > 1 && DIM < 3)
    {
        EXCEPTION("Trying to create a 3D SemMesh with DIM < 3");
    }

    mTotalSubcellularElementsPerCell = numSubCellularElementsPerCellAcross*+ numSubCellularElementsPerCellUp * numSubCellularElementsPerCellDeep;

    /**
     * Equlibrium distance between subcellular elements is given in Sandersius and Newman 2008.
     * doi:10.1088/1478-3975/5/1/015002.
     *
     * The packing density of spheres (circles) given is from Sandersius et al.
     * doi:10.1088/1478-3975/8/4/045007
     */

    double packing_density_of_spheres_2d = M_PI / (2.0 * sqrt(3.0));
    double packing_density_of_spheres_3d = M_PI / (3.0 * sqrt(2.0));


    double equlibrium_distance;
    switch (DIM)
    {
        case 2:
        {
            equlibrium_distance = 2.0 * sqrt(packing_density_of_spheres_2d / mTotalSubcellularElementsPerCell);
            break;
        }
        case 3:
        {
            equlibrium_distance = 2.0 * pow((packing_density_of_spheres_3d / mTotalSubcellularElementsPerCell), 1.0/3.0);
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    mNodeEquilibriumDistance = equlibrium_distance;

    // Create a cube of nodes for each cell.
    std::vector<PottsElement<DIM>* > elements;
    std::vector<Node<DIM>* > nodes;

    unsigned node_index_counter = 0;
    unsigned element_counter = 0;
    for (unsigned cells_deep = 0; cells_deep < numCellsDeep; cells_deep++)                                          // Cells in the z direction
    {
        for (unsigned cells_up = 0; cells_up < numCellsUp; cells_up++)                                                  // Cells in the y direction
        {
            for (unsigned cells_across = 0; cells_across < numCellsAcross; cells_across++)                                  // Cells in the x direction
            {
                std::vector<Node<DIM>* > local_nodes;
                for (unsigned nodes_deep = 0; nodes_deep < numSubCellularElementsPerCellDeep; nodes_deep++)                     // Nodes (subcellular elements) in the z direction
                {
                    for (unsigned nodes_up = 0; nodes_up < numSubCellularElementsPerCellUp; nodes_up++)                             // Nodes (subcellular elements) in the y direction
                    {
                        for (unsigned nodes_across = 0; nodes_across < numSubCellularElementsPerCellAcross; nodes_across++)             // Nodes (subcellular elements) in the x direction
                        {
                            Node<DIM>* p_node = new Node<DIM>(node_index_counter++,
                                                                false,
                                                                ((cells_across * numSubCellularElementsPerCellAcross + nodes_across) * equlibrium_distance),
                                                                ((cells_up * numSubCellularElementsPerCellUp + nodes_up) * equlibrium_distance),
                                                                ((cells_deep * numSubCellularElementsPerCellDeep + nodes_deep) * equlibrium_distance));
                            local_nodes.push_back(p_node);
                            nodes.push_back(p_node);
                        }
                    }
                }
                elements.push_back(new PottsElement<DIM>(element_counter++, local_nodes));
            }
        }
    }

    mpMesh = new SemMesh<DIM>(nodes, elements);
}

template<unsigned DIM>
SemMeshGenerator<DIM>::~SemMeshGenerator()
{
    delete mpMesh;
}

template<unsigned DIM>
double SemMeshGenerator<DIM>::GetEquilibriumDistance()
{
    return mNodeEquilibriumDistance;
}

template<unsigned DIM>
SemMesh<DIM>* SemMeshGenerator<DIM>::GetMesh()
{
    return mpMesh;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class SemMeshGenerator<1>;
template class SemMeshGenerator<2>;
template class SemMeshGenerator<3>;
