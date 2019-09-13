
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

#include "HoneycombMeshGenerator.hpp"

#include <boost/foreach.hpp>
#include "TrianglesMeshReader.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "MathsCustomFunctions.hpp"
#include "ChasteSyscalls.hpp"

HoneycombMeshGenerator::HoneycombMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts, double scaleFactor)
  : mpMesh(nullptr),
    mMeshFilename("mesh"),
    mDomainWidth(numNodesAlongWidth*scaleFactor),
    mNumCellWidth(numNodesAlongWidth), //*1 because cells are considered to be size one
    mNumCellLength(numNodesAlongLength)
{
    // The code below won't work in parallel
    assert(PetscTools::IsSequential());

    // An older version of the constructor might allow the wrong argument through to the scale factor
    assert(scaleFactor > 0.0);

    // Get a unique temporary foldername
    std::stringstream pid;
    pid << getpid();
    OutputFileHandler output_file_handler("2D_temporary_honeycomb_mesh_" + pid.str());
    std::string output_dir = output_file_handler.GetOutputDirectoryFullPath();

    unsigned num_nodes_along_width = mNumCellWidth;
    unsigned num_nodes_along_length = mNumCellLength;
    double horizontal_spacing = mDomainWidth / (double)num_nodes_along_width;
    double vertical_spacing = (sqrt(3.0)/2)*horizontal_spacing;

    // This line is needed to define ghost nodes later
    mDomainDepth = (double)(num_nodes_along_length) * vertical_spacing;

    // Take account of ghost nodes
    num_nodes_along_width +=  2*ghosts;
    num_nodes_along_length += 2*ghosts;

    unsigned num_nodes            = num_nodes_along_width*num_nodes_along_length;
    unsigned num_elem_along_width = num_nodes_along_width-1;
    unsigned num_elem_along_length = num_nodes_along_length-1;
    unsigned num_elem             = 2*num_elem_along_width*num_elem_along_length;
    unsigned num_edges            = 3*num_elem_along_width*num_elem_along_length + num_elem_along_width + num_elem_along_length;

    double x0 = -horizontal_spacing*ghosts;
    double y0 = -vertical_spacing*ghosts;

    mBottom = -vertical_spacing*ghosts;
    mTop = mBottom + vertical_spacing*(num_nodes_along_length-1);

    // Write node file
    out_stream p_node_file = output_file_handler.OpenOutputFile(mMeshFilename+".node");
    (*p_node_file) << std::scientific;
    (*p_node_file) << std::setprecision(20);
    (*p_node_file) << num_nodes << "\t2\t0\t1" << std::endl;

    unsigned node = 0;
    for (unsigned i=0; i<num_nodes_along_length; i++)
    {
        for (unsigned j=0; j<num_nodes_along_width; j++)
        {
            if (i<ghosts || i>=(ghosts+mNumCellLength))
            {
                mGhostNodeIndices.insert(node);
            }
            else if (j < ghosts || j >= (ghosts+mNumCellWidth))
            {
                mGhostNodeIndices.insert(node);
            }
            unsigned boundary = 0;
            if (i==0 || i==num_nodes_along_length-1 || j==0 || j==num_nodes_along_width-1)
            {
                boundary = 1;
            }

            double x = x0 + horizontal_spacing*((double)j + 0.25*(1.0+ SmallPow(-1.0,i+1)));
            double y = y0 + vertical_spacing*(double)i;

            // Avoid floating point errors which upset OffLatticeSimulation
            if ((y<0.0) && (y>-1e-12))
            {
                // Difficult to cover - just corrects floating point errors that have occurred from time to time!
                // LCOV_EXCL_START
                y = 0.0;
                // LCOV_EXCL_STOP
            }

            (*p_node_file) << node++ << "\t" << x << "\t" << y << "\t" << boundary << std::endl;
        }
    }
    p_node_file->close();

    // Write element file and edge file
    out_stream p_elem_file = output_file_handler.OpenOutputFile(mMeshFilename+".ele");
    (*p_elem_file) << std::scientific;

    out_stream p_edge_file = output_file_handler.OpenOutputFile(mMeshFilename+".edge");
    (*p_node_file) << std::scientific;

    (*p_elem_file) << num_elem << "\t3\t0" << std::endl;
    (*p_edge_file) << num_edges << "\t1" << std::endl;

    unsigned elem = 0;
    unsigned edge = 0;
    for (unsigned i=0; i<num_elem_along_length; i++)
    {
        for (unsigned j=0; j < num_elem_along_width; j++)
        {
            unsigned node0 =     i*num_nodes_along_width + j;
            unsigned node1 =     i*num_nodes_along_width + j+1;
            unsigned node2 = (i+1)*num_nodes_along_width + j;

            if (i%2 != 0)
            {
                node2 = node2 + 1;
            }

            (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;

            unsigned horizontal_edge_is_boundary_edge = 0;
            unsigned vertical_edge_is_boundary_edge = 0;
            if (i==0)
            {
                horizontal_edge_is_boundary_edge = 1;
            }
            if (j==0 && i%2==0)
            {
                vertical_edge_is_boundary_edge = 1;
            }

            (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << horizontal_edge_is_boundary_edge << std::endl;
            (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node2 << "\t" << 0 << std::endl;
            (*p_edge_file) << edge++ << "\t" << node2 << "\t" << node0 << "\t" << vertical_edge_is_boundary_edge << std::endl;

            node0 = i*num_nodes_along_width + j + 1;

            if (i%2 != 0)
            {
                node0 = node0 - 1;
            }
            node1 = (i+1)*num_nodes_along_width + j+1;
            node2 = (i+1)*num_nodes_along_width + j;

            (*p_elem_file) << elem++ << "\t" << node0 << "\t" << node1 << "\t" << node2 << std::endl;
        }
    }

    for (unsigned i=0; i<num_elem_along_length; i++)
    {
        unsigned node0, node1;

        if (i%2==0)
        {
             node0 = (i+1)*num_nodes_along_width - 1;
             node1 = (i+2)*num_nodes_along_width - 1;
        }
        else
        {
            node0 = (i+1)*num_nodes_along_width;
            node1 = (i)*num_nodes_along_width;
        }
        (*p_edge_file) << edge++ << "\t" << node0 << "\t" << node1 << "\t" << 1 << std::endl;
    }

    for (unsigned j=0; j<num_elem_along_width; j++)
    {
        unsigned node0 = num_nodes_along_width*(num_nodes_along_length-1) + j;
        unsigned node1 = num_nodes_along_width*(num_nodes_along_length-1) + j+1;
        (*p_edge_file) << edge++ << "\t" << node1 << "\t" << node0 << "\t" << 1 << std::endl;
    }

    p_elem_file->close();
    p_edge_file->close();

    // Having written the mesh to file, now construct it using TrianglesMeshReader.
    // Nested scope so the reader closes files before we try deleting them below the scope.
    {
        TrianglesMeshReader<2,2> mesh_reader(output_dir + mMeshFilename);
        mpMesh = new MutableMesh<2,2>;
        mpMesh->ConstructFromMeshReader(mesh_reader);
    }

    // Delete the temporary folder
    output_file_handler.FindFile("").Remove();

    // The original files have been deleted, it is better if the mesh object forgets about them
    mpMesh->SetMeshHasChangedSinceLoading();
}

HoneycombMeshGenerator::~HoneycombMeshGenerator()
{
    delete mpMesh;
}

MutableMesh<2,2>* HoneycombMeshGenerator::GetMesh()
{
    return mpMesh;
}

std::vector<unsigned> HoneycombMeshGenerator::GetCellLocationIndices()
{
    std::vector<unsigned> location_indices;

    for (unsigned i=0; i<mpMesh->GetNumNodes(); i++)
    {
        if (mGhostNodeIndices.find(i)==mGhostNodeIndices.end())
        {
            location_indices.push_back(i);
        }
    }
    return location_indices;
}

MutableMesh<2,2>* HoneycombMeshGenerator::GetCircularMesh(double radius)
{
    if (!mGhostNodeIndices.empty())
    {
        EXCEPTION("Cannot call GetCircularMesh on a HoneycombMesh with ghost nodes");
    }

    // Centre the mesh at (0,0)
    c_vector<double,2> centre = zero_vector<double>(2);
    for (unsigned i=0; i<mpMesh->GetNumNodes(); i++)
    {
        centre += mpMesh->GetNode(i)->rGetLocation();
    }
    centre /= (double)mpMesh->GetNumNodes();

    mpMesh->Translate(-centre[0], -centre[1]);

    // Iterate over nodes, deleting any that lie more than the specified radius from (0,0)
    for (unsigned i=0; i<mpMesh->GetNumAllNodes(); i++)
    {
        if (norm_2(mpMesh->GetNode(i)->rGetLocation()) >= radius)
        {
            mpMesh->DeleteNodePriorToReMesh(i);
        }
        else
        {
            // Jiggle the data
            c_vector<double,2>& r_location = mpMesh->GetNode(i)->rGetModifiableLocation();
            c_vector<double,2> shift;
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            double max_jiggle = radius*5e-6;
            shift[0] = max_jiggle*(p_gen->ranf()-0.5);
            shift[1] = max_jiggle*(p_gen->ranf()-0.5);
            r_location += shift;
        }
    }

    // Remesh
    NodeMap map(mpMesh->GetNumNodes());
    mpMesh->ReMesh(map);

    return mpMesh;
}

double HoneycombMeshGenerator::GetDomainDepth()
{
    return mDomainDepth;
}

double HoneycombMeshGenerator::GetDomainWidth()
{
    return mDomainWidth;
}
