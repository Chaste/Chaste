/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "VoronoiPrism3dVertexMeshGenerator.hpp"

#include "MonolayerVertexMeshGenerator.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#if BOOST_VERSION >= 105200

VoronoiPrism3dVertexMeshGenerator::VoronoiPrism3dVertexMeshGenerator(unsigned numElementsX,
                                                                     unsigned numElementsY,
                                                                     double elementHeightZ,
                                                                     unsigned numRelaxationSteps,
                                                                     double elementTargetApicalArea)
        : mpMesh(NULL),
          mNumElementsX(numElementsX),
          mNumElementsY(numElementsY),
          mElementHeightZ(elementHeightZ),
          mNumRelaxationSteps(numRelaxationSteps),
          mElementTargetApicalArea(elementTargetApicalArea),
          mMaxExpectedNumSidesPerPolygon(17)
{
    if ((mNumElementsX < 2) || (mNumElementsY < 2))
    {
        EXCEPTION("Need at least 2 by 2 cells");
    }
    if (mElementHeightZ <= 0.0)
    {
        EXCEPTION("Specified element height must be strictly positive");
    }
    if (mElementTargetApicalArea <= 0.0)
    {
        EXCEPTION("Specified target apical area must be strictly positive");
    }

    this->GenerateVoronoiMesh();
}

VoronoiPrism3dVertexMeshGenerator::~VoronoiPrism3dVertexMeshGenerator()
{
    if (mpMesh)
    {
        delete mpMesh;
    }
}

void VoronoiPrism3dVertexMeshGenerator::GenerateVoronoiMesh()
{
    VoronoiVertexMeshGenerator generator2(mNumElementsX, mNumElementsY, mNumRelaxationSteps,
                                          mElementTargetApicalArea);
    MutableVertexMesh<2, 2>* p_mesh2 = generator2.GetMesh();
    MonolayerVertexMeshGenerator generator("", false);
    mpMesh = generator.MakeMeshUsing2dMesh(*p_mesh2, mElementHeightZ);
}

MutableVertexMesh<3, 3>* VoronoiPrism3dVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}

MutableVertexMesh<3, 3>* VoronoiPrism3dVertexMeshGenerator::GetMeshAfterReMesh()
{
    mpMesh->ReMesh();
    return mpMesh;
}

std::vector<double> VoronoiPrism3dVertexMeshGenerator::GetPolygonDistribution()
{
    assert(mpMesh != NULL);

    // Number of elements in the mesh
    unsigned num_elems = mpMesh->GetNumElements();

    // Store the number of each class of polygons
    std::vector<unsigned> num_polygons(mMaxExpectedNumSidesPerPolygon - 2, 0);

    // Container to return the polygon distribution
    std::vector<double> polygon_dist;

    // Loop over elements in the mesh to get the number of each class of polygon
    for (unsigned elem_idx = 0; elem_idx < num_elems; ++elem_idx)
    {
        unsigned num_nodes_this_elem = mpMesh->GetElement(elem_idx)->GetNumNodes();

        // All polygons are assumed to have 3, 4, 5, ..., mMaxExpectedNumSidesPerPolygon sides
        // and since there should be an upper node for every lower node, it should be even number
        assert(num_nodes_this_elem % 2 == 0); ///\todo check if the pairs is really lower and upper after simulation #2850
        assert(num_nodes_this_elem > 2);
        assert(num_nodes_this_elem / 2 <= mMaxExpectedNumSidesPerPolygon);

        // Increment correct place in counter - triangles in place 0, squares in 1, etc
        num_polygons[num_nodes_this_elem / 2 - 3]++;
    }

    // Loop over the vector of polygon numbers and calculate the distribution vector to return
    unsigned elems_accounted_for = 0;
    for (unsigned polygon = 0; polygon < num_polygons.size(); polygon++)
    {
        elems_accounted_for += num_polygons[polygon];

        polygon_dist.push_back(static_cast<double>(num_polygons[polygon]) / static_cast<double>(num_elems));

        // Only fill the vector of polygon distributions to the point where there are none of higher class in the mesh
        if (elems_accounted_for == num_elems)
        {
            break;
        }
    }

    return polygon_dist;
}

double VoronoiPrism3dVertexMeshGenerator::GetApicalAreaCoefficientOfVariation()
{
    assert(mpMesh != NULL);

    // Number of elements in the mesh, and check there are at least two
    unsigned num_elems = mpMesh->GetNumElements();
    assert(num_elems > 1);

    double var = 0.0;

    // Loop over elements in the mesh to get the contributions to the variance
    for (unsigned elem_idx = 0; elem_idx < num_elems; elem_idx++)
    {
        double deviation = mpMesh->GetVolumeOfElement(elem_idx) / mElementHeightZ - mElementTargetApicalArea;
        var += deviation * deviation;
    }

    var /= static_cast<double>(num_elems - 1);

    return sqrt(var) / mElementTargetApicalArea;
}

void VoronoiPrism3dVertexMeshGenerator::SetMaxExpectedNumSidesPerPolygon(unsigned maxExpectedNumSidesPerPolygon)
{
    mMaxExpectedNumSidesPerPolygon = maxExpectedNumSidesPerPolygon;
}

unsigned VoronoiPrism3dVertexMeshGenerator::GetMaxExpectedNumSidesPerPolygon()
{
    return mMaxExpectedNumSidesPerPolygon;
}

#endif // BOOST_VERSION >= 105200
#if BOOST_VERSION < 105200

/**
 * This is a fake class to suppress coverage warnings. To get the real class
 * you must build with Boost version 1.52 or above.
 */
VoronoiPrism3dVertexMeshGenerator::VoronoiPrism3dVertexMeshGenerator()
{
    EXCEPTION("This is a dummy class. Build with Boost version 1.52 or above for functionality.");
}

#endif // BOOST_VERSION < 105200
