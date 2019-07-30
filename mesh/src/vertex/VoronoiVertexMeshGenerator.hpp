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

#ifndef VORONOIVERTEXMESHGENERATOR_HPP_
#define VORONOIVERTEXMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "MutableVertexMesh.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "RandomNumberGenerator.hpp"

#include <boost/version.hpp>

#if BOOST_VERSION >= 105200
#include "boost/polygon/voronoi.hpp"

/**
 * Mesh generator that creates a 2D Voronoi tessellation using a number
 * of Lloyd's relaxation steps (http://en.wikipedia.org/wiki/Lloyd%27s_algorithm).
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class VoronoiVertexMeshGenerator
{
    friend class TestVoronoiVertexMeshGenerator;

protected:

    /** A pointer to the mesh that this class creates. */
    MutableVertexMesh<2,2>* mpMesh;

    /** A pointer to a toroidal vertex mesh (allowing for doubly periodic boundaries), whose creation is optional. */
    Toroidal2dVertexMesh* mpTorMesh;

    /** The number of elements requested across the mesh. */
    unsigned mNumElementsX;

    /** The number of elements requested up the mesh. */
    unsigned mNumElementsY;

    /** The total number of elements requested in the mesh. */
    unsigned mTotalNumElements;

    /** The number of Lloyd's relaxation steps requested in the Voronoi iteration. */
    unsigned mNumRelaxationSteps;

    /** The number of elements requested across or up the mesh, whichever is greater. */
    unsigned mMaxNumElems;

    /** The maximum integer position in the x-direction that a seed point could occupy after discretisation.  */
    int mMaxIntX;

    /** The maximum integer position in the y-direction that a seed point could occupy after discretisation. */
    int mMaxIntY;

    /** The requested average target area of elements in the mesh. */
    double mElementTargetArea;

    /** The scaling factor necessary to ensure that the x-component of seed locations lies between 0.0 and 1.0. */
    double mMultiplierInX;

    /** The scaling factor necessary to ensure that the y-component of seed locations lies between 0.0 and 1.0. */
    double mMultiplierInY;

    /** A floating point representation of MAX_INT/2 */
    double mSamplingMultiplier;

    /** The floating-point tolerance used for checking equality of node locations. */
    double mTol;

    /**
     * The maximum expected number of sides any polygon will have.
     * Based on experience, we initialise this to the default value 17 in the constructor.
     * It is theoretically possible that a polygon may be generated with more sides.
     */
    unsigned mMaxExpectedNumSidesPerPolygon;

    /**
     * Helper method for the constructor.
     *
     * Produces a vector of random seed points, with some validation, that lie in the rectangle
     * [0.0, mMultiplierInX] x [0.0, mMultiplierInY]. These seeds are used for the initial
     * Voronoi diagram construction, before any steps of Lloyd's relaxation.
     *
     * @return A vector of seed points
     */
    std::vector<c_vector<double, 2> > GetInitialPointLocations();

    /**
     * Helper method for the constructor.
     *
     * @return A vector of centroids corresponding to the elements currently held in mpMesh
     */
    std::vector<c_vector<double, 2> > GetElementCentroidsFromMesh();

    /**
     * Helper method for the constructor.
     *
     * Takes seed locations and updates mpMesh with Nodes and VertexElements corresponding to a
     * Voronoi diagram derived from the seed points. Boundary effects are eliminated by considering a 3x3 tessellation
     * of the seed points and only keeping those Voronoi cells corresponding to seed points in the centre tile.
     * Each seed point must lie in [0.0, mMultiplierInX] x [0.0, mMultiplierInY] for the tessellation to be valid.
     *
     * @param rSeedLocations A vector of seed locations for Voronoi cells
     */
    void CreateVoronoiTessellation(std::vector<c_vector<double, 2> >& rSeedLocations);

    /**
     * Helper method for the constructor.
     *
     * Validates the input parameters, and sets up remaining member variables.
     */
    void ValidateInputAndSetMembers();

    /**
     * Helper method for the constructor.
     *
     * Validates the initial random seed locations, making minute adjustments to seed locations if two points would be
     * mapped to the same integer lattice point after discretisation. This functionality is in a separate method to allow effective
     * unit testing.
     *
     * @param rSeedLocations The vector of random seed locations
     */
    void ValidateSeedLocations(std::vector<c_vector<double, 2> >& rSeedLocations);

    /**
     * Helper method for GetToroidalMesh().
     *
     * Takes the first boundary node, A, and checks it for congruence in position with all other boundary nodes. If a
     * congruent location is found, then A replaces the congruent node and this function returns true.  If none are
     * found, A is tagged as no longer being a boundary node.
     *
     * @param pMesh a pointer to the MutableVertexMesh we are operating on
     * @param width the periodic width of the mesh
     * @param height the periodic height of the mesh
     */
    bool CheckForCongruentNodes(MutableVertexMesh<2,2>* pMesh, double width, double height);

public:

    /**
     * Constructor.
     *
     * @param numElementsX  The number of elements requested across the mesh
     * @param numElementsY  The number of elements requested up the mesh
     * @param numRelaxationSteps  The number of Lloyd's Relaxation steps in the Voronoi iteration
     * @param elementTargetArea The requested average target area of elements in the mesh (defaults to 1.0)
     */
    VoronoiVertexMeshGenerator(unsigned numElementsX,
                               unsigned numElementsY,
                               unsigned numRelaxationSteps,
                               double elementTargetArea=1.0);

    /**
     * Null constructor for derived classes to call.
     */
    VoronoiVertexMeshGenerator()
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~VoronoiVertexMeshGenerator();

    /**
     * Helper function for the constructor: allows a new mesh to be generated by the user if they don't like the
     * current one.
     */
    void GenerateVoronoiMesh();

    /**
     * @return A pointer to a 2D mutable vertex mesh
     */
    virtual MutableVertexMesh<2,2>* GetMesh();

    /**
     * @return A pointer to a 2D mutable vertex mesh, after ReMesh() has been called to remove short edges
     */
    virtual MutableVertexMesh<2,2>* GetMeshAfterReMesh();

    /**
     * @return A pointer to a 2D toroidal vertex mesh with periodic boundaries
     */
    virtual Toroidal2dVertexMesh* GetToroidalMesh();

    /**
     * @return A vector representing the polygon distribution of the generated mesh, triangles upwards
     */
    std::vector<double> GetPolygonDistribution();

    /**
     * Computes the coefficient of variation of the areas of elements in the mesh, defined to be the sample standard
     * deviation in area divided by the mean area.
     *
     * @return The coefficient of variation of the area of elements in the mesh
     */
    double GetAreaCoefficientOfVariation();

    /**
     * Call GenerateVoronoiMesh().
     *
     * Allows the user to make a new mesh, for instance if trying to generate a specific polygon distribution.
     */
    void RefreshSeedsAndRegenerateMesh();

    /**
     * Set mMaxNumSidesPerPolygon.
     *
     * @param maxExpectedNumSidesPerPolygon the new value of mMaxNumSidesPerPolygon
     */
    void SetMaxExpectedNumSidesPerPolygon(unsigned maxExpectedNumSidesPerPolygon);

    /**
     * @return mMaxNumSidesPerPolygon
     */
    unsigned GetMaxExpectedNumSidesPerPolygon();
};

#endif // BOOST_VERSION >= 105200

#if BOOST_VERSION < 105200

/**
 * This is a fake class to suppress coverage warnings. To get the real class
 * you must build with Boost version 1.52 or above.
 */
class VoronoiVertexMeshGenerator
{
public:

    /**
     * Fake constructor.
     */
    VoronoiVertexMeshGenerator();
};

#endif // BOOST_VERSION < 105200

#endif /*VORONOIVERTEXMESHGENERATOR_HPP_*/
