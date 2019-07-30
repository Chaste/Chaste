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


#ifndef AIRWAY_GENERATOR_HPP_
#define AIRWAY_GENERATOR_HPP_

#include "AirwayGeneration.hpp"

#include <deque>
#include <set>

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"

#if ((VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkSelectEnclosedPoints.h"

#include "vtkCellLocator.h"

/**
 * Airway Generator
 *
 * This class is an implementation of the airway tree growing algorithm described in
 * Tawhai et. al. 2004. J Appl Physiol. This class only handles growing airways
 * into a single  contiguous volume, such as a single lobe. For generation of a complete airway for a pair of lungs
 * see `MultiLobeAirwayGenerator`.
 */
class AirwayGenerator
{
    friend class TestAirwayGenerator;
    friend class TestMultiLobeAirwayGenerator;

public:

    /** The maximum number of generations the generator will create */
    static const unsigned MAX_GENERATIONS = 30;

    /**
     * Constructor
     *
     * @param lobeSurface A vtkPolyData containing a closed surface
     * @param minBranchLength minimum length a branch can have before it is considered terminal
     * @param pointLimit The minimum number of points an apex can contain
     * @param angleLimit The maximum branch angle, in degrees
     * @param branchingFraction The fraction a branch extends towards the centre of the point cloud
     */
    AirwayGenerator(vtkSmartPointer<vtkPolyData> lobeSurface,
                    double minBranchLength = 1.2,
                    unsigned pointLimit = 5,
                    double angleLimit = 90.0,
                    double branchingFraction = 0.5,
                    bool pointDistanceLimit = false);

    /** Destructor */
    ~AirwayGenerator();

    /**
     * Creates a cloud of regularly spaced points within the lobe surface given an
     * approximate number of required points
     *
     * @param rApproxPoints The number of (approximate) number of points to generate
     * @return A vtkPolyData containing the points
     */
    vtkSmartPointer<vtkPolyData> CreatePointCloudUsingTargetPoints(const unsigned& rApproxPoints);

    /**
     * Creates a cloud of regularly spaced points within the lobe surface given an
     * approximate volume for the space occupied around the points.
     *
     * @param rApproxVolume The approximate volume associated with the points in mm^3
     * @return A vtkPolyData containing the points
     */
    vtkSmartPointer<vtkPolyData> CreatePointCloudUsingTargetVolume(const double& rApproxVolume);

    /**
     * Creates a cloud of regularly spaced points within the lobe surface given an
     * approximate spacing between points
     *
     * @param rPointSpacing The spacing required between the points
     * @return A vtkPolyData containing the points
     */
    vtkSmartPointer<vtkPolyData> CreatePointCloud(const double& rPointSpacing);


    /**
     * Returns the global point cloud being used by the generator
     *
     * @return A vtkPolyData containing the points
     */
    vtkSmartPointer<vtkPolyData> GetPointCloud();

    /**
     * Returns a pointer to the Airway tree
     *
     * @return a pointer to the airway tree
     */
    vtkSmartPointer<vtkPolyData> GetAirwayTree();

    /**
     * Returns the center of mass of the active point cloud
     *
     * @param pointCloud A vtkPolyData containing the points to find the centre of mass of.
     * @param centre A preallocated array that the coordinates of the centre of mass will be written to
     */
    void GetCentreOfMass(vtkSmartPointer<vtkPolyData> pointCloud,
                         double centre[3]);

    /**
     * Splits a point cloud using a clipping plane
     *
     * @param pointCloud The set of points to be split
     * @param rNormal The normal of a plane to split the point cloud with
     * @param rOrigin A point on a plane to split the point cloud with
     * @param insideOut Should the cloud be split in the direction of the plane normal or the anti-normal
     *
     * @return The split point cloud
     */
    vtkSmartPointer<vtkPolyData> SplitPointCloud(vtkSmartPointer<vtkPolyData> pointCloud,
                                                 double rNormal[3],
                                                 double rOrigin[3],
                                                 bool insideOut = false);

    /**
     * Returns the angle between a potential new branch and an existing branch
     *
     * @param rCentre The point (centre of point cloud) towards which the apex wants to grow
     * @param rCurrentApex The apex representing a potential
     *
     * @return The angle of the branch
     */
    double CalculateBranchAngle(double rCentre[3], Apex& rCurrentApex);

    /**
     * Adds an initial growth apex to the system
     *
     * @param rStartLocation The spatial location where the apex will be located
     * @param rOriginalDirection The direction vector of the branch this apex is at the end of
     * @param rOriginalDirection The direction vector of the parent branch of this apex
     * @param rRadius The radius of the parent branch of this apex
     * @param rGeneration The generation of the branch this apex will produce
     */
    void AddInitialApex(double rStartLocation[3], double rOriginalDirection[3],  double rParentDirection[3], const double& rRadius, const unsigned& rGeneration);

    /**
     * Used internally to add a child growth apex to the system
     *
     * @param rStartId The ID of the initial point to grow the apex from
     * @param rOriginalDirection The direction vector of the parent branch of this apex
     * @param rParentDirection The direction vector of the parent branch of this apex
     * @param rGeneration The generation of the branch this apex will produce
     * @param pointCloud The point cloud this apex will try to grow into
     */
    void AddApex(const unsigned& rStartId, double rStartLocation[3], double rOriginalDirection[3], double rParentDirection[3], const unsigned& rGeneration);

    /**
     * Returns all the generation data structures
     *
     * @return All the generation data structures
     */
    std::deque<AirwayGeneration>& GetGenerations();

    /**
     * Returns a set of growth point IDs that have been invalidated.
     */
    std::set<unsigned>& GetInvalidIds();

    /**
     * Grows the given apex.
     *
     * @param rApex The apex to grow
     */
    void GrowApex(Apex& rApex);

    /**
     * Inserts a new branch into the airway tree
     *
     * @param pPointCloud The point cloud to grow towards
     * @param startId The id of the start point
     * @param originalDirection The direction of the parent of this branch
     * @param endLocation The location the branch is finally grown towards is written to this array
     * @return The ID of the inserted branch or -1 if insertion failed due to being outside the host volume
     */
    vtkIdType InsertBranch(vtkSmartPointer<vtkPolyData> pPointCloud,
                           unsigned startId,
                           double originalDirection[3],
                           double endLocation[3]);

    /**
     * Sets the closest generation seed point to a given location to be invalid
     *
     * @param point The location, does not have to be exactly equal to the seed point
     * @param searchCloud The cloud of points to find the nearest point in.
     */
    void InvalidateClosestPoint(double point[3], vtkSmartPointer<vtkPolyData> searchCloud);

    /**
     * Generates the entire tree
     */
    void Generate();

    /**
     * Checks the branch angle of a proposed growth centre and adjusts it within tolerance, if necessary
     *
     * @param startId The id of the branch start point
     * @param originalDirection The direction of the parent of the branch
     * @param centre The proposed point cloud centre that the apex will grow into (this will be updated with a corrected centre, if needed)
     */
    void CheckBranchAngleLengthAndAdjust(unsigned startId, double originalDirection[3], double centre[3]);

    /**
     * Calculates the Horsfield order of branches.
     */
    void CalculateHorsfieldOrder();

    /**
     * Calculates the radii of the branches
     *
     * Note, CalculateHorsfieldOrder must be called first
     *
     * @param rDiameterRatio The ratio in diameter between successive generations (see Horsfield 1978/1981)
     * @param rMaxDiameter The diameter of a branch of the highest order
     */
    void CalculateRadii(const double& rDiameterRatio);

    /**
     * Marks all points belonging to an initial apex
     *
     * A vtkDoubleArray named 'start_id' is associated to the airway tree. Any initial points are given a value of 1.0,
     * other points have the value 0.0
     */
    void MarkStartIds();

    /**
     * Calculates if a point is inside the lobe surface definition or not.
     *
     * @param point The point to test against the lobe surface
     * @return true if the point is inside the surface
     */
    bool IsInsideLobeSurface(double point[3]);

    /**
     * Calculates the distance of a point from the lobe surface definition.
     *
     * @param point The point to check the distance of
     * @return The distance of the point from the surface
     */
    double DistanceFromLobeSurface(double point[3]);

    /**
     * Computes the volume of the lobe.
     *
     * @return The volume of the lobe
     */
    double CalculateLobeVolume();

    /**
     * Writes the generated airway trees out as separate files, one per segment, for use in 3D-1D coupled ventilation models
     *
     * The files are written out as both vtkUnstructuredGrids and triangles/tetgen meshes.
     *
     * @param rOutputDirectory The directory to output the generated mesh to
     * @param rOutputFileNameRoot The root name of generated mesh files
     */
    void WriteDecomposedAirways(std::string rOutputDirectory, std::string rOutputFileNameRoot);

private:

    /** An enclosed surface representing the lobe geometry */
    vtkSmartPointer<vtkPolyData> mLobeSurface;

    /** A poly data representing the airway tree */
    vtkSmartPointer<vtkPolyData> mAirwayTree;

    /** The cloud of seed points */
    vtkSmartPointer<vtkPolyData> mSeedPointCloud;

    /** A point locator for easy access to the seed point cloud */
    vtkSmartPointer<vtkPointLocator> mSeedPointLocator;

    /** All the growth generations */
    std::deque<AirwayGeneration> mGenerations;

    /** A set of seed point Ids that are no longer valid for use as growth points */
    std::set<unsigned> mInvalidIds;

    /** The minimum length a branch can have before it is considered terminal */
    const double mLengthLimit;

    /** The minimum number of points a apex can have before it is considered terminal */
    const unsigned mPointLimit;

    /** The maximum branch angle, in degrees */
    const double mAngleLimit;

    /** The fraction a branch extends towards the centre of the point cloud */
    const double mBranchingFraction;

    /** The 'ratio' by which generated branch diameters decrease */
    double mDiameterRatio;

    /** A points selector to determine if we are inside the host volume */
    vtkSmartPointer<vtkSelectEnclosedPoints> mPointSelector;

    /** A list of initial node indices for use when calculating Horsfield order. */
    std::vector<unsigned> mStartIds;

    /** The initial radii corresponding to the initial node indices */
    std::vector<double> mStartRadii;

    /**
     * Private method to help recursively calculate the Horsfield order of the tree
     *
     * Order is calculated using a depth-first post-order traversal of the tree.
     *
     * @param pointId The point to be processed
     * @param pOrder Data structure to store orders in
     * @param rProcessedPoints VTK makes it difficult to visit the neighbours of a point, this datastructure is used to ensure the tree isn't traversed backwards
     *
     * @return The Horsfield order of this point
     */
    unsigned HorsfieldProcessPoint(unsigned pointId, vtkSmartPointer<vtkDoubleArray> pOrder, std::deque<unsigned>& rProcessedPoints);

    /**
     * Private method to help recursively calculate the radii of the tree
     *
     * @param pointId The point to be processed
     * @param startId The ID of the root of the tree corresponding to this point
     * @param rProcessedPoints VTK makes it difficult to visit the neighbours of a point, this datastructure is used to ensure the tree isn't traversed backwards
     */
    void RadiiProcessPoint(unsigned pointId, unsigned startId, std::deque<unsigned>& rProcessedPoints);
};

#endif //( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

#endif //CHASTE_VTK

#if !(defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6))
/**
 * This is a fake class to suppress coverage warnings. To get the real class
 * you must build with VTK of version 5.6 or above.
 */
class AirwayGenerator //This is here to suppress coverage warnings on machines that do not have vtk 5.6 or higher
{
public:
    /**
     * Fake constructor.
     */
    AirwayGenerator()
    {
        std::cout << "Dummy airway generator class for coverage" << std::endl;
    }
};
#endif //No VTK

#endif // AIRWAY_GENERATOR_HPP_
