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


#ifndef MULTI_LOBE_AIRWAY_GENERATOR_HPP_
#define MULTI_LOBE_AIRWAY_GENERATOR_HPP_

#include <utility>
#include <iostream>

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"

#include "TetrahedralMesh.hpp"
#include "AirwayGenerator.hpp"
#include "LungTools.hpp"

#if  ((VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)


/**
 * Multi Lobe Airway Generator
 *
 * This class wraps `AirwayGenerator` to facilitate easy generation of airway trees for a set of lungs with
 * multiple lobes. The class assumes that the major airways have been obtained from images and are available
 * to provide start points for airways generation.
 */
class MultiLobeAirwayGenerator
{
    friend class TestMultiLobeAirwayGenerator;

public:
    /**
     * Constructor
     *
     * @param rAirwaysMesh Mesh of the major airways. The generated airways will be appended to this mesh
     * @param pointDistanceLimit Specifies whether to limit the size of terminal generations using a heuristic limit
     */
    MultiLobeAirwayGenerator(TetrahedralMesh<1,3>& rAirwaysMesh, bool pointDistanceLimit = false);

    /**
     * Destructor
     */
    ~MultiLobeAirwayGenerator();

    /**
     * Allows the user to add a lung lobe from a file
     *
     * This method assumes that the file is a .stl file.
     *
     * @param rFileName The name of the STL file representing the lobe surface
     * @param lungLocation Is this a right or a left lung lobe?
     */
    void AddLobe(const std::string& rFileName, LungLocation lungLocation);

    /**
     * Adds a lobe to the mesh generator
     *
     * @param pLobeSurface A VTK poly data containing the surface of the lobe
     * @param lungLocation Is this a right or a left lung lobe?
     */
    void AddLobe(vtkSmartPointer<vtkPolyData> pLobeSurface, LungLocation lungLocation);

    /**
     * Counts the number of lobes added in each lung location
     *
     * @param lungLocation Is this a right or a left lung lobe?
     * @return The number of lobes in this lung location
     */
    unsigned GetNumLobes(LungLocation lungLocation);

    /**
     * Assigns the end points of the major airway tree to initial growth apices in the correct lobes
     */
    void AssignGrowthApices();

    /**
     * Distributes the correct number of points to the clouds of the individual lobe generators
     */
    void DistributePoints();

    /**
     * Sets the target number of generation points for each lung.
     *
     * Note this is PER LUNG, not PER PAIR OF LUNGS.
     *
     * Either this method should be used or SetPointVolume should be used, not both
     *
     * @param rNumberOfPointsPerLung The target number of generation points for each lung
     */
    void SetNumberOfPointsPerLung(const unsigned& rNumberOfPointsPerLung );

    /**
     * Sets the target volume of generation points.
     *
     * Either this method should be used or SetNumberOfPointsPerLung should be used, not both
     *
     * @param rVolume The target volume associated with each generation point
     */
    void SetPointVolume(const double& rVolume );

    /**
     * Sets the minimum branch length for generated airways.
     *
     * @param rMinimumbranchLength The minimum branch length for generated airways
     */
    void SetMinimumBranchLength(const double& rMinimumbranchLength);

    /**
     * Sets the minimum number of points in a point cloud
     *
     * @param rPointLimit The minimum number of points in a point cloud
     */
    void SetPointLimit(const unsigned& rPointLimit);

    /**
     * Sets the maximum angle with the parent branch for generated branches
     *
     * @param rAngleLimit The maximum angle with the parent branch for generated branches
     */
    void SetAngleLimit(const double& rAngleLimit);

    /**
     * Sets the fraction of the distance towards the centre of the point cloud assigned to a new branch's length
     *
     * @param rBranchingFraction The branch fraction
     */
    void SetBranchingFraction(const double& rBranchingFraction);

    /**
     * The scaling factor for generated airway diameters
     *
     * @param rDiameterRatio Scaling factor for generated airway diameters
     */
    void SetDiameterRatio(const double& rDiameterRatio);

    /**
     * Generates a complete airway tree and writes it to the given file name.
     *
     * @param rOutputDirectory The directory to output the generated mesh to
     * @param rOutputFileNameRoot The root name of generated mesh files
     */
    void Generate(std::string rOutputDirectory, std::string rOutputFileNameRoot);

private:

    /** A mesh containing the major airways.  */
    TetrahedralMesh<1,3>& mAirwaysMesh;

    /** A vector containing the left lobes associated with this generator */
    std::vector<std::pair<AirwayGenerator*, LungLocation> > mLobeGenerators;

    /** The target number of points per lung */
    unsigned mNumberOfPointsPerLung;

    /** The target Volume for points in lung */
    double mPointVolume;

    /** The minimum branch length for generated airways */
    double mMinimumBranchLength;

    /** The minimum number of points in a point cloud */
    unsigned mPointLimit;

    /** The maximum angle with the parent branch for generated branches */
    double mAngleLimit;

    /** The fraction of the distance towards the centre of the point cloud assigned to a new branch's length */
    double mBranchingFraction;

    /** The scaling factor for generated airway diameters */
    double mDiameterRatio;

    /** A flag to turn on the point distance limit heuristic */
    bool mPointDistanceLimit;
};

#endif // (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6

#endif //CHASTE_VTK

#if !(defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6))

/**
 * This is a fake class to suppress coverage warnings. To get the real class
 * you must build with VTK of version 5.6 or above.
 */
class MultiLobeAirwayGenerator //This is here to suppress coverage warnings on machines that do not have vtk 5.6 or higher
{
public:
    /**
     * Fake constructor.
     */
    MultiLobeAirwayGenerator()
    {
        std::cout << "Dummy multi lobe airway generator class for coverage" << std::endl;
    }
};
#endif //No VTK

#endif // MULTI_LOBE_AIRWAY_GENERATOR_HPP_
