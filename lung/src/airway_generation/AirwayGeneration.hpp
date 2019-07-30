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


#ifndef AIRWAY_GENERATION_HPP_
#define AIRWAY_GENERATION_HPP_

#include <deque>
#include <set>
#include <iostream>

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"

#if ((VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

/*
 * This block is here to prevent a conflict when including vtkPointLocator.h:
 *     From VTK 7.0, there is a variable HZ in vtkPointLocator.h which is #def'd in linux header asm-generic/param.h,
 *     so we temporarily #undef it prior to the include.  See #2883 for details.
 */
#if (VTK_MAJOR_VERSION >= 7)
#pragma push_macro("HZ")
#undef HZ
#include <vtkPointLocator.h>
#pragma pop_macro("HZ")
#else // (VTK_MAJOR_VERSION < 7)
#include <vtkPointLocator.h>
#endif // (VTK_MAJOR_VERSION >= 7)

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkSelectEnclosedPoints.h"

/**
 * Auxiliary class representing a growth apex in the system
 *
 * A growth apex contains all the information required to grow a branch
 */
class Apex
{
public:
    /** The location of this growth apex */
    double mCurrentLocation[3];

    /** The direction of this growth apex */
    double mOriginalDirection[3];

    /** The direction of the parent branch of this growth apex */
    double mParentDirection[3];

    /** The ID of the the point in the airway tree to grow this apex from */
    vtkIdType mStartId;

    /** The generation number of this apex */
    unsigned mGeneration;

    /** The point cloud that this apex will attempt to grow into. */
    vtkSmartPointer<vtkPolyData> mPointCloud;
};

/**
 * Utility class representing a complete generation of growth apices
 */
class AirwayGeneration
{
public:
    /** Constructor */
    AirwayGeneration(unsigned generationNumber);

    /**
     * Add an apex to this generation
     *
     * @param startId The Id number of the start point of this growth apex
     * @param originalDirection The current direction of this growth apex
     * @param parentDirection The direction of the parent branch of this growth apex
     */
    void AddApex(unsigned startId, double currentLocation[3], double currentDirection[3], double parentDirection[3]);

    /**
     * Distributes current growth points to the growth apices
     *
     * @param pAllGrowthPoints Polydata containing the growth points
     * @param invalidIds A set of invalid growth point ids
     */
    void DistributeGrowthPoints(vtkSmartPointer<vtkPolyData> pAllGrowthPoints, std::set<unsigned>& invalidIds);

    /** Returns the apices associated with this generation */
    std::deque<Apex>& GetApices();

    /** Sets a maximum radius for points to be assigned to an apex */
    void SetDistributionRadius(double radius);

    /** Gets maximum radius for points to be assigned to an apex */
    double GetDistributionRadius();

private:
    /** The generation number of this generation object. */
    unsigned mGenerationNumber;

    /** All the growth apices associated with this generation. */
    std::deque<Apex> mApices;

    /** The maximum radius within which a point can be distributed to an apex */
    double mDistributionRadius;
};

#endif //( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

#endif //CHASTE_VTK

#if !(defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6))

/**
 * This is a fake class to suppress coverage warnings. To get the real class
 * you must build with VTK of version 5.6 or above.
 */
class AirwayGeneration // This is here to suppress coverage warnings on machines that do not have vtk 5.6 or higher
{
public:
    /**
     * Fake constructor.
     */
    AirwayGeneration()
    {
        std::cout << "Dummy airway generation class for coverage" << std::endl;
    }
};
#endif // No VTK

#endif // AIRWAY_GENERATION_HPP_
