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

#include "AirwayGeneration.hpp"

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <cassert>

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"
#include "vtkPolyVertex.h"

#if ((VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

AirwayGeneration::AirwayGeneration(unsigned generation_number) : mGenerationNumber(generation_number), mDistributionRadius(DBL_MAX)
{}

void AirwayGeneration::AddApex(unsigned startId, double currentLocation[3], double originalDirection[3], double parentBranch[3])
{
    Apex apex;
    apex.mStartId = startId;
    std::copy(currentLocation, currentLocation+3, apex.mCurrentLocation);
    std::copy(originalDirection, originalDirection+3, apex.mOriginalDirection);
    std::copy(parentBranch, parentBranch+3, apex.mParentDirection);

    apex.mPointCloud = vtkSmartPointer<vtkPolyData>::New();
    apex.mPointCloud->SetPoints(vtkSmartPointer<vtkPoints>::New());
    apex.mGeneration = mGenerationNumber;

    mApices.push_back(apex);
}

void AirwayGeneration::DistributeGrowthPoints(vtkSmartPointer<vtkPolyData> pAllGrowthPoints, std::set<unsigned>& invalidIds)
{
    if (mApices.size() == 0) // Nothing to do if we don't have any apices
    {
        return;
    }

    // Create a polydata containing the current apices
    vtkSmartPointer<vtkPolyData> apex_poly_data = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> apex_points = vtkSmartPointer<vtkPoints>::New();
    apex_poly_data->SetPoints(apex_points);
    apex_points->SetNumberOfPoints(mApices.size());

    for (std::deque<Apex>::iterator iter = mApices.begin();
        iter != mApices.end();
        ++iter)
    {
        apex_points->SetPoint(iter - mApices.begin(), iter->mCurrentLocation);
    }

    // Build a point locator to find apices
    vtkSmartPointer<vtkPointLocator> apex_locator = vtkSmartPointer<vtkPointLocator>::New();
    apex_locator->SetDataSet(apex_poly_data);
    apex_locator->BuildLocator();

    for (int index = 0; index < pAllGrowthPoints->GetNumberOfPoints(); ++index)
    {
        if (!invalidIds.count(index)) // Only copy over valid ids
        {
            // Find closest apex
            double growth_point[3];
            pAllGrowthPoints->GetPoints()->GetPoint(index, growth_point);

            double dist;
            int closest_apex_id = apex_locator->FindClosestPointWithinRadius(mDistributionRadius, growth_point, dist);

            if (closest_apex_id != -1) // No point within the radius
            {
                assert(sqrt(dist) <= mDistributionRadius);

                // Copy this point into the apex
                mApices[closest_apex_id].mPointCloud->GetPoints()->InsertNextPoint(growth_point);
            }
        }
    }

    // Some vtk filters get confused by individual points, so we add each point to a poly_vertex to keep them happy
    for (std::deque<Apex>::iterator apex_iter = GetApices().begin();
         apex_iter != GetApices().end();
         ++apex_iter)
    {
        vtkSmartPointer<vtkPolyVertex> poly_vertex = vtkSmartPointer<vtkPolyVertex>::New();
        poly_vertex->GetPointIds()->SetNumberOfIds(apex_iter->mPointCloud->GetNumberOfPoints());

        for (int i = 0; i < apex_iter->mPointCloud->GetNumberOfPoints(); ++i)
        {
            poly_vertex->GetPointIds()->SetId(i, i);
        }
        apex_iter->mPointCloud->Allocate(1,1);
        apex_iter->mPointCloud->InsertNextCell(poly_vertex->GetCellType(), poly_vertex->GetPointIds());
    }
}

std::deque<Apex>& AirwayGeneration::AirwayGeneration::GetApices()
{
    return mApices;
}

void AirwayGeneration::SetDistributionRadius(double radius)
{
    mDistributionRadius = radius;
}

double AirwayGeneration::GetDistributionRadius()
{
    return mDistributionRadius;
}

#endif //( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

#endif //CHASTE_VTK
