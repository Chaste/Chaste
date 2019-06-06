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

#include "AirwayGenerator.hpp"

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <sstream>

#include "VtkMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "TetrahedralMesh.hpp"
#include "OutputFileHandler.hpp"
#include <UblasIncludes.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"

#if ((VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

#include "vtkGeometryFilter.h"
#include "vtkDoubleArray.h"
#include "vtkPlane.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyVertex.h"
#include "vtkMassProperties.h"
#include "vtkMath.h"
#include "vtkTableBasedClipDataSet.h"
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMassProperties.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkAppendFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkSelectEnclosedPoints.h"

AirwayGenerator::AirwayGenerator(vtkSmartPointer<vtkPolyData> LobeSurface,
                                 double branchLengthLimit,
                                 unsigned pointLimit,
                                 double angleLimit,
                                 double branchingFraction,
                                 bool pointDistanceLimit) : mLobeSurface(LobeSurface),
                                                                  mAirwayTree(vtkSmartPointer<vtkPolyData>::New()),
                                                                  mLengthLimit(branchLengthLimit),
                                                                  mPointLimit(pointLimit),
                                                                  mAngleLimit(angleLimit),
                                                                  mBranchingFraction(branchingFraction),
                                                                  mPointSelector(vtkSmartPointer<vtkSelectEnclosedPoints>::New())

{
    mAirwayTree->SetPoints(vtkSmartPointer<vtkPoints>::New());
    mAirwayTree->SetLines(vtkSmartPointer<vtkCellArray>::New());

    // Set up a filter to check if a point is inside the surface or not
    mPointSelector->CheckSurfaceOn();
    mPointSelector->Initialize(mLobeSurface);
    mPointSelector->SetTolerance(1e-6);

    // Calculate the bounding box of the lobe to allow progressive point reassignment limit to be used
    c_vector<double, 6> lobe_bounds;
    mLobeSurface->GetBounds(lobe_bounds.data());
    double lobe_bound_size = norm_2(subslice(lobe_bounds, 0,2,3) - subslice(lobe_bounds, 1,2,3));

    for (unsigned generation_number = 0; generation_number < MAX_GENERATIONS; ++generation_number)
    {
        AirwayGeneration gen(generation_number);

        if (pointDistanceLimit)
        {
            // Heuristic formula to determine distribution radius (see #2275)
            double scale_distance_limit = lobe_bound_size/30; //30 is the expected total number of generations
            gen.SetDistributionRadius(std::max(lobe_bound_size - scale_distance_limit*generation_number, 5.0));
        }
        mGenerations.push_back(gen);
    }
}

AirwayGenerator::~AirwayGenerator()
{
    mPointSelector->Complete();
}

vtkSmartPointer<vtkPolyData> AirwayGenerator::CreatePointCloudUsingTargetPoints(const unsigned& rApproxPoints)
{
    //Determine the spacing of the points being generated
    vtkSmartPointer<vtkMassProperties> mass_properties = vtkSmartPointer<vtkMassProperties>::New();
#if VTK_MAJOR_VERSION >= 6
    mass_properties->SetInputData(mLobeSurface);
#else
    mass_properties->SetInput(mLobeSurface);
#endif

    double point_spacing = std::pow(mass_properties->GetVolume()/rApproxPoints, 1.0/3.0);

    return CreatePointCloud(point_spacing);
}

vtkSmartPointer<vtkPolyData> AirwayGenerator::CreatePointCloudUsingTargetVolume(const double& rApproxDensity)
{
    double point_spacing = std::pow(rApproxDensity, 1.0/3.0);
    return CreatePointCloud(point_spacing);
}

vtkSmartPointer<vtkPolyData> AirwayGenerator::CreatePointCloud(const double& rPointSpacing)
{
    double bounds[6];
    mLobeSurface->GetBounds(bounds);

    unsigned xi_max = std::ceil((bounds[1] - bounds[0])/rPointSpacing);
    unsigned yi_max = std::ceil((bounds[3] - bounds[2])/rPointSpacing);
    unsigned zi_max = std::ceil((bounds[5] - bounds[4])/rPointSpacing);

    // Generate the point cloud
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    double node_index = 0;

    for (unsigned xi = 0; xi < xi_max; ++xi)
    {
        for (unsigned yi = 0; yi < yi_max; ++yi)
        {
            for (unsigned zi = 0; zi < zi_max; ++zi)
            {
                double x = bounds[0] + xi*rPointSpacing;
                double y = bounds[2] + yi*rPointSpacing;
                double z = bounds[4] + zi*rPointSpacing;
                // Repeat test for "isinside" because, prior to VTK 8.2, it may fail when the point is on surface (#3002)
                if (mPointSelector->IsInsideSurface(x, y, z) || mPointSelector->IsInsideSurface(x, y, z) || mPointSelector->IsInsideSurface(x, y, z))
                {
                    points->InsertPoint(node_index, x, y, z);
                    node_index += 1;
                }
            }
        }
    }

    // Some vtk filters get confused by individual points, so we add each point to a poly_vertex to keep them happy
    vtkSmartPointer<vtkPolyVertex> poly_vertex = vtkSmartPointer<vtkPolyVertex>::New();
    poly_vertex->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());

    for (int i = 0; i < points->GetNumberOfPoints(); ++i)
    {
        poly_vertex->GetPointIds()->SetId(i, i);
    }

    mSeedPointCloud = vtkSmartPointer<vtkPolyData>::New();
    mSeedPointCloud->SetPoints(points);
    mSeedPointCloud->Allocate(1,1);
    mSeedPointCloud->InsertNextCell(poly_vertex->GetCellType(), poly_vertex->GetPointIds());

    mSeedPointLocator = vtkSmartPointer<vtkPointLocator>::New();
    mSeedPointLocator->SetDataSet(mSeedPointCloud);
    mSeedPointLocator->BuildLocator();

    return mSeedPointCloud;
}

vtkSmartPointer<vtkPolyData> AirwayGenerator::GetPointCloud()
{
    return mSeedPointCloud;
}

vtkSmartPointer<vtkPolyData> AirwayGenerator::GetAirwayTree()
{
    return mAirwayTree;
}

void AirwayGenerator::GetCentreOfMass(vtkSmartPointer<vtkPolyData> pointCloud,
                                      double centre[3])
{
    assert(pointCloud->GetNumberOfPoints() > 0);

    double point[3];

    centre[0] = 0.0;
    centre[1] = 0.0;
    centre[2] = 0.0;

    for (int i = 0; i < pointCloud->GetNumberOfPoints(); ++i)
    {
        pointCloud->GetPoint(i, point);
        for (unsigned j = 0; j < 3; ++j)
        {
            centre[j] += point[j];
        }
    }

    for (unsigned j = 0; j < 3; ++j)
    {
        centre[j] /= pointCloud->GetNumberOfPoints();
    }
}

vtkSmartPointer<vtkPolyData> AirwayGenerator::SplitPointCloud(vtkSmartPointer<vtkPolyData> pointCloud,
                                                              double rNormal[3],
                                                              double rOrigin[3],
                                                              bool insideOut)
{
    vtkSmartPointer<vtkPlane> clip_plane = vtkSmartPointer<vtkPlane>::New();
    clip_plane->SetNormal(rNormal);
    clip_plane->SetOrigin(rOrigin);

    vtkSmartPointer<vtkTableBasedClipDataSet > clipper = vtkSmartPointer<vtkTableBasedClipDataSet >::New();
#if VTK_MAJOR_VERSION >= 6
    clipper->SetInputData(pointCloud);
#else
    clipper->SetInputConnection(pointCloud->GetProducerPort());
#endif
    clipper->SetClipFunction(clip_plane);
    clipper->GenerateClippedOutputOff();
    clipper->SetInsideOut(insideOut);
    clipper->Update();

    vtkSmartPointer<vtkGeometryFilter> filter = vtkSmartPointer<vtkGeometryFilter>::New();
    filter->SetInputConnection(clipper->GetOutputPort());
    filter->Update();

    return filter->GetOutput();
}

void AirwayGenerator::AddInitialApex(double rStartLocation[3],
                                     double rOriginalDirection[3],
                                     double rParentDirection[3],
                                     const double& rRadius,
                                     const unsigned& rGeneration)
{
    // Create the vtk point and add to the airway tree
    const unsigned start_id = mAirwayTree->GetPoints()->InsertNextPoint(rStartLocation);
    AddApex(start_id, rStartLocation, rOriginalDirection, rParentDirection, rGeneration);

    // The starting point id is stored to allow Horsfield order to be calculated correctly
    mStartIds.push_back(start_id);
    mStartRadii.push_back(rRadius); // Radius must be retained for resulting radii calculations
}

void AirwayGenerator::AddApex(const unsigned& rStartId,
                              double rStartLocation[3],
                              double rOriginalDirection[3],
                              double rParentDirection[3],
                              const unsigned& rGeneration)
{
    if (rGeneration < MAX_GENERATIONS)
    {
        mGenerations[rGeneration].AddApex(rStartId, rStartLocation, rOriginalDirection, rParentDirection);
    }
    else
    {
        std::stringstream message;
        message << "Error: Airway generation can only generate up to ";
        message << MAX_GENERATIONS;
        message << " generations.";
        EXCEPTION(message.str());
    }
}

std::deque<AirwayGeneration>& AirwayGenerator::GetGenerations()
{
    return mGenerations;
}

std::set<unsigned>& AirwayGenerator::GetInvalidIds()
{
    return mInvalidIds;
}

void AirwayGenerator::GrowApex(Apex& rApex)
{
    if (rApex.mPointCloud->GetNumberOfPoints() == 0)
    {
        return; //Can't grow an apex without any points associated to it
    }

    // Determine the current point cloud centre of mass
    double centre[3];
    GetCentreOfMass(rApex.mPointCloud, centre);

    // Determine the splitting plane: normal to the direction of the existing branch and the direction to the centre of the point cloud
    double centre_direction[3];
    vtkMath::Subtract(centre, rApex.mCurrentLocation, centre_direction);
    vtkMath::Normalize(centre_direction);

    double normal[3];
    vtkMath::Cross(rApex.mOriginalDirection, centre_direction, normal);

    // If the current direction and the vector to the centre are colinear then the normal is not well defined.
    // Tawhai 2004 does not define what should happen in this case. We fall back to the approach in Tawhai 2000:
    // Use the plane made by the branch and its parent as the splitting plane.
    if (vtkMath::Norm(normal) < 1e-10) //This tolerance is arbitrary, need to think about what is appropriate here.
    {
        vtkMath::Cross(rApex.mOriginalDirection, rApex.mParentDirection, normal);
    }

    assert(vtkMath::Norm(normal) > 1e-10);
    vtkMath::Normalize(normal);

    // Split the point cloud twice, creating new apices if needed. Should this be bundled up into a method/loop?
    {
        vtkSmartPointer<vtkPolyData> cloud = SplitPointCloud(rApex.mPointCloud,
                                                             normal,
                                                             centre,
                                                             false);

        if (cloud->GetNumberOfPoints() > 0)
        {
            double end_location[3];
            vtkIdType end_id = InsertBranch(cloud, rApex.mStartId, rApex.mOriginalDirection, end_location);

            if (std::sqrt(vtkMath::Distance2BetweenPoints(mAirwayTree->GetPoints()->GetPoint(rApex.mStartId), end_location)) > mLengthLimit
                && (unsigned)cloud->GetNumberOfPoints() > mPointLimit
                && end_id != -1)
            {
                double end_direction[3];
                vtkMath::Subtract(end_location, rApex.mCurrentLocation, end_direction);
                vtkMath::Normalize(end_direction);
                AddApex(end_id, end_location, end_direction, rApex.mOriginalDirection, rApex.mGeneration + 1);
            }
            else
            {
                InvalidateClosestPoint(end_location, cloud);
            }
        }
    }

    {
        vtkSmartPointer<vtkPolyData> cloud = SplitPointCloud(rApex.mPointCloud,
                                                             normal,
                                                             centre,
                                                             true);

        if (cloud->GetNumberOfPoints() > 0)
        {
            double end_location[3];
            vtkIdType end_id = InsertBranch(cloud, rApex.mStartId, rApex.mOriginalDirection, end_location);

            if (std::sqrt(vtkMath::Distance2BetweenPoints(mAirwayTree->GetPoints()->GetPoint(rApex.mStartId), end_location)) > mLengthLimit
                && (unsigned)cloud->GetNumberOfPoints() > mPointLimit
                && end_id != -1)
            {
                double end_direction[3];
                vtkMath::Subtract(end_location, rApex.mCurrentLocation, end_direction);
                vtkMath::Normalize(end_direction);
                AddApex(end_id, end_location, end_direction, rApex.mOriginalDirection, rApex.mGeneration + 1);
            }
            else
            {
                InvalidateClosestPoint(end_location, cloud);
            }
        }
    }
}

vtkIdType AirwayGenerator::InsertBranch(vtkSmartPointer<vtkPolyData> pPointCloud,
                                        unsigned startId,
                                        double originalDirection[3],
                                        double endLocation[3])
{
    GetCentreOfMass(pPointCloud, endLocation);

    CheckBranchAngleLengthAndAdjust(startId, originalDirection, endLocation);

    // If the branch point isn't inside the surface then terminate
    // Repeat test for "isinside" because, prior to VTK 8.2, it may fail when the point is on surface (#3002)
    if (!mPointSelector->IsInsideSurface(endLocation) && !mPointSelector->IsInsideSurface(endLocation) && !mPointSelector->IsInsideSurface(endLocation))
    {
        return -1;
    }

    // Create the new branch in mAirwayTree
    vtkIdType new_id = mAirwayTree->GetPoints()->InsertNextPoint(endLocation);

    vtkIdType pt_ids[2];
    pt_ids[0] = startId;
    pt_ids[1] = new_id;

    mAirwayTree->InsertNextCell(VTK_LINE, 2, pt_ids);

    return new_id;
}

void AirwayGenerator::InvalidateClosestPoint(double point[3], vtkSmartPointer<vtkPolyData> searchCloud)
{
    assert(searchCloud->GetNumberOfPoints() > 0);

    vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
    locator->SetDataSet(searchCloud);
    locator->BuildLocator();

    unsigned search_cloud_id = locator->FindClosestPoint(point);

    unsigned invalid_id = mSeedPointLocator->FindClosestPoint(searchCloud->GetPoints()->GetPoint(search_cloud_id));

    assert(!mInvalidIds.count(invalid_id)); //Need to modify the point splitter to ensure this.
    mInvalidIds.insert(invalid_id);
}

void AirwayGenerator::Generate()
{
    for (std::deque<AirwayGeneration>::iterator gen_iter = mGenerations.begin();
         gen_iter != mGenerations.end();
         ++gen_iter)
    {
        gen_iter->DistributeGrowthPoints(mSeedPointCloud, mInvalidIds);

        for (std::deque<Apex>::iterator apex_iter = gen_iter->GetApices().begin();
            apex_iter != gen_iter->GetApices().end();
            ++apex_iter)
        {
            GrowApex(*apex_iter);
        }
    }
}

void AirwayGenerator::CheckBranchAngleLengthAndAdjust(unsigned startId, double originalDirection[3], double centre[3])
{
    // Calculate vector from apex start to the centre
    double new_direction[3];
    double start_point[3];
    mAirwayTree->GetPoints()->GetPoint(startId, start_point);
    vtkMath::Subtract(centre, start_point, new_direction);

    // Record the branching length
    double branch_length = vtkMath::Norm(new_direction);

    vtkMath::Normalize(new_direction);

    double old_direction[3];
    std::copy(originalDirection, originalDirection+3, old_direction);
    vtkMath::Normalize(old_direction);

    // Determine branch angle & rotate if needed
    const double dot_product = vtkMath::Dot(new_direction, old_direction);
    const double branch_angle = std::acos(dot_product/(vtkMath::Norm(new_direction)*vtkMath::Norm(old_direction)));

    if (branch_angle > (M_PI*mAngleLimit/180.0))
    {
        // Rotate centre to within tolerance
        // Instead of forming an axis of rotation and a rotation matrix,
        // two cross products are taken and used to project direction
        // onto the correct angle.
        double perpendicular_one[3];
        vtkMath::Cross(new_direction, old_direction, perpendicular_one);
        vtkMath::Normalize(perpendicular_one);
        double perpendicular_two[3];
        vtkMath::Cross(old_direction, perpendicular_one, perpendicular_two);
        vtkMath::Normalize(perpendicular_two);

        vtkMath::MultiplyScalar(old_direction, std::cos(M_PI*mAngleLimit/180.0));
        vtkMath::MultiplyScalar(perpendicular_two, std::sin(M_PI*mAngleLimit/180.0));
        vtkMath::Add(perpendicular_two, old_direction, new_direction);
    }

    // Reduce the branch length as required
    branch_length = branch_length*mBranchingFraction;
    vtkMath::MultiplyScalar(new_direction, branch_length);

    // Update the centre
    centre[0] = 0, centre[1] = 0, centre[2] = 0;
    vtkMath::Add(start_point, new_direction, centre);
}

unsigned AirwayGenerator::HorsfieldProcessPoint(unsigned pointId, vtkSmartPointer<vtkDoubleArray> pOrder, std::deque<unsigned>& rProcessedPoints)
{
    rProcessedPoints.push_back(pointId);

    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> cell_ids = vtkSmartPointer<vtkIdList>::New();

    unsigned max_child_order = 1; // Terminal branches are defined to have order one
    unsigned num_children = 0;

    // Loop over the cells, finding connected points
    mAirwayTree->GetPointCells(pointId, cell_ids);
    for (int i = 0; i < cell_ids->GetNumberOfIds(); ++i)
    {
        point_ids->Reset();
        mAirwayTree->GetCellPoints(cell_ids->GetId(i), point_ids);

        for (int j = 0; j < point_ids->GetNumberOfIds(); ++j)
        {
            unsigned new_point_id = point_ids->GetId(j);

            if (!(std::find(rProcessedPoints.begin(), rProcessedPoints.end(), new_point_id) != rProcessedPoints.end()))
            {
                unsigned new_order = HorsfieldProcessPoint(new_point_id, pOrder, rProcessedPoints);

                if (new_order > max_child_order)
                {
                    max_child_order = new_order;
                }
                ++num_children;
            }
        }
    }

    // Horsfield always increments if there is more than one child
    if (num_children > 1)
    {
        ++max_child_order;
    }

    pOrder->SetValue(pointId, max_child_order);

    return max_child_order;
}

void AirwayGenerator::CalculateHorsfieldOrder()
{
    std::deque<unsigned> processed_points;

    vtkSmartPointer<vtkDoubleArray> order = vtkSmartPointer<vtkDoubleArray>::New();
    order->SetName("horsfield_order");

    order->SetNumberOfValues(mAirwayTree->GetNumberOfPoints());

    for (std::vector<unsigned>::iterator iter = mStartIds.begin();
         iter != mStartIds.end();
         ++iter)
    {
        HorsfieldProcessPoint(*iter, order, processed_points);
    }

    mAirwayTree->GetPointData()->AddArray(order);
}

void AirwayGenerator::MarkStartIds()
{
    vtkSmartPointer<vtkDoubleArray> start_ids = vtkSmartPointer<vtkDoubleArray>::New();
    start_ids->SetName("start_id");
    start_ids->SetNumberOfValues(mAirwayTree->GetNumberOfPoints());

    for (int i = 0; i < mAirwayTree->GetNumberOfPoints(); ++i)
    {
        start_ids->SetValue(i, 0.0);
    }

    for (std::vector<unsigned>::iterator iter = mStartIds.begin();
        iter != mStartIds.end();
        ++iter)
    {
        start_ids->SetValue(*iter, 1.0);
    }

    mAirwayTree->GetPointData()->AddArray(start_ids);
}

void AirwayGenerator::CalculateRadii(const double& rDiameterRatio)
{
    mDiameterRatio = rDiameterRatio;

    std::deque<unsigned> processed_points;

    vtkSmartPointer<vtkDoubleArray> radius = vtkSmartPointer<vtkDoubleArray>::New();
    radius->SetName("radius");
    radius->SetNumberOfValues(mAirwayTree->GetNumberOfPoints());

    mAirwayTree->GetPointData()->AddArray(radius);

    std::vector<unsigned>::iterator id_iter;
    std::vector<double>::iterator radius_iter;
    for (id_iter = mStartIds.begin(), radius_iter = mStartRadii.begin();
        id_iter != mStartIds.end() && radius_iter != mStartRadii.end();
        ++id_iter, ++radius_iter)
    {
        radius->SetValue(*id_iter, *radius_iter);
        RadiiProcessPoint(*id_iter, *id_iter, processed_points);
    }
}

void AirwayGenerator::RadiiProcessPoint(unsigned pointId, unsigned startId, std::deque<unsigned>& rProcessedPoints)
{
    rProcessedPoints.push_back(pointId);

    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> cell_ids = vtkSmartPointer<vtkIdList>::New();

    // Loop over the cells, finding connected points
    mAirwayTree->GetPointCells(pointId, cell_ids);
    for (int i = 0; i < cell_ids->GetNumberOfIds(); ++i)
    {
        point_ids->Reset();
        mAirwayTree->GetCellPoints(cell_ids->GetId(i), point_ids);

        for (int j = 0; j < point_ids->GetNumberOfIds(); ++j)
        {
            unsigned new_point_id = point_ids->GetId(j);

            if (!(std::find(rProcessedPoints.begin(), rProcessedPoints.end(), new_point_id) != rProcessedPoints.end()))
            {
                RadiiProcessPoint(new_point_id, startId, rProcessedPoints);
            }
        }
    }

    vtkSmartPointer<vtkDoubleArray> orders = (vtkDoubleArray*) mAirwayTree->GetPointData()->GetArray("horsfield_order");
    vtkSmartPointer<vtkDoubleArray> radii = (vtkDoubleArray*) mAirwayTree->GetPointData()->GetArray("radius");

    double current_order = orders->GetTuple1(pointId);
    double max_order = orders->GetTuple1(startId);
    double max_diameter = 2*radii->GetTuple1(startId);

    double radius = std::pow(10.0, std::log10(mDiameterRatio)*(current_order - max_order) + std::log10(max_diameter))*0.5;

    radii->SetValue(pointId, radius);
}

bool AirwayGenerator::IsInsideLobeSurface(double point[3])
{
    // Repeat test for "isinside" because, prior to VTK 8.2, it may fail when the point is on surface (#3002)
    return (mPointSelector->IsInsideSurface(point) || mPointSelector->IsInsideSurface(point) || mPointSelector->IsInsideSurface(point));
}

double AirwayGenerator::DistanceFromLobeSurface(double point[3])
{
    vtkSmartPointer<vtkCellLocator> mLobeSurfaceLocator = vtkSmartPointer<vtkCellLocator>::New();
    mLobeSurfaceLocator->SetDataSet(mLobeSurface);
    mLobeSurfaceLocator->BuildLocator();

    double closest_point[3];
    vtkIdType cell_id;
    int sub_id;
    double dist2;
    mLobeSurfaceLocator->FindClosestPoint(point, closest_point, cell_id, sub_id, dist2);

    return sqrt(dist2);
}

double AirwayGenerator::CalculateLobeVolume()
{
    vtkSmartPointer<vtkMassProperties> mass_properties = vtkSmartPointer<vtkMassProperties>::New();
#if VTK_MAJOR_VERSION >= 6
    mass_properties->SetInputData(mLobeSurface);
#else
    mass_properties->SetInput(mLobeSurface);
#endif

    return mass_properties->GetVolume();
}

void AirwayGenerator::WriteDecomposedAirways(std::string rOutputDirectory, std::string rOutputFileNameRoot)
{
    // Use a vtk connectivity filter to separate the airway tree
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivity_filter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
#if VTK_MAJOR_VERSION >= 6
    connectivity_filter->SetInputData(mAirwayTree);
#else
    connectivity_filter->SetInput(mAirwayTree);
#endif
    connectivity_filter->Update();
    connectivity_filter->SetExtractionModeToSpecifiedRegions();

    OutputFileHandler output(rOutputDirectory, false);

    // Loop over the airway components writing out each one in turn
    for (int connected_region = 0; connected_region < connectivity_filter->GetNumberOfExtractedRegions(); ++connected_region)
    {
        connectivity_filter->InitializeSpecifiedRegionList();
        connectivity_filter->AddSpecifiedRegion(connected_region);
        connectivity_filter->Update();

        std::stringstream vtu_file_name;
        vtu_file_name << output.GetOutputDirectoryFullPath() << rOutputFileNameRoot << "_" << connected_region << ".vtu";

        vtkSmartPointer<vtkCleanPolyData> poly_clean = vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION >= 6
        poly_clean->SetInputConnection(connectivity_filter->GetOutputPort());
#else
        poly_clean->SetInput(connectivity_filter->GetOutput());
#endif
        poly_clean->Update();

        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        grid->DeepCopy(poly_clean->GetOutput());
        grid->SetCells(VTK_LINE, poly_clean->GetOutput()->GetLines());

        // In extremely rare cases, the vtkAppendFilter can incorrectly merge two points that are close but aren't actually coincident.
        // This can leave repeated line elements, these are filtered out here.
        const vtkSmartPointer<vtkUnstructuredGrid>& appended_grid = grid;
        vtkSmartPointer<vtkUnstructuredGrid> filtered_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        filtered_grid->SetPoints(appended_grid->GetPoints());
        filtered_grid->SetCells(VTK_LINE, vtkSmartPointer<vtkCellArray>::New());

        std::set<std::set<int> > cell_set;
        std::set<std::set<int> >::iterator cell_set_iter;

        for (int i = 0; i < appended_grid->GetNumberOfCells(); ++i)
        {
            vtkSmartPointer<vtkLine> line = (vtkLine*) appended_grid->GetCell(i);

            vtkSmartPointer<vtkIdList> line_ids = line->GetPointIds();

            std::set<int> line_ids_set;
            line_ids_set.insert(line_ids->GetId(0));
            line_ids_set.insert(line_ids->GetId(1));

            if (cell_set.find(line_ids_set) == cell_set.end())
            {
                filtered_grid->InsertNextCell(VTK_LINE, line_ids);
                cell_set.insert(line_ids_set);
            }
        }

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtu_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        vtu_writer->SetFileName(vtu_file_name.str().c_str());
#if VTK_MAJOR_VERSION >= 6
        vtu_writer->SetInputData(filtered_grid);
#else
        vtu_writer->SetInput(filtered_grid);
#endif
        vtu_writer->Write();

        // Load the vtu in to a Chaste mesh and serialize out in triangles/tetgen format
        VtkMeshReader<1,3> combined_mesh_reader(filtered_grid);
        TetrahedralMesh<1,3> combined_mesh;
        combined_mesh.ConstructFromMeshReader(combined_mesh_reader);

        // Insert data attributes in the vtu file as node attributes in the mesh
        for (TetrahedralMesh<1,3>::NodeIterator node_iter = combined_mesh.GetNodeIteratorBegin();
             node_iter != combined_mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            double radius = appended_grid->GetPointData()->GetArray("radius")->GetTuple1(node_iter->GetIndex());
            node_iter->AddNodeAttribute(radius);
        }

        std::stringstream triangles_file_name;
        triangles_file_name << rOutputFileNameRoot << "_" << connected_region;

        TrianglesMeshWriter<1,3> combined_mesh_writer(rOutputDirectory, triangles_file_name.str(), false);
        combined_mesh_writer.WriteFilesUsingMesh(combined_mesh);
    }
}

#endif //( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

#endif //CHASTE_VTK
