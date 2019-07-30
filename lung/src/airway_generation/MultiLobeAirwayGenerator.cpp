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


#include <boost/foreach.hpp>
#include <map>

#include "MultiLobeAirwayGenerator.hpp"
#include "VtkMeshWriter.hpp"
#include "VtkMeshReader.hpp"
#include "CmguiMeshWriter.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"
#include "vtkAppendFilter.h"
#include "vtkSTLReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkCellArray.h"

#if ((VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

MultiLobeAirwayGenerator::MultiLobeAirwayGenerator(TetrahedralMesh<1,3>& rAirwaysMesh, bool pointDistanceLimit) :
                                                                                         mAirwaysMesh(rAirwaysMesh),
                                                                                         mNumberOfPointsPerLung(0),
                                                                                         mPointVolume(-1),
                                                                                         mPointDistanceLimit(pointDistanceLimit)

{
}

MultiLobeAirwayGenerator::~MultiLobeAirwayGenerator()
{
    // Deallocate the individual airway generators
    typedef std::pair<AirwayGenerator*, LungLocation> pair_type;
    for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
         generators_iter != mLobeGenerators.end();
         ++generators_iter)
    {
        delete generators_iter->first;
    }
}

void MultiLobeAirwayGenerator::AddLobe(const std::string& rFileName, LungLocation lungLocation)
{
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(rFileName.c_str());
    reader->Update();

    AddLobe(reader->GetOutput(), lungLocation);
}

void MultiLobeAirwayGenerator::AddLobe(vtkSmartPointer<vtkPolyData> pLobeSurface, LungLocation lungLocation)
{
    AirwayGenerator* p_lobe_generator = new AirwayGenerator(pLobeSurface,
                                                            mMinimumBranchLength,
                                                            mPointLimit,
                                                            mAngleLimit,
                                                            mBranchingFraction,
                                                            mPointDistanceLimit);

    mLobeGenerators.push_back(std::make_pair(p_lobe_generator, lungLocation));
}


unsigned MultiLobeAirwayGenerator::GetNumLobes(LungLocation lungLocation)
{
    unsigned lobe_count = 0u;

    typedef std::pair<AirwayGenerator*, LungLocation> pair_type;
    for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
         generators_iter != mLobeGenerators.end();
         ++generators_iter)
    {
        if (generators_iter->second == lungLocation)
        {
            lobe_count++;
        }
    }

    return lobe_count;
}

void MultiLobeAirwayGenerator::AssignGrowthApices()
{
    typedef std::pair<AirwayGenerator*, LungLocation> pair_type;

    // Create the initial growth apices for each major airways end point
    for (TetrahedralMesh<1,3>::NodeIterator iter = mAirwaysMesh.GetNodeIteratorBegin();
         iter != mAirwaysMesh.GetNodeIteratorEnd();
         ++iter)
    {
        double radius = iter->rGetNodeAttributes()[0];
        unsigned is_terminal = iter->rGetNodeAttributes()[1];

        if (is_terminal == 1 && iter->GetIndex() != 0) // Airway end point
        {
            // Determine branch direction
            Element<1,3>* p_branch = mAirwaysMesh.GetElement(*(iter->ContainingElementsBegin())); //Assume that there is only one branch attached to this node, should probably add asserts here...

            c_vector<double, 3> direction = p_branch->GetNodeLocation(0) - p_branch->GetNodeLocation(1);

            if (p_branch->GetNode(1)->rGetNodeAttributes()[1] == 1)
            {
                direction = -direction;
            }
            direction = direction/norm_2(direction);

            // Add apex to the correct generator
            double origin[3];
            origin[0] = iter->rGetLocation()[0];
            origin[1] = iter->rGetLocation()[1];
            origin[2] = iter->rGetLocation()[2];
            double dir[3];
            dir[0] = direction[0];
            dir[1] = direction[1];
            dir[2] = direction[2];

            bool end_point_assigned = false;

            for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
                 generators_iter != mLobeGenerators.end();
                 ++generators_iter)
            {
                if (generators_iter->first->IsInsideLobeSurface(origin))
                {
                    end_point_assigned = true;
                    vtkSmartPointer<vtkPolyData> cloud = vtkSmartPointer<vtkPolyData>::New();
                    cloud->SetPoints(vtkSmartPointer<vtkPoints>::New());

                    //\todo This needs to be given the correct generation number.
                    // Unfortunately, this isn't available in the data and will have to be calculated.
                    generators_iter->first->AddInitialApex(origin, dir, dir, radius, 0);
                }
            }

            if (!end_point_assigned) // If the end point is not contained in any lobe assign it to the nearest
            {
                double dist_min = DBL_MAX;
                AirwayGenerator* p_generator = nullptr;

                for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
                     generators_iter != mLobeGenerators.end();
                     ++generators_iter)
                {
                    double dist = generators_iter->first->DistanceFromLobeSurface(origin);
                    if (dist < dist_min)
                    {
                        dist_min = dist;
                        p_generator = generators_iter->first;
                    }
                }
                p_generator->AddInitialApex(origin, dir, dir, radius, 0);
            }
        }
    }
}

void MultiLobeAirwayGenerator::DistributePoints()
{
    if (mNumberOfPointsPerLung == 0u && mPointVolume == -1)
    {
        EXCEPTION("Must call SetNumberOfPointsPerLung or SetPointVolume before distributing points.");
    }

    if (mNumberOfPointsPerLung > 0u && mPointVolume == -1)
    {
        // Calculate the total volume of each lung
        std::map<LungLocation, double> lung_volumes;
        lung_volumes[LEFT] = 0.0;
        lung_volumes[RIGHT] = 0.0;

        typedef std::pair<AirwayGenerator*, LungLocation> pair_type;
        for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
             generators_iter != mLobeGenerators.end();
             ++generators_iter)
        {
            lung_volumes[generators_iter->second] += generators_iter->first->CalculateLobeVolume();
        }

        // Distribute the initial point clouds
        for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
             generators_iter != mLobeGenerators.end();
             ++generators_iter)
        {
            unsigned num_points =  static_cast<unsigned>((generators_iter->first->CalculateLobeVolume() / lung_volumes[generators_iter->second]) * mNumberOfPointsPerLung + 0.5);
            generators_iter->first->CreatePointCloudUsingTargetPoints(num_points);
        }
    }
    else if (mNumberOfPointsPerLung == 0u && mPointVolume > 0.0)
    {
        typedef std::pair<AirwayGenerator*, LungLocation> pair_type;
        for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
             generators_iter != mLobeGenerators.end();
             ++generators_iter)
        {
            generators_iter->first->CreatePointCloudUsingTargetVolume(mPointVolume);
        }
    }
    else
    {
        EXCEPTION("Both SetNumberOfPointsPerLung and SetPointVolume called. Please use one or the other.");
    }
}

void MultiLobeAirwayGenerator::SetNumberOfPointsPerLung(const unsigned& rNumberOfPointsPerLung)
{
    mNumberOfPointsPerLung = rNumberOfPointsPerLung;
}

void MultiLobeAirwayGenerator::SetPointVolume(const double& rVolume )
{
    mPointVolume = rVolume;
}

void MultiLobeAirwayGenerator::SetMinimumBranchLength(const double& rMinimumbranchLength)
{
    mMinimumBranchLength = rMinimumbranchLength;
}

void MultiLobeAirwayGenerator::SetPointLimit(const unsigned& rPointLimit)
{
    mPointLimit = rPointLimit;
}

void MultiLobeAirwayGenerator::SetAngleLimit(const double& rAngleLimit)
{
    mAngleLimit = rAngleLimit;
}

void MultiLobeAirwayGenerator::SetBranchingFraction(const double& rBranchingFraction)
{
    mBranchingFraction = rBranchingFraction;
}

void  MultiLobeAirwayGenerator::SetDiameterRatio(const double& rDiameterRatio)
{
    mDiameterRatio = rDiameterRatio;
}

void MultiLobeAirwayGenerator::Generate(std::string rOutputDirectory, std::string rBaseName)
{
    vtkSmartPointer<vtkAppendFilter> append_filter = vtkSmartPointer<vtkAppendFilter>::New();

    // Merge points cannot be set in vtk5.6 but defaults to on
    // In vtk5.8 and higher we must explicitly set it to be on
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 8) || VTK_MAJOR_VERSION >= 6)
    append_filter->MergePointsOn();
#endif

    // Merge in the major airways, the mesh has to be converted to a vtk unstructured grid first
    // We use the Chaste VtkMeshWriter to write the mesh to disk then load it back in as a vtu.
    ///\todo It would be much cleaner to do this in memory - might require a major refactor of VtkMeshWriter
    VtkMeshWriter<1,3> mesh_writer(rOutputDirectory, "major_airways", false);

    std::vector<double> order; //The Horsfield order of nodes in the tree
    std::vector<double> radii; //The radius of nodes in the tree
    std::vector<double> start_ids; //Indication of node type in the tree; imaging: 2, generated start: 1, generated: 0

    for (TetrahedralMesh<1,3>::NodeIterator node_iter = mAirwaysMesh.GetNodeIteratorBegin();
         node_iter != mAirwaysMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        order.push_back(20.0);
        radii.push_back(node_iter->rGetNodeAttributes()[0]); ///\todo Magic number, allow the user to specify the attribute index

        if (node_iter->rGetNodeAttributes()[1] == 0.0)
        {
            start_ids.push_back(2.0);
        }
        else if (node_iter->rGetNodeAttributes()[1] == 1.0)
        {
            start_ids.push_back(1.0);
        }
        else
        {
            EXCEPTION("The second node attribute in the major airways mesh must be 0.0 for non-terminal or 1.0 for terminal.");
        }
    }

    mesh_writer.AddPointData("horsfield_order", order);
    mesh_writer.AddPointData("radius", radii);
    mesh_writer.AddPointData("start_id", start_ids);
    mesh_writer.WriteFilesUsingMesh(mAirwaysMesh);

    OutputFileHandler output(rOutputDirectory, false);
    std::string major_airways_file_name = output.GetOutputDirectoryFullPath() + "major_airways.vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridReader> major_airways_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    major_airways_reader->SetFileName(major_airways_file_name.c_str());
    major_airways_reader->Update();

#if VTK_MAJOR_VERSION >= 6
    append_filter->AddInputData(major_airways_reader->GetOutput());
#else
    append_filter->AddInput(major_airways_reader->GetOutput());
#endif

    // Loop over generators, generate and merge the results
    typedef std::pair<AirwayGenerator*, LungLocation> pair_type;
    for (std::vector<pair_type>::iterator generators_iter = mLobeGenerators.begin();
         generators_iter != mLobeGenerators.end();
         ++generators_iter)
    {
        generators_iter->first->Generate();
        generators_iter->first->CalculateHorsfieldOrder();
        generators_iter->first->CalculateRadii(mDiameterRatio);
        generators_iter->first->MarkStartIds();

#if VTK_MAJOR_VERSION >= 6
        append_filter->AddInputData(generators_iter->first->GetAirwayTree());
#else
        append_filter->AddInput(generators_iter->first->GetAirwayTree());
#endif
    }

    append_filter->Update();

    // In extremely rare cases, the append_filter can incorrectly merge two points that are close but aren't actually coincident.
    // This can leave repeated line elements, these are filtered out here.
    vtkSmartPointer<vtkUnstructuredGrid> appended_grid = append_filter->GetOutput();
    vtkSmartPointer<vtkUnstructuredGrid> filtered_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    filtered_grid->SetPoints(appended_grid->GetPoints());

#if VTK_MAJOR_VERSION == 6
    /* Note that there is a bug in some versions of VTK 6 which involves filtering for
     * duplicate points and attributes but not altering the size of the attribute vectors.
     * This leads to attribute vectors which contain unintialised values.  See
     * "BUG #15746: Fixes issues with Merge Points in vtkAppendFilter"
     */
    unsigned num_points = appended_grid->GetNumberOfPoints();
    appended_grid->GetPointData()->GetArray("radius")->SetNumberOfTuples( num_points );
    appended_grid->GetPointData()->GetArray("horsfield_order")->SetNumberOfTuples( num_points );
    appended_grid->GetPointData()->GetArray("start_id")->SetNumberOfTuples( num_points );
#endif

    filtered_grid->GetPointData()->AddArray(appended_grid->GetPointData()->GetArray("radius"));
    filtered_grid->GetPointData()->AddArray(appended_grid->GetPointData()->GetArray("horsfield_order"));
    filtered_grid->GetPointData()->AddArray(appended_grid->GetPointData()->GetArray("start_id"));
    filtered_grid->SetCells(VTK_LINE, vtkSmartPointer<vtkCellArray>::New());

    std::set<std::set<int> > cell_set;

    for (unsigned i = 0; i < (unsigned) appended_grid->GetNumberOfCells(); ++i)
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

    // Write the merged grid to disk
    std::string output_file_name = output.GetOutputDirectoryFullPath() + rBaseName + ".vtu";

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> combined_vtu_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    combined_vtu_writer->SetFileName(output_file_name.c_str());
#if VTK_MAJOR_VERSION >= 6
    combined_vtu_writer->SetInputData(filtered_grid);
#else
    combined_vtu_writer->SetInput(filtered_grid);
#endif
    combined_vtu_writer->Write();

    // Load the vtu in to a Chaste mesh and serialize out in triangles/tetgen format
    VtkMeshReader<1,3> combined_mesh_reader(filtered_grid);
    TetrahedralMesh<1,3> combined_mesh;
    combined_mesh.ConstructFromMeshReader(combined_mesh_reader);

    // Insert data attributes in the vtu file as node attributes in the mesh
    for (TetrahedralMesh<1,3>::NodeIterator node_iter = combined_mesh.GetNodeIteratorBegin();
         node_iter != combined_mesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        double radius = append_filter->GetOutput()->GetPointData()->GetArray("radius")->GetTuple1(node_iter->GetIndex());
        node_iter->AddNodeAttribute(radius);

        double start_id = append_filter->GetOutput()->GetPointData()->GetArray("start_id")->GetTuple1(node_iter->GetIndex());
        node_iter->AddNodeAttribute(start_id);
    }

    TrianglesMeshWriter<1,3> combined_mesh_writer(rOutputDirectory, rBaseName, false);
    combined_mesh_writer.WriteFilesUsingMesh(combined_mesh);

    CmguiMeshWriter<1,3> cmgui_writer(rOutputDirectory, rBaseName, false);
    cmgui_writer.WriteFilesUsingMesh(combined_mesh);
}

#endif // (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6

#endif //CHASTE_VTK
