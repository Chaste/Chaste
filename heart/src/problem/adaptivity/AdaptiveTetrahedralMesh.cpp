/*

Copyright (C) Fujitsu Laboratories of Europe, 2009-2010

*/

/*

Copyright (c) 2005-2013, University of Oxford.
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


#define COVERAGE_IGNORE /// \todo #2367 We no longer test adaptivity library
#ifdef CHASTE_ADAPTIVITY

#include "AdaptiveTetrahedralMesh.hpp"
#include "HeartConfig.hpp"
#include "OutputFileHandler.hpp"

AdaptiveTetrahedralMesh::AdaptiveTetrahedralMesh() :
    mNumNodes(0),
    mNumElements(0),
    mNumLocalNodes(0),
    mAdaptSuccess(false),
    mpDiscreteGeometryConstraints(NULL),
    mpErrorMeasure(NULL),
    mpAdapt(NULL),
    mGoodEdgeRange(0.0),
    mBadEdgeCriterion(0.0),
    mVerbose(false)

{
        mpVtkUnstructuredGrid = vtkUnstructuredGrid::New();
}

AdaptiveTetrahedralMesh::~AdaptiveTetrahedralMesh()
{
    Reset();//Delete multiple-use pointers
    mpVtkUnstructuredGrid->Delete();
}

void AdaptiveTetrahedralMesh::ConstructFromVtuFile(std::string fileName)
{
    vtkXMLUnstructuredGridReader *p_vtk_reader = vtkXMLUnstructuredGridReader::New();
    p_vtk_reader->SetFileName(fileName.c_str());
    p_vtk_reader->Update();

    mpVtkUnstructuredGrid->DeepCopy(p_vtk_reader->GetOutput());
    mpVtkUnstructuredGrid->Update();

    p_vtk_reader->Delete();

    mNumNodes = mpVtkUnstructuredGrid->GetNumberOfPoints();
    mNumElements = mpVtkUnstructuredGrid->GetNumberOfCells();
    mNumLocalNodes = mNumNodes;       // Since each process has a complete copy of the vtkUnstructuredGrid
}

void AdaptiveTetrahedralMesh::ConstructFromMesh(AbstractTetrahedralMesh<3,3>* rMesh)
{
    vtkPoints *p_pts = vtkPoints::New();
    p_pts->SetDataTypeToDouble();
    for (unsigned i=0; i<rMesh->GetNumNodes(); i++)
    {
        p_pts->InsertPoint(i, rMesh->GetNode(i)->rGetLocation().data());
    }

    mpVtkUnstructuredGrid->Allocate(rMesh->GetNumNodes(), rMesh->GetNumNodes());
    mpVtkUnstructuredGrid->SetPoints(p_pts);
    mpVtkUnstructuredGrid->Update();
    p_pts->Delete();    // Reference counted

    for (unsigned i=0; i<rMesh->GetNumElements(); i++)
    {
        vtkTetra *p_tetra = vtkTetra::New();
        for (int j = 0; j < 4; ++j)
        {
            p_tetra->GetPointIds()->SetId(j, rMesh->GetElement(i)->GetNodeGlobalIndex(j));
        }
        mpVtkUnstructuredGrid->InsertNextCell(p_tetra->GetCellType(), p_tetra->GetPointIds());
        mpVtkUnstructuredGrid->Update();
        p_tetra->Delete();    // Reference counted
    }

    mNumNodes = mpVtkUnstructuredGrid->GetNumberOfPoints();
    mNumElements = mpVtkUnstructuredGrid->GetNumberOfCells();
    mNumLocalNodes = mNumNodes;       // Since each process has a complete copy of the vtkUnstructuredGrid
}

void AdaptiveTetrahedralMesh::ConstructFromDistributedMesh(DistributedTetrahedralMesh<3,3>* rMesh)
{
    vtkPoints *p_pts = vtkPoints::New();
    p_pts->SetDataTypeToDouble();

    std::vector<unsigned> global_node_numbers, halo_nodes;
    std::map<unsigned, unsigned> global_to_local_index_map;

//    vtkUnsignedIntArray *p_scalars = vtkUnsignedIntArray::New();
//    p_scalars->SetName("GlobalNodeNumbers");

    unsigned index = 0;
    for (DistributedTetrahedralMesh<3,3>::NodeIterator it=rMesh->GetNodeIteratorBegin();
         it != rMesh->GetNodeIteratorEnd();
         ++it)
    {
        p_pts->InsertPoint(index, it->rGetLocation().data());
        global_node_numbers.push_back(it->GetIndex());
        global_to_local_index_map[it->GetIndex()] = index;
        index++;
    }

    rMesh->GetHaloNodeIndices(halo_nodes);
    for(unsigned i=0; i<halo_nodes.size(); i++)
    {
        global_node_numbers.push_back(halo_nodes[i]);
        p_pts->InsertPoint(index, rMesh->GetNodeOrHaloNode(halo_nodes[i])->rGetLocation().data());
        global_to_local_index_map[halo_nodes[i]] = index;
        index++;
    }

    AddPointData("GlobalNodeNumbers", global_node_numbers);

    mpVtkUnstructuredGrid->Allocate(global_node_numbers.size(), global_node_numbers.size());
    mpVtkUnstructuredGrid->SetPoints(p_pts);
    mpVtkUnstructuredGrid->Update();
    p_pts->Delete();    // Reference counted

    for (DistributedTetrahedralMesh<3,3>::ElementIterator it=rMesh->GetElementIteratorBegin();
         it != rMesh->GetElementIteratorEnd();
         ++it)
    {
        vtkTetra *p_tetra = vtkTetra::New();
        for (int j = 0; j < 4; ++j)
        {
            p_tetra->GetPointIds()->SetId(j, global_to_local_index_map[it->GetNodeGlobalIndex(j)]);
        }
        mpVtkUnstructuredGrid->InsertNextCell(p_tetra->GetCellType(), p_tetra->GetPointIds());
        mpVtkUnstructuredGrid->Update();
        p_tetra->Delete();    // Reference counted
    }

    mNumNodes = rMesh->GetNumNodes();
    mNumElements = rMesh->GetNumElements();
    mNumLocalNodes = mpVtkUnstructuredGrid->GetNumberOfPoints() - halo_nodes.size();
}

void AdaptiveTetrahedralMesh::AddPointData(std::string dataName, std::vector<double> dataPayload)
{
    vtkDoubleArray *p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    mpVtkUnstructuredGrid->GetPointData()->AddArray(p_scalars);
    mpVtkUnstructuredGrid->Update();
    p_scalars->Delete(); // Reference counted
}

void AdaptiveTetrahedralMesh::AddPointData(std::string dataName, std::vector<unsigned> dataPayload)
{
    vtkUnsignedIntArray *p_scalars = vtkUnsignedIntArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    mpVtkUnstructuredGrid->GetPointData()->AddArray(p_scalars);
    mpVtkUnstructuredGrid->Update();
    p_scalars->Delete(); // Reference counted
}

void AdaptiveTetrahedralMesh::RemoveArray(std::string dataName)
{
    mpVtkUnstructuredGrid->GetPointData()->RemoveArray(dataName.c_str());
}

void AdaptiveTetrahedralMesh::WriteMeshToFile(std::string directory, std::string fileName)

{
    std::string vtk_file_name = directory + fileName;

    vtkXMLUnstructuredGridWriter *vtk_writer = vtkXMLUnstructuredGridWriter::New();
    //Uninitialised stuff arises (see #1079), but you can remove
    //valgrind problems by removing compression:
    // **** REMOVE WITH CAUTION *****
    vtk_writer->SetCompressor(NULL);
    // **** REMOVE WITH CAUTION *****
    vtk_writer->SetFileName(vtk_file_name.c_str());
    vtk_writer->SetInput(mpVtkUnstructuredGrid);
    vtk_writer->Write();
    vtk_writer->Delete();
}

void AdaptiveTetrahedralMesh::WriteMeshToDistributedFile(std::string directory, std::string fileName)
{
    std::string vtk_file_name = directory + fileName;

    vtkXMLPUnstructuredGridWriter *vtk_writer = vtkXMLPUnstructuredGridWriter::New();
    //Uninitialised stuff arises (see #1079), but you can remove
    //valgrind problems by removing compression:
    // **** REMOVE WITH CAUTION *****
    vtk_writer->SetCompressor(NULL);
    // **** REMOVE WITH CAUTION *****
    vtk_writer->SetFileName(vtk_file_name.c_str());
    vtk_writer->SetDataModeToBinary();

    vtk_writer->SetNumberOfPieces(PetscTools::GetNumProcs());
    vtk_writer->SetGhostLevel(1);
    vtk_writer->SetStartPiece(PetscTools::GetMyRank());
    vtk_writer->SetEndPiece(PetscTools::GetMyRank());

    vtk_writer->SetInput(mpVtkUnstructuredGrid);
    vtk_writer->Write();
    vtk_writer->Delete();
}

vtkUnstructuredGrid* AdaptiveTetrahedralMesh::GetVtkUnstructuredGrid()
{
    return mpVtkUnstructuredGrid;
}

void AdaptiveTetrahedralMesh::SetAdaptCriterion(double range, double criterion)
{
    mGoodEdgeRange = range;
    mBadEdgeCriterion = criterion;
}

unsigned AdaptiveTetrahedralMesh::GetNumNodes()
{
    return mNumNodes;
}

unsigned AdaptiveTetrahedralMesh::GetNumLocalNodes()
{
    return mNumLocalNodes;
}

unsigned AdaptiveTetrahedralMesh::GetNumLocalAndHaloNodes()
{
    return mpVtkUnstructuredGrid->GetNumberOfPoints();
}

unsigned AdaptiveTetrahedralMesh::GetNumElements()
{
    return mNumElements;
}

unsigned AdaptiveTetrahedralMesh::GetNumLocalElements()
{
    return mpVtkUnstructuredGrid->GetNumberOfCells();
}

unsigned AdaptiveTetrahedralMesh::GetNumSurfaceElements()
{
    return sids.size();
}

void AdaptiveTetrahedralMesh::CalculateSENListAndSids(double coplanarTolerance)
{
    DiscreteGeometryConstraints constraints;

    if (mVerbose) constraints.verbose_on();
    constraints.set_coplanar_tolerance( coplanarTolerance );
    constraints.set_volume_input(mpVtkUnstructuredGrid);

    constraints.get_coplanar_ids(sids);
    constraints.get_surface(SENList);
    assert(sids.size()*3==SENList.size());
}

void AdaptiveTetrahedralMesh::GetGeometryConstraints()
{
    assert(mpDiscreteGeometryConstraints==NULL);
    mpDiscreteGeometryConstraints = new DiscreteGeometryConstraints;

    if (mVerbose) mpDiscreteGeometryConstraints->verbose_on();
    mpDiscreteGeometryConstraints->set_surface_input(mpVtkUnstructuredGrid, SENList, sids);
    mpDiscreteGeometryConstraints->get_constraints(max_len);
//    mpDiscreteGeometryConstraints->write_vtk(std::string("sids.vtu"));
}

void AdaptiveTetrahedralMesh::CalculateErrorMetric()
{
    double error           = HeartConfig::Instance()->GetTargetErrorForAdaptivity();
    double sigma           = HeartConfig::Instance()->GetSigmaForAdaptivity();
    double max_edge_length = HeartConfig::Instance()->GetMaxEdgeLengthForAdaptivity();
    double min_edge_length = HeartConfig::Instance()->GetMinEdgeLengthForAdaptivity();
    double gradation       = HeartConfig::Instance()->GetGradationForAdaptivity();
    unsigned max_nodes     = HeartConfig::Instance()->GetMaxNodesForAdaptivity();
    assert(mpErrorMeasure == NULL);
    mpErrorMeasure = new ErrorMeasure;

    if (mVerbose) mpErrorMeasure->verbose_on();
    mpErrorMeasure->set_input(mpVtkUnstructuredGrid);
    mpErrorMeasure->add_field("Vm", error, false, sigma);
    mpErrorMeasure->set_max_length(max_edge_length);
    mpErrorMeasure->set_max_length(&(max_len[0]), mpVtkUnstructuredGrid->GetNumberOfPoints());
    mpErrorMeasure->set_min_length(min_edge_length);
    mpErrorMeasure->apply_gradation(gradation);
    mpErrorMeasure->set_max_nodes(max_nodes);

    mpErrorMeasure->diagnostics();
}

double AdaptiveTetrahedralMesh::GetEdgeLengthDistribution(double range)
{
    //This is a temporary adaptivity class
    Adaptivity an_adapter;
    an_adapter.set_from_vtk(mpVtkUnstructuredGrid, true);
     return an_adapter.edgeLengthDistribution(range);
}

void AdaptiveTetrahedralMesh::Adapt()
{
    mpAdapt = new Adaptivity;
    mAdaptSuccess = false;

    if (mVerbose) mpAdapt->verbose_on();
    mpAdapt->set_from_vtk(mpVtkUnstructuredGrid, true);
    mpAdapt->set_adapt_sweeps(HeartConfig::Instance()->GetNumberOfAdaptiveSweeps());
    mpAdapt->set_surface_mesh(SENList);
    mpAdapt->set_surface_ids(sids);

    if ( mpAdapt->edgeLengthDistribution( mGoodEdgeRange ) > mBadEdgeCriterion )
    {
        mpAdapt->adapt();

        mpAdapt->get_surface_ids(sids);
        mpAdapt->get_surface_mesh(SENList);

        vtkUnstructuredGrid *p_new_vtk_unstructured_grid = mpAdapt->get_adapted_vtu();
        mpVtkUnstructuredGrid->DeepCopy( p_new_vtk_unstructured_grid );
        mpVtkUnstructuredGrid->Update();
        p_new_vtk_unstructured_grid->Delete();

        mAdaptSuccess = true;
    }
    else
    {
        /// \todo I can't see the point of this, since the VTK mesh is discarded (but not garbage collected)
        //mpAdapt->get_adapted_vtu()->Initialize();
    }

    mpVtkUnstructuredGrid->GetPointData()->RemoveArray("metric");
    mpVtkUnstructuredGrid->Squeeze();

    // Need to free these pointers as fresh instances are required for the next adapt (else odd things happen)
    Reset();

    // Update private member variables - note: these assume that the adapt happens sequentially
    mNumNodes = mpVtkUnstructuredGrid->GetNumberOfPoints();
    mNumElements = mpVtkUnstructuredGrid->GetNumberOfCells();
    mNumLocalNodes = mNumNodes;       // Since each process has a complete copy of the vtkUnstructuredGrid
}

void AdaptiveTetrahedralMesh::AdaptMesh()
{
    GetGeometryConstraints();

    CalculateErrorMetric();

    RemoveArray("mean_desired_lengths");
    RemoveArray("desired_lengths");
//    SetAdaptCriterion( mGoodEdgeRange, mBadEdgeCriterion );
    Adapt();
}

void AdaptiveTetrahedralMesh::Reset()
{
    delete mpDiscreteGeometryConstraints;
    delete mpErrorMeasure;
    delete mpAdapt;
    mpDiscreteGeometryConstraints=NULL;
    mpErrorMeasure=NULL;
    mpAdapt=NULL;
}

bool AdaptiveTetrahedralMesh::GetAdaptSuccess()
{
    return mAdaptSuccess;
}

void AdaptiveTetrahedralMesh::MakeVerbose(bool verbose)
{
    mVerbose = verbose;
}

#endif // CHASTE_ADAPTIVITY
#undef COVERAGE_IGNORE /// \todo #2367 We no longer test adaptivity library
