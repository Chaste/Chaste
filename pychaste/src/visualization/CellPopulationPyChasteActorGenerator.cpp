/*

Copyright (c) 2005-2024, University of Oxford.
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

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkGlyph2D.h>
#include <vtkCubeAxesActor2D.h>
#include <vtkImageData.h>
#include <vtkGeometryFilter.h>
#include <vtkTubeFilter.h>
#include <vtkExtractEdges.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkPolygon.h>
#include <vtkIdList.h>
#include <vtkFeatureEdges.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkVoxel.h>
#include <vtkPixel.h>
#include <vtkThreshold.h>
#include <vtkMarchingCubes.h>
#include <vtkMarchingSquares.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkFeatureEdges.h>
#include "UblasIncludes.hpp"
#include "UblasVectorInclude.hpp"
#include "Exception.hpp"
#include "CellLabel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellPopulationPyChasteActorGenerator.hpp"
#include "Debug.hpp"


template<unsigned DIM>
CellPopulationPyChasteActorGenerator<DIM>::CellPopulationPyChasteActorGenerator()
    : AbstractPyChasteActorGenerator<DIM>(),
      mpCellPopulation(),
      mShowMutableMeshEdges(false),
      mShowVoronoiMeshEdges(true),
      mShowPottsMeshEdges(false),
      mShowPottsMeshOutlines(false),
      mColorByCellType(false),
      mColorByCellData(false),
      mColorByCellMutationState(false),
      mColorByCellLabel(false),
      mShowCellCentres(false),
      mColorCellByUserDefined(false)
{

}

template<unsigned DIM>
CellPopulationPyChasteActorGenerator<DIM>::~CellPopulationPyChasteActorGenerator()
{

}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::AddCaBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer)
{
    auto p_potts_grid = vtkSmartPointer<vtkImageData>::New();
    auto p_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();

    boost::shared_ptr<CaBasedCellPopulation<DIM> > p_ca_population =
        boost::dynamic_pointer_cast<CaBasedCellPopulation<DIM> >(mpCellPopulation);

    if(p_ca_population && mShowPottsMeshEdges)
    {
        auto p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->GetData()->SetName("Vertex positions");

        unsigned counter = 0;
        c_vector<double, DIM> old_loc = zero_vector<double>(DIM);
        double spacing = 1.0;
        for (typename PottsMesh<DIM>::NodeIterator node_iter = p_ca_population->rGetMesh().GetNodeIteratorBegin();
             node_iter != p_ca_population->rGetMesh().GetNodeIteratorEnd(); ++node_iter)
        {
            c_vector<double, DIM> current_item = node_iter->rGetLocation();
            if (DIM == 3)
            {
                p_points->InsertNextPoint(current_item[0], current_item[1], current_item[2]);
            }
            else if (DIM == 2)
            {
                p_points->InsertNextPoint(current_item[0], current_item[1], 0.0);
            }
            else // (DIM == 1)
            {
                p_points->InsertNextPoint(current_item[0], 0.0, 0.0);
            }

            if(counter == 0)
            {
                old_loc = current_item;
            }
            if(counter == 1)
            {
                spacing = norm_2(current_item - old_loc);
            }
            counter ++;
        }

        auto p_temp_polydata = vtkSmartPointer<vtkPolyData>::New();
        p_temp_polydata->SetPoints(p_points);

        double bounds[6];
        p_temp_polydata->GetBounds(bounds);
        auto p_ca_image = vtkSmartPointer<vtkImageData>::New();
        p_ca_image->SetDimensions(std::floor((bounds[1]-bounds[0])/spacing) + 1,
                std::floor((bounds[3]-bounds[2])/spacing) + 1,
                std::floor((bounds[5]-bounds[4])/spacing) + 1);
        p_ca_image->SetOrigin(bounds[0], bounds[2], bounds[4]);
        p_ca_image->SetSpacing(spacing, spacing, spacing);

        p_geom_filter->SetInputData(p_ca_image);

        auto p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_mapper->SetInputConnection(p_geom_filter->GetOutputPort());

        auto p_volume_actor = vtkSmartPointer<vtkActor>::New();
        p_volume_actor->SetMapper(p_mapper);
        p_volume_actor->GetProperty()->SetEdgeVisibility(this->mShowEdges);
        p_volume_actor->GetProperty()->SetLineWidth(this->mEdgeSize);
        p_volume_actor->GetProperty()->SetOpacity(0.6);
        pRenderer->AddActor(p_volume_actor);
    }
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::AddPottsBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer)
{
    auto p_potts_grid = vtkSmartPointer<vtkImageData>::New();

    boost::shared_ptr<PottsBasedCellPopulation<DIM> > p_potts_population =
            boost::dynamic_pointer_cast<PottsBasedCellPopulation<DIM> >(mpCellPopulation);

    if(p_potts_population && mShowPottsMeshEdges)
    {
        auto p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->GetData()->SetName("Vertex positions");

        unsigned counter = 0;
        c_vector<double, DIM> old_loc = zero_vector<double>(DIM);
        double spacing = 1.0;

        for (typename PottsMesh<DIM>::NodeIterator node_iter = p_potts_population->rGetMesh().GetNodeIteratorBegin();
             node_iter != p_potts_population->rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            c_vector<double, DIM> current_item = node_iter->rGetLocation();
            if (DIM == 3)
            {
                p_points->InsertNextPoint(current_item[0], current_item[1], current_item[2]);
            }
            else if (DIM == 2)
            {
                p_points->InsertNextPoint(current_item[0], current_item[1], 0.0);
            }
            else // (DIM == 1)
            {
                p_points->InsertNextPoint(current_item[0], 0.0, 0.0);
            }
            if(counter == 0)
            {
                old_loc = current_item;
            }
            if(counter == 1)
            {
                spacing = norm_2(current_item - old_loc);
            }
            counter ++;
        }

        auto p_temp_polydata = vtkSmartPointer<vtkPolyData>::New();
        p_temp_polydata->SetPoints(p_points);

        double bounds[6];
        p_temp_polydata->GetBounds(bounds);
        p_potts_grid = vtkSmartPointer<vtkImageData>::New();

        // Important: We color VTK cells, not points. We add a VTK cell for each Chaste node.
        p_potts_grid->SetDimensions(std::floor((bounds[1]-bounds[0])/spacing) + 2,
                std::floor((bounds[3]-bounds[2])/spacing) + 2,
                std::floor((bounds[5]-bounds[4])/spacing) + 2);
        p_potts_grid->SetOrigin(bounds[0]-spacing/2.0, bounds[2]-spacing/2.0, bounds[4]-spacing/2.0);
        p_potts_grid->SetSpacing(spacing, spacing, spacing);

        auto p_element_ids = vtkSmartPointer<vtkDoubleArray>::New();
        p_element_ids->SetNumberOfTuples(p_potts_grid->GetNumberOfPoints());
        p_element_ids->SetName("Cell Id");

        auto p_element_base_ids = vtkSmartPointer<vtkDoubleArray>::New();
        p_element_base_ids->SetNumberOfTuples(p_potts_grid->GetNumberOfPoints());
        p_element_base_ids->SetName("Cell Base Id");
        for(unsigned idx=0; idx<p_potts_grid->GetNumberOfPoints(); idx++)
        {
            p_element_ids->SetTuple1(idx, -1.0);
        }

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
             cell_iter != mpCellPopulation->End(); ++cell_iter)
        {
            PottsElement<DIM>* p_element = p_potts_population->GetElementCorrespondingToCell(*cell_iter);

            for(unsigned idx=0; idx<p_element->GetNumNodes(); idx++)
            {
                unsigned node_index = p_element->GetNode(idx)->GetIndex();

                if(mColorByCellType)
                {
                    p_element_ids->InsertNextTuple1((*cell_iter)->GetCellProliferativeType()->GetColour());
                }
                else if(mColorByCellData && !this->mDataLabel.empty())
                {
                    std::vector<std::string> keys = (*cell_iter)->GetCellData()->GetKeys();
                    if (std::find(keys.begin(), keys.end(), this->mDataLabel) != keys.end())
                    {
                        p_element_ids->InsertNextTuple1((*cell_iter)->GetCellData()->GetItem(this->mDataLabel));
                    }
                    else
                    {
                        p_element_ids->InsertNextTuple1(0.0);
                    }
                }

                else if(mColorByCellMutationState)
                {
                    double mutation_state = (*cell_iter)->GetMutationState()->GetColour();
                    CellPropertyCollection collection = (*cell_iter)->rGetCellPropertyCollection();
                    CellPropertyCollection label_collection = collection.GetProperties<CellLabel>();

                    if (label_collection.GetSize() == 1)
                    {
                        boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(label_collection.GetProperty());
                        mutation_state = p_label->GetColour();
                    }
                    p_element_ids->InsertNextTuple1(mutation_state);
                }
                else if(mColorByCellLabel)
                {
                    double label = 0.0;
                    if ((*cell_iter)->template HasCellProperty<CellLabel>())
                    {
                        CellPropertyCollection collection = (*cell_iter)->rGetCellPropertyCollection().template  GetProperties<CellLabel>();
                        boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
                        label = p_label->GetColour();
                    }
                    p_element_ids->InsertNextTuple1(label);
                }
                else
                {
                    p_element_ids->SetTuple1(node_index, double((*cell_iter)->GetCellId()+1));

                }
                p_element_base_ids->SetTuple1(node_index, double((*cell_iter)->GetCellId()+1));
            }
        }
        p_potts_grid->GetCellData()->SetScalars(p_element_ids);
        p_potts_grid->GetCellData()->AddArray(p_element_ids);
        p_potts_grid->GetCellData()->AddArray(p_element_base_ids);

        auto p_scaled_ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
        if(!mColorByCellData)
        {
            double range[2];
            p_potts_grid->GetCellData()->GetArray("Cell Id")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpDiscreteColorTransferFunction->GetColor(double(idx)/255.0, color);
                p_scaled_ctf->AddRGBPoint(double(idx)*range[1]/255.0, color[0], color[1], color[2]);
            }
        }

        auto p_geometry_filter_pre = vtkSmartPointer<vtkGeometryFilter>::New();
        p_geometry_filter_pre->SetInputData(p_potts_grid);

        auto p_threshold = vtkSmartPointer<vtkThreshold>::New();
        p_threshold->SetInputConnection(p_geometry_filter_pre->GetOutputPort());

#if VTK_MAJOR_VERSION < 9
        p_threshold->ThresholdByUpper(0.0);
#else
        p_threshold->SetUpperThreshold(0.0);
#endif
        p_threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Cell Id");

        auto p_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();
        p_geom_filter->SetInputConnection(p_threshold->GetOutputPort());

        auto p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_mapper->SetInputConnection(p_geom_filter->GetOutputPort());
        p_mapper->SetLookupTable(p_scaled_ctf);
        p_mapper->ScalarVisibilityOn();
        p_mapper->SelectColorArray("Cell Id");
        p_mapper->SetScalarModeToUseCellData();
        p_mapper->SetColorModeToMapScalars();

        auto p_volume_actor = vtkSmartPointer<vtkActor>::New();
        p_volume_actor->SetMapper(p_mapper);
        //p_volume_actor->GetProperty()->SetEdgeVisibility(this->mShowEdges);
        p_volume_actor->GetProperty()->SetLineWidth(this->mEdgeSize);
        p_volume_actor->GetProperty()->SetOpacity(this->mVolumeOpacity);
        if(mColorCellByUserDefined)
        {
            p_volume_actor->GetProperty()->SetColor(this->mPointColor[0], this->mPointColor[1], this->mPointColor[2]);
        }
        pRenderer->AddActor(p_volume_actor);

        if(mShowPottsMeshOutlines)
        {
            auto p_bounds = vtkSmartPointer<vtkPolyData>::New();

            for(unsigned idx=0; idx<p_potts_grid->GetNumberOfCells(); idx++)
            {
                auto p_local_threshold = vtkSmartPointer<vtkThreshold>::New();
                p_local_threshold->SetInputData(p_geom_filter->GetOutput());

#if VTK_MAJOR_VERSION < 9
                p_local_threshold->ThresholdBetween(p_element_base_ids->GetTuple1(idx),
                                                    p_element_base_ids->GetTuple1(idx));
#else
                p_local_threshold->SetLowerThreshold(p_element_base_ids->GetTuple1(idx));
                p_local_threshold->SetUpperThreshold(p_element_base_ids->GetTuple1(idx));
#endif

                p_local_threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Cell Base Id");

                auto p_local_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();
                p_local_geom_filter->SetInputConnection(p_local_threshold->GetOutputPort());

                auto p_features = vtkSmartPointer<vtkFeatureEdges>::New();
                p_features->SetInputConnection(p_local_geom_filter->GetOutputPort());
                p_features->Update();

                auto p_append = vtkSmartPointer<vtkAppendPolyData>::New();
                p_append->AddInputData(p_bounds);
                p_append->AddInputData(p_features->GetOutput());
                p_append->Update();
                p_bounds = p_append->GetOutput();
            }

            auto p_mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
            p_mapper2->SetInputData(p_bounds);

            auto p_volume_actor2 = vtkSmartPointer<vtkActor>::New();
            p_volume_actor2->SetMapper(p_mapper2);
            p_volume_actor2->GetProperty()->SetEdgeVisibility(this->mShowEdges);
            p_volume_actor2->GetProperty()->SetLineWidth(this->mEdgeSize);
            p_volume_actor2->GetProperty()->SetColor(this->mEdgeColor[0], this->mEdgeColor[1], this->mEdgeColor[2]);
            p_volume_actor2->GetProperty()->SetOpacity(this->mVolumeOpacity);
            pRenderer->AddActor(p_volume_actor2);
        }
    }
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::AddActor(vtkSmartPointer<vtkRenderer> pRenderer)
{
    if(!mpCellPopulation)
    {
        return;
    }

    // Show cell centres if requested
    if(mShowCellCentres || boost::dynamic_pointer_cast<CaBasedCellPopulation<DIM> >(mpCellPopulation) ||
            boost::dynamic_pointer_cast<NodeBasedCellPopulation<DIM> >(mpCellPopulation))
    {
        auto p_points = vtkSmartPointer<vtkPoints>::New();
        auto p_cell_color_reference_data = vtkSmartPointer<vtkDoubleArray>::New();
        p_cell_color_reference_data->SetName("CellColors");
        auto p_polydata = vtkSmartPointer<vtkPolyData>::New();

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
             cell_iter != mpCellPopulation->End(); ++cell_iter)
        {
            c_vector<double, DIM> centre = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            if(DIM==3)
            {
                p_points->InsertNextPoint(centre[0], centre[1], centre[2]);
            }
            else
            {
                p_points->InsertNextPoint(centre[0], centre[1], 0.0);
            }

            if(mColorByCellType)
            {
                p_cell_color_reference_data->InsertNextTuple1((*cell_iter)->GetCellProliferativeType()->GetColour());
            }
            else if(mColorByCellData && !this->mDataLabel.empty())
            {
                std::vector<std::string> keys = (*cell_iter)->GetCellData()->GetKeys();
                if (std::find(keys.begin(), keys.end(), this->mDataLabel) != keys.end())
                {
                    p_cell_color_reference_data->InsertNextTuple1((*cell_iter)->GetCellData()->GetItem(this->mDataLabel));
                }
                else
                {
                    p_cell_color_reference_data->InsertNextTuple1(0.0);
                }
            }
            else if(mColorByCellMutationState)
            {
                double mutation_state = (*cell_iter)->GetMutationState()->GetColour();

                CellPropertyCollection collection = (*cell_iter)->rGetCellPropertyCollection();
                CellPropertyCollection label_collection = collection.GetProperties<CellLabel>();

                if (label_collection.GetSize() == 1)
                {
                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(label_collection.GetProperty());
                    mutation_state = p_label->GetColour();
                }
                p_cell_color_reference_data->InsertNextTuple1(mutation_state);
            }
            else if(mColorByCellLabel)
            {
                double label = 0.0;
                if ((*cell_iter)->template HasCellProperty<CellLabel>())
                {
                    CellPropertyCollection collection = (*cell_iter)->rGetCellPropertyCollection().template GetProperties<CellLabel>();
                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
                    label = p_label->GetColour();
                }
                p_cell_color_reference_data->InsertNextTuple1(label);
            }
            else
            {
                p_cell_color_reference_data->InsertNextTuple1((*cell_iter)->GetCellId());
            }
        }
        p_polydata->SetPoints(p_points);
        p_polydata->GetPointData()->AddArray(p_cell_color_reference_data);

        auto p_scaled_ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
        if(!mColorByCellData)
        {
            double range[2];
            p_polydata->GetPointData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpDiscreteColorTransferFunction->GetColor((255.0-double(idx))/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }
        else
        {
            double range[2];
            p_polydata->GetPointData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpColorTransferFunction->GetColor(double(idx)/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }

        auto p_spheres = vtkSmartPointer<vtkSphereSource>::New();
        p_spheres->SetRadius(this->mPointSize);
        p_spheres->SetPhiResolution(16);
        p_spheres->SetThetaResolution(16);

        auto p_glyph = vtkSmartPointer<vtkGlyph3D>::New();
        p_glyph->SetInputData(p_polydata);
        p_glyph->SetSourceConnection(p_spheres->GetOutputPort());
        p_glyph->ClampingOff();
        p_glyph->SetScaleModeToScaleByScalar();
        p_glyph->SetScaleFactor(1.0);
        p_glyph->Update();

        auto p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_mapper->SetInputData(p_glyph->GetOutput());
        p_mapper->SetLookupTable(p_scaled_ctf);
        p_mapper->ScalarVisibilityOn();
        p_mapper->SelectColorArray("CellColors");
        p_mapper->SetScalarModeToUsePointFieldData();
        p_mapper->SetColorModeToMapScalars();

        auto p_actor = vtkSmartPointer<vtkActor>::New();
        p_actor->SetMapper(p_mapper);
        p_actor->GetProperty()->SetOpacity(this->mVolumeOpacity);
        if(mColorCellByUserDefined)
        {
            p_actor->GetProperty()->SetColor(this->mPointColor[0], this->mPointColor[1], this->mPointColor[2]);
        }
        pRenderer->AddActor(p_actor);

        if(!this->mDataLabel.empty() && this->mShowScaleBar)
        {
            this->mpScaleBar->SetLookupTable(p_scaled_ctf);
            this->mpScaleBar->SetTitle(this->mDataLabel.c_str());
            pRenderer->AddActor(this->mpScaleBar);
        }
    }

    if(boost::dynamic_pointer_cast<MeshBasedCellPopulation<DIM> >(mpCellPopulation) && (mShowMutableMeshEdges || mShowVoronoiMeshEdges))
    {
        AddMeshBasedCellPopulationActor(pRenderer);
    }
    else if (boost::dynamic_pointer_cast<VertexBasedCellPopulation<DIM> >(mpCellPopulation) && mShowVoronoiMeshEdges)
    {
        AddVertexBasedCellPopulationActor(pRenderer);
    }
    else if (boost::dynamic_pointer_cast<PottsBasedCellPopulation<DIM> >(mpCellPopulation) && (mShowPottsMeshEdges || mShowPottsMeshOutlines))
    {
        AddPottsBasedCellPopulationActor(pRenderer);
    }
    else if (boost::dynamic_pointer_cast<CaBasedCellPopulation<DIM> >(mpCellPopulation) && mShowPottsMeshEdges)
    {
        AddCaBasedCellPopulationActor(pRenderer);
    }
    else if (boost::dynamic_pointer_cast<ImmersedBoundaryCellPopulation<DIM> >(mpCellPopulation))
    {
        AddImmersedBoundaryCellPopulationActor(pRenderer);
    }
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::AddVertexBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer)
{
    boost::shared_ptr<VertexBasedCellPopulation<DIM> > p_cell_population = boost::dynamic_pointer_cast<VertexBasedCellPopulation<DIM> >(mpCellPopulation);

    if(!p_cell_population)
    {
        EXCEPTION("Could not cast mesh to Vertex Based type.");
    }

    if(mShowVoronoiMeshEdges)
    {
        auto p_voronoi_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        auto p_cell_color_reference_data = vtkSmartPointer<vtkDoubleArray>::New();
        p_cell_color_reference_data->SetName("CellColors");

        auto p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->GetData()->SetName("Vertex positions");
        for (unsigned node_num=0; node_num<p_cell_population->rGetMesh().GetNumNodes(); node_num++)
        {
            c_vector<double, DIM> position = p_cell_population->rGetMesh().GetNode(node_num)->rGetLocation();
            if (DIM==2)
            {
                p_points->InsertPoint(node_num, position[0], position[1], 0.0);
            }
            else
            {
                p_points->InsertPoint(node_num, position[0], position[1], position[2]);
            }
        }
        p_voronoi_grid->SetPoints(p_points);

        for (typename VertexMesh<DIM,DIM>::VertexElementIterator iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
             iter != p_cell_population->rGetMesh().GetElementIteratorEnd(); ++iter)
        {
            vtkSmartPointer<vtkCell> p_cell;
            if (DIM == 2)
            {
                p_cell = vtkSmartPointer<vtkPolygon>::New();
            }
            else
            {
                p_cell = vtkSmartPointer<vtkConvexPointSet>::New();
            }
            vtkSmartPointer<vtkIdList> p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(iter->GetNumNodes());
            for (unsigned j=0; j<iter->GetNumNodes(); ++j)
            {
                p_cell_id_list->SetId(j, iter->GetNodeGlobalIndex(j));
            }
            p_voronoi_grid->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);

            unsigned element_index = iter->GetIndex();
            CellPtr p_biological_cell = p_cell_population->GetCellUsingLocationIndex(element_index);

            if(mColorByCellType)
            {
                p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellProliferativeType()->GetColour());
            }

            else if(mColorByCellData && !this->mDataLabel.empty())
            {
                std::vector<std::string> keys = p_biological_cell->GetCellData()->GetKeys();
                if (std::find(keys.begin(), keys.end(), this->mDataLabel) != keys.end())
                {
                    p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellData()->GetItem(this->mDataLabel));
                }
                else
                {
                    p_cell_color_reference_data->InsertNextTuple1(0.0);
                }
            }

            else if(mColorByCellMutationState)
            {
                double mutation_state = p_biological_cell->GetMutationState()->GetColour();

                CellPropertyCollection collection = p_biological_cell->rGetCellPropertyCollection();
                CellPropertyCollection label_collection = collection.GetProperties<CellLabel>();

                if (label_collection.GetSize() == 1)
                {
                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(label_collection.GetProperty());
                    mutation_state = p_label->GetColour();
                }
                p_cell_color_reference_data->InsertNextTuple1(mutation_state);
            }
            else if(mColorByCellLabel)
            {
                double label = 0.0;
                if (p_biological_cell->HasCellProperty<CellLabel>())
                {
                    CellPropertyCollection collection = p_biological_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
                    label = p_label->GetColour();
                }
                p_cell_color_reference_data->InsertNextTuple1(label);
            }
            else
            {
                p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellId());
            }
        }

        p_voronoi_grid->GetCellData()->AddArray(p_cell_color_reference_data);
        p_voronoi_grid->GetCellData()->SetScalars(p_cell_color_reference_data);

        auto p_scaled_ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
        if(!mColorByCellData)
        {
            double range[2];
            p_voronoi_grid->GetCellData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpDiscreteColorTransferFunction->GetColor((255.0-double(idx))/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }
        else
        {
            double range[2];
            p_voronoi_grid->GetCellData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpColorTransferFunction->GetColor(double(idx)/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }
        p_scaled_ctf->Build();

        auto p_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();
        p_geom_filter->SetInputData(p_voronoi_grid);

        auto p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_mapper->SetInputConnection(p_geom_filter->GetOutputPort());
        p_mapper->SetLookupTable(p_scaled_ctf);
        p_mapper->ScalarVisibilityOn();
        p_mapper->SelectColorArray("CellColors");
        p_mapper->SetScalarModeToUseCellData();
        p_mapper->SetColorModeToMapScalars();

        auto p_actor = vtkSmartPointer<vtkActor>::New();
        p_actor->SetMapper(p_mapper);
        p_actor->GetProperty()->SetOpacity(this->mVolumeOpacity);
        if(mColorCellByUserDefined)
        {
            p_actor->GetProperty()->SetColor(this->mPointColor[0], this->mPointColor[1], this->mPointColor[2]);
        }
        pRenderer->AddActor(p_actor);

        auto p_voronoi_extract_edges = vtkSmartPointer<vtkFeatureEdges>::New();
        p_voronoi_extract_edges->SetInputConnection(p_geom_filter->GetOutputPort());
        p_voronoi_extract_edges->SetFeatureEdges(false);
        p_voronoi_extract_edges->SetBoundaryEdges(true);
        p_voronoi_extract_edges->SetManifoldEdges(true);
        p_voronoi_extract_edges->SetNonManifoldEdges(false);

        auto p_voronoi_tubes = vtkSmartPointer<vtkTubeFilter>::New();
        p_voronoi_tubes->SetInputConnection(p_voronoi_extract_edges->GetOutputPort());
        p_voronoi_tubes->SetRadius(this->mEdgeSize);
        p_voronoi_tubes->SetNumberOfSides(12);

        auto p_voronoi_tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_voronoi_tube_mapper->SetInputConnection(p_voronoi_tubes->GetOutputPort());
        p_voronoi_tube_mapper->ScalarVisibilityOff();

        auto p_voronoi_tube_actor = vtkSmartPointer<vtkActor>::New();
        p_voronoi_tube_actor->SetMapper(p_voronoi_tube_mapper);
        p_voronoi_tube_actor->GetProperty()->SetColor(this->mEdgeColor[0], this->mEdgeColor[1], this->mEdgeColor[2]);
        pRenderer->AddActor(p_voronoi_tube_actor);

        if(!this->mDataLabel.empty() && this->mShowScaleBar)
        {
            this->mpScaleBar->SetLookupTable(p_scaled_ctf);
            this->mpScaleBar->SetTitle(this->mDataLabel.c_str());
            pRenderer->AddActor(this->mpScaleBar);
        }
    }
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::AddImmersedBoundaryCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer)
{
    boost::shared_ptr<ImmersedBoundaryCellPopulation<DIM> > p_cell_population = boost::dynamic_pointer_cast<ImmersedBoundaryCellPopulation<DIM> >(mpCellPopulation);

    if(!p_cell_population)
    {
        EXCEPTION("Could not cast mesh to Immersed Boundary type.");
    }

    if(mShowVoronoiMeshEdges)
    {
        auto p_voronoi_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        auto p_cell_color_reference_data = vtkSmartPointer<vtkDoubleArray>::New();
        p_cell_color_reference_data->SetName("CellColors");

        auto p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->GetData()->SetName("Vertex positions");
        for (unsigned node_num=0; node_num<p_cell_population->rGetMesh().GetNumNodes(); node_num++)
        {
            c_vector<double, DIM> position = p_cell_population->rGetMesh().GetNode(node_num)->rGetLocation();
            if (DIM==2)
            {
                p_points->InsertPoint(node_num, position[0], position[1], 0.0);
            }
            else
            {
                p_points->InsertPoint(node_num, position[0], position[1], position[2]);
            }
        }
        p_voronoi_grid->SetPoints(p_points);

        for (typename ImmersedBoundaryMesh<DIM,DIM>::ImmersedBoundaryElementIterator iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
             iter != p_cell_population->rGetMesh().GetElementIteratorEnd(); ++iter)
        {
            vtkSmartPointer<vtkCell> p_cell;
            if (DIM == 2)
            {
                p_cell = vtkSmartPointer<vtkPolygon>::New();
            }
            else
            {
                p_cell = vtkSmartPointer<vtkConvexPointSet>::New();
            }
            vtkSmartPointer<vtkIdList> p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(iter->GetNumNodes());
            for (unsigned j=0; j<iter->GetNumNodes(); ++j)
            {
                p_cell_id_list->SetId(j, iter->GetNodeGlobalIndex(j));
            }
            p_voronoi_grid->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);

            unsigned element_index = iter->GetIndex();
            CellPtr p_biological_cell = p_cell_population->GetCellUsingLocationIndex(element_index);

            if(mColorByCellType)
            {
                p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellProliferativeType()->GetColour());
            }

            else if(mColorByCellData && !this->mDataLabel.empty())
            {
                std::vector<std::string> keys = p_biological_cell->GetCellData()->GetKeys();
                if (std::find(keys.begin(), keys.end(), this->mDataLabel) != keys.end())
                {
                    p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellData()->GetItem(this->mDataLabel));
                }
                else
                {
                    p_cell_color_reference_data->InsertNextTuple1(0.0);
                }
            }

            else if(mColorByCellMutationState)
            {
                double mutation_state = p_biological_cell->GetMutationState()->GetColour();

                CellPropertyCollection collection = p_biological_cell->rGetCellPropertyCollection();
                CellPropertyCollection label_collection = collection.GetProperties<CellLabel>();

                if (label_collection.GetSize() == 1)
                {
                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(label_collection.GetProperty());
                    mutation_state = p_label->GetColour();
                }
                p_cell_color_reference_data->InsertNextTuple1(mutation_state);
            }
            else if(mColorByCellLabel)
            {
                double label = 0.0;
                if (p_biological_cell->HasCellProperty<CellLabel>())
                {
                    CellPropertyCollection collection = p_biological_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
                    boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
                    label = p_label->GetColour();
                }
                p_cell_color_reference_data->InsertNextTuple1(label);
            }
            else
            {
                p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellId());
            }
        }

        p_voronoi_grid->GetCellData()->AddArray(p_cell_color_reference_data);
        p_voronoi_grid->GetCellData()->SetScalars(p_cell_color_reference_data);

        auto p_scaled_ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
        if(!mColorByCellData)
        {
            double range[2];
            p_voronoi_grid->GetCellData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpDiscreteColorTransferFunction->GetColor((255.0-double(idx))/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }
        else
        {
            double range[2];
            p_voronoi_grid->GetCellData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpColorTransferFunction->GetColor(double(idx)/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }
        p_scaled_ctf->Build();

        auto p_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();
        p_geom_filter->SetInputData(p_voronoi_grid);

        auto p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_mapper->SetInputConnection(p_geom_filter->GetOutputPort());
        p_mapper->SetLookupTable(p_scaled_ctf);
        p_mapper->ScalarVisibilityOn();
        p_mapper->SelectColorArray("CellColors");
        p_mapper->SetScalarModeToUseCellData();
        p_mapper->SetColorModeToMapScalars();

        auto p_actor = vtkSmartPointer<vtkActor>::New();
        p_actor->SetMapper(p_mapper);
        p_actor->GetProperty()->SetOpacity(this->mVolumeOpacity);
        if(mColorCellByUserDefined)
        {
            p_actor->GetProperty()->SetColor(this->mPointColor[0], this->mPointColor[1], this->mPointColor[2]);
        }
        pRenderer->AddActor(p_actor);

        auto p_voronoi_extract_edges = vtkSmartPointer<vtkFeatureEdges>::New();
        p_voronoi_extract_edges->SetInputConnection(p_geom_filter->GetOutputPort());
        p_voronoi_extract_edges->SetFeatureEdges(false);
        p_voronoi_extract_edges->SetBoundaryEdges(true);
        p_voronoi_extract_edges->SetManifoldEdges(true);
        p_voronoi_extract_edges->SetNonManifoldEdges(false);

        auto p_voronoi_tubes = vtkSmartPointer<vtkTubeFilter>::New();
        p_voronoi_tubes->SetInputConnection(p_voronoi_extract_edges->GetOutputPort());
        p_voronoi_tubes->SetRadius(this->mEdgeSize);
        p_voronoi_tubes->SetNumberOfSides(12);

        auto p_voronoi_tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_voronoi_tube_mapper->SetInputConnection(p_voronoi_tubes->GetOutputPort());
        p_voronoi_tube_mapper->ScalarVisibilityOff();

        auto p_voronoi_tube_actor = vtkSmartPointer<vtkActor>::New();
        p_voronoi_tube_actor->SetMapper(p_voronoi_tube_mapper);
        p_voronoi_tube_actor->GetProperty()->SetColor(this->mEdgeColor[0], this->mEdgeColor[1], this->mEdgeColor[2]);
        pRenderer->AddActor(p_voronoi_tube_actor);

        if(!this->mDataLabel.empty() && this->mShowScaleBar)
        {
            this->mpScaleBar->SetLookupTable(p_scaled_ctf);
            this->mpScaleBar->SetTitle(this->mDataLabel.c_str());
            pRenderer->AddActor(this->mpScaleBar);
        }
    }
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::AddMeshBasedCellPopulationActor(vtkSmartPointer<vtkRenderer> pRenderer)
{
    boost::shared_ptr<MeshBasedCellPopulation<DIM> > p_cell_population = boost::dynamic_pointer_cast<MeshBasedCellPopulation<DIM> >(mpCellPopulation);
    boost::shared_ptr<MeshBasedCellPopulationWithGhostNodes<DIM> > p_cell_population_with_ghost =
            boost::dynamic_pointer_cast<MeshBasedCellPopulationWithGhostNodes<DIM> >(mpCellPopulation);


    if(!p_cell_population)
    {
        EXCEPTION("Could not cast mesh to MeshBased type.");
    }

    // Add the voronoi mesh
    if(mShowVoronoiMeshEdges)
    {
        p_cell_population->CreateVoronoiTessellation();
        auto p_voronoi_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        auto p_cell_color_reference_data = vtkSmartPointer<vtkDoubleArray>::New();
        p_cell_color_reference_data->SetName("CellColors");

        if(p_cell_population->GetVoronoiTessellation() != NULL)
        {
            auto p_points = vtkSmartPointer<vtkPoints>::New();
            p_points->GetData()->SetName("Vertex positions");
            for (unsigned node_num=0; node_num<p_cell_population->GetVoronoiTessellation()->GetNumNodes(); node_num++)
            {
                c_vector<double, DIM> position = p_cell_population->GetVoronoiTessellation()->GetNode(node_num)->rGetLocation();
                if (DIM==2)
                {
                    p_points->InsertPoint(node_num, position[0], position[1], 0.0);
                }
                else
                {
                    p_points->InsertPoint(node_num, position[0], position[1], position[2]);
                }
            }
            p_voronoi_grid->SetPoints(p_points);

            for (typename VertexMesh<DIM,DIM>::VertexElementIterator iter = p_cell_population->GetVoronoiTessellation()->GetElementIteratorBegin();
                 iter != p_cell_population->GetVoronoiTessellation()->GetElementIteratorEnd(); ++iter)
            {
                vtkSmartPointer<vtkCell> p_cell;
                if (DIM == 2)
                {
                    p_cell = vtkSmartPointer<vtkPolygon>::New();
                }
                else
                {
                    p_cell = vtkSmartPointer<vtkConvexPointSet>::New();
                }
                vtkSmartPointer<vtkIdList> p_cell_id_list = p_cell->GetPointIds();
                p_cell_id_list->SetNumberOfIds(iter->GetNumNodes());
                for (unsigned j=0; j<iter->GetNumNodes(); ++j)
                {
                    p_cell_id_list->SetId(j, iter->GetNodeGlobalIndex(j));
                }
                p_voronoi_grid->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);

                unsigned node_index =
                        p_cell_population->GetVoronoiTessellation()->GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(iter->GetIndex());

                bool is_ghost_node = false;
                if(p_cell_population_with_ghost)
                {
                    is_ghost_node = p_cell_population_with_ghost->IsGhostNode(node_index);
                }

                if(!is_ghost_node)
                {
                    CellPtr p_biological_cell = p_cell_population->GetCellUsingLocationIndex(node_index);

                    if(mColorByCellType)
                    {
                        p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellProliferativeType()->GetColour());
                    }

                    else if(mColorByCellData && !this->mDataLabel.empty())
                    {
                        std::vector<std::string> keys = p_biological_cell->GetCellData()->GetKeys();
                        if (std::find(keys.begin(), keys.end(), this->mDataLabel) != keys.end())
                        {
                            p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellData()->GetItem(this->mDataLabel));
                        }
                        else
                        {
                            p_cell_color_reference_data->InsertNextTuple1(0.0);
                        }
                    }
                    else if(mColorByCellMutationState)
                    {
                        double mutation_state = p_biological_cell->GetMutationState()->GetColour();

                        CellPropertyCollection collection = p_biological_cell->rGetCellPropertyCollection();
                        CellPropertyCollection label_collection = collection.GetProperties<CellLabel>();

                        if (label_collection.GetSize() == 1)
                        {
                            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(label_collection.GetProperty());
                            mutation_state = p_label->GetColour();
                        }
                        p_cell_color_reference_data->InsertNextTuple1(mutation_state);
                    }
                    else if(mColorByCellLabel)
                    {
                        double label = 0.0;
                        if (p_biological_cell->HasCellProperty<CellLabel>())
                        {
                            CellPropertyCollection collection = p_biological_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
                            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
                            label = p_label->GetColour();
                        }
                        p_cell_color_reference_data->InsertNextTuple1(label);
                    }
                    else
                    {
                        p_cell_color_reference_data->InsertNextTuple1(p_biological_cell->GetCellId());
                    }
                }
                else
                {
                    p_cell_color_reference_data->InsertNextTuple1(-1.0);
                }
            }
        }

        p_voronoi_grid->GetCellData()->AddArray(p_cell_color_reference_data);
        p_voronoi_grid->GetCellData()->SetScalars(p_cell_color_reference_data);

        auto p_scaled_ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
        if(!mColorByCellData)
        {
            double range[2];
            p_voronoi_grid->GetCellData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpDiscreteColorTransferFunction->GetColor((255.0-double(idx))/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }
        else
        {
            double range[2];
            p_voronoi_grid->GetCellData()->GetArray("CellColors")->GetRange(range);
            for(unsigned idx=0; idx<255; idx++)
            {
                double color[3];
                this->mpColorTransferFunction->GetColor(double(idx)/255.0, color);
                p_scaled_ctf->AddRGBPoint(range[0] + double(idx)*(range[1]-range[0])/255.0, color[0], color[1], color[2]);
            }
        }
        p_scaled_ctf->Build();

        auto p_threshold = vtkSmartPointer<vtkThreshold>::New();
        p_threshold->SetInputData(p_voronoi_grid);
#if VTK_MAJOR_VERSION < 9
        p_threshold->ThresholdByUpper(0.0);
#else
        p_threshold->SetUpperThreshold(0.0);
#endif
        p_threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "CellColors");

        auto p_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();
        p_geom_filter->SetInputConnection(p_threshold->GetOutputPort());

        auto p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_mapper->SetInputConnection(p_geom_filter->GetOutputPort());
        p_mapper->SetLookupTable(p_scaled_ctf);
        p_mapper->ScalarVisibilityOn();
        //p_grid_mapper->InterpolateScalarsBeforeMappingOn();
        p_mapper->SelectColorArray("CellColors");
        p_mapper->SetScalarModeToUseCellData();
        p_mapper->SetColorModeToMapScalars();

        auto p_actor = vtkSmartPointer<vtkActor>::New();
        p_actor->SetMapper(p_mapper);
        p_actor->GetProperty()->SetOpacity(this->mVolumeOpacity);
        if(mColorCellByUserDefined)
        {
            p_actor->GetProperty()->SetColor(this->mPointColor[0], this->mPointColor[1], this->mPointColor[2]);
        }
        pRenderer->AddActor(p_actor);

        auto p_voronoi_extract_edges = vtkSmartPointer<vtkFeatureEdges>::New();
        p_voronoi_extract_edges->SetInputConnection(p_geom_filter->GetOutputPort());
        p_voronoi_extract_edges->SetFeatureEdges(false);
        p_voronoi_extract_edges->SetBoundaryEdges(true);
        p_voronoi_extract_edges->SetManifoldEdges(true);
        p_voronoi_extract_edges->SetNonManifoldEdges(false);

        auto p_voronoi_tubes = vtkSmartPointer<vtkTubeFilter>::New();
        p_voronoi_tubes->SetInputConnection(p_voronoi_extract_edges->GetOutputPort());
        p_voronoi_tubes->SetRadius(this->mEdgeSize);
        p_voronoi_tubes->SetNumberOfSides(12);

        auto p_voronoi_tube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_voronoi_tube_mapper->SetInputConnection(p_voronoi_tubes->GetOutputPort());
        p_voronoi_tube_mapper->ScalarVisibilityOff();

        auto p_voronoi_tube_actor = vtkSmartPointer<vtkActor>::New();
        p_voronoi_tube_actor->SetMapper(p_voronoi_tube_mapper);
        p_voronoi_tube_actor->GetProperty()->SetColor(this->mEdgeColor[0], this->mEdgeColor[1], this->mEdgeColor[2]);
        pRenderer->AddActor(p_voronoi_tube_actor);

        if(!this->mDataLabel.empty() && this->mShowScaleBar)
        {
            this->mpScaleBar->SetLookupTable(p_scaled_ctf);
            this->mpScaleBar->SetTitle(this->mDataLabel.c_str());
            pRenderer->AddActor(this->mpScaleBar);
        }
    }

    if(mShowMutableMeshEdges)
    {
        // Do the mutable mesh
        //Make the local mesh into a VtkMesh
        auto p_mutable_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        auto p_points = vtkSmartPointer<vtkPoints>::New();
        p_points->GetData()->SetName("Vertex positions");

        for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = p_cell_population->rGetMesh().GetNodeIteratorBegin();
             node_iter != p_cell_population->rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            c_vector<double, DIM> current_item = node_iter->rGetLocation();
            if (DIM == 3)
            {
                p_points->InsertNextPoint(current_item[0], current_item[1], current_item[2]);
            }
            else if (DIM == 2)
            {
                p_points->InsertNextPoint(current_item[0], current_item[1], 0.0);
            }
            else // (DIM == 1)
            {
                p_points->InsertNextPoint(current_item[0], 0.0, 0.0);
            }
        }

        p_mutable_grid->SetPoints(p_points);

        for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
             elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
             ++elem_iter)
        {

            vtkSmartPointer<vtkCell> p_cell;
            if (DIM == 3)
            {
                p_cell = vtkSmartPointer<vtkTetra>::New();
            }
            else if (DIM == 2)
            {
                p_cell = vtkSmartPointer<vtkTriangle>::New();
            }
            else //(DIM == 1)
            {
                p_cell = vtkSmartPointer<vtkLine>::New();
            }
            vtkSmartPointer<vtkIdList> p_cell_id_list = p_cell->GetPointIds();
            for (unsigned j = 0; j < DIM+1; ++j)
            {
                unsigned global_node_index = elem_iter->GetNodeGlobalIndex(j);
                p_cell_id_list->SetId(j, global_node_index);
            }
            p_mutable_grid->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
        }

        auto p_mutable_geom_filter = vtkSmartPointer<vtkGeometryFilter>::New();
        p_mutable_geom_filter->SetInputData(p_mutable_grid);

        auto p_extract_edges = vtkSmartPointer<vtkFeatureEdges>::New();
        p_extract_edges->SetInputConnection(p_mutable_geom_filter->GetOutputPort());
        p_extract_edges->SetFeatureEdges(false);
        p_extract_edges->SetBoundaryEdges(true);
        p_extract_edges->SetManifoldEdges(true);
        p_extract_edges->SetNonManifoldEdges(false);

        auto p_mutable_tubes = vtkSmartPointer<vtkTubeFilter>::New();
        p_mutable_tubes->SetInputConnection(p_extract_edges->GetOutputPort());
        p_mutable_tubes->SetRadius(0.02);
        p_mutable_tubes->SetNumberOfSides(12);

        auto p_mutable_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        p_mutable_mapper->SetInputConnection(p_mutable_tubes->GetOutputPort());

        auto p_mutable_actor = vtkSmartPointer<vtkActor>::New();
        p_mutable_actor->SetMapper(p_mutable_mapper);
        p_mutable_actor->GetProperty()->SetColor(1,1,1);
        pRenderer->AddActor(p_mutable_actor);
    }
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation)
{
    this->mpCellPopulation = pCellPopulation;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetShowMutableMeshEdges(bool showEdges)
{
    mShowMutableMeshEdges = showEdges;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetShowVoronoiMeshEdges(bool showEdges)
{
    mShowVoronoiMeshEdges = showEdges;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetColorByUserDefined(bool colorByCellUserDefined)
{
    mColorCellByUserDefined = colorByCellUserDefined;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetShowPottsMeshEdges(bool showEdges)
{
    mShowPottsMeshEdges = showEdges;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetColorByCellMutationState(bool colorByCellMutationState)
{
    mColorByCellMutationState = colorByCellMutationState;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetColorByCellLabel(bool colorByCellLabel)
{
    mColorByCellLabel = colorByCellLabel;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetShowPottsMeshOutlines(bool showEdges)
{
    mShowPottsMeshOutlines = showEdges;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetColorByCellType(bool colorByCellType)
{
    mColorByCellType = colorByCellType;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetColorByCellData(bool colorByCellData)
{
    mColorByCellData = colorByCellData;
}

template<unsigned DIM>
void CellPopulationPyChasteActorGenerator<DIM>::SetShowCellCentres(bool showCentres)
{
    mShowCellCentres = showCentres;
}

template class CellPopulationPyChasteActorGenerator<2>;
template class CellPopulationPyChasteActorGenerator<3>;
