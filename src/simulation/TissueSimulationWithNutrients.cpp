/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#include "TissueSimulationWithNutrients.hpp"
#include "MeshBasedTissueWithGhostNodes.hpp"
#include "SimpleDataWriter.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "SimpleLinearEllipticAssembler.hpp"
#include "TissueSimulationWithNutrientsAssembler.hpp"
#include "CellwiseData.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "TrianglesMeshReader.hpp"

template<unsigned DIM>
TissueSimulationWithNutrients<DIM>::TissueSimulationWithNutrients(AbstractTissue<DIM>& rTissue,
                                   std::vector<AbstractForce<DIM>*> forceCollection,
                                   AbstractLinearEllipticPde<DIM,DIM>* pPde,
                                   AveragedSinksPde<DIM>* pAveragedSinksPde,
                                   bool deleteTissueAndForceCollection,
                                   bool initialiseCells)
    : TissueSimulation<DIM>(rTissue,
                            forceCollection,
                            deleteTissueAndForceCollection,
                            initialiseCells),
      mNutrientSolution(NULL),
      mpPde(pPde),
      mpAveragedSinksPde(pAveragedSinksPde),
      mWriteAverageRadialNutrientResults(false),
      mWriteDailyAverageRadialNutrientResults(false),
      mNumRadialIntervals(0), // 'unset' value
      mpCoarseNutrientMesh(NULL)
{
    // We must be using a mesh-based tissue
    assert(dynamic_cast<MeshBasedTissue<DIM>*>(&(this->mrTissue)) != NULL);

    // We must not have any ghost nodes
    assert(dynamic_cast<MeshBasedTissueWithGhostNodes<DIM>*>(&(this->mrTissue)) == NULL);
}

template<unsigned DIM>
TissueSimulationWithNutrients<DIM>::~TissueSimulationWithNutrients()
{
    if (mNutrientSolution)
    {
        VecDestroy(mNutrientSolution);
    }
    if (mpCoarseNutrientMesh)
    {
        delete mpCoarseNutrientMesh;
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetPde(AbstractLinearEllipticPde<DIM,DIM>* pPde)
{
    mpPde = pPde;
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetAveragedSinksPde(AveragedSinksPde<DIM>* pAveragedSinksPde)
{
    mpAveragedSinksPde = pAveragedSinksPde;
}

template<unsigned DIM>
Vec TissueSimulationWithNutrients<DIM>::GetNutrientSolution()
{
    return mNutrientSolution;
}

//////////////////////////////////////////////////////////////////////////////
//                          Setup/AfterSolve methods                        //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::WriteVisualizerSetupFile()
{
    for (unsigned i=0; i<this->mForceCollection.size(); i++)
    {
        if (dynamic_cast<AbstractTwoBodyInteractionForce<DIM>*>(this->mForceCollection[i]))
        {
            double cutoff = (static_cast<AbstractTwoBodyInteractionForce<DIM>*>(this->mForceCollection[i]))->GetCutoffPoint();
            *(this->mpSetupFile) << "Cutoff\t" << cutoff << "\n";
        }
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetupSolve()
{
    if (mpCoarseNutrientMesh!=NULL)
    {
        InitialiseCoarseNutrientMesh();
    }
    if (this->mrTissue.Begin() != this->mrTissue.End())
    {
        SetupWriteNutrient();
        double current_time = SimulationTime::Instance()->GetTime();
        WriteNutrient(current_time);
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetupWriteNutrient()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/", false);
    if (PetscTools::AmMaster())
    {
        mpNutrientResultsFile = output_file_handler.OpenOutputFile("results.viznutrient");
        *this->mpSetupFile << "Nutrient \n";
        if (mWriteAverageRadialNutrientResults)
        {
            mpAverageRadialNutrientResultsFile = output_file_handler.OpenOutputFile("radial_dist.dat");
        }
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::UseCoarseNutrientMesh(double coarseGrainScaleFactor)
{
    assert(mpAveragedSinksPde);
    CreateCoarseNutrientMesh(coarseGrainScaleFactor);
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::CreateCoarseNutrientMesh(double coarseGrainScaleFactor)
{
    EXCEPTION("This method is only implemented in 2D");
}

/**
 * The CreateCoarseNutrientMesh method is currently only implemented in 2D, hence there
 * are two definitions to this method (one templated and one not).
 *
 * @param coarseGrainScaleFactor the ratio of the width of the coarse nutrient mesh to the initial width of the tissue
 */
template<>
void TissueSimulationWithNutrients<2>::CreateCoarseNutrientMesh(double coarseGrainScaleFactor)
{
    // Create coarse nutrient mesh (can use a larger mesh if required, e.g. disk_984_elements)
    TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
    mpCoarseNutrientMesh = new TetrahedralMesh<2,2>;
    mpCoarseNutrientMesh->ConstructFromMeshReader(mesh_reader);

    // Find centre of tissue
    c_vector<double,2> centre_of_tissue = zero_vector<double>(2);
    for (AbstractTissue<2>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        centre_of_tissue += this->mrTissue.GetLocationOfCellCentre(*cell_iter);
    }
    centre_of_tissue /= this->mrTissue.GetNumRealCells();

    // Find max radius of tissue
    double max_tissue_radius = 0.0;
    for (AbstractTissue<2>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        double radius = norm_2(centre_of_tissue - this->mrTissue.GetLocationOfCellCentre(*cell_iter));
        if (radius > max_tissue_radius)
        {
            max_tissue_radius = radius;
        }
    }

    // Find centre of coarse nutrient mesh
    c_vector<double,2> centre_of_nutrient_mesh = zero_vector<double>(2);

    for (unsigned i=0; i<mpCoarseNutrientMesh->GetNumNodes(); i++)
    {
        centre_of_nutrient_mesh += mpCoarseNutrientMesh->GetNode(i)->rGetLocation();
    }
    centre_of_nutrient_mesh /= mpCoarseNutrientMesh->GetNumNodes();

    // Find max radius of coarse nutrient mesh
    double max_mesh_radius = 0.0;
    for (unsigned i=0; i<mpCoarseNutrientMesh->GetNumNodes(); i++)
    {
        double radius = norm_2(centre_of_nutrient_mesh - mpCoarseNutrientMesh->GetNode(i)->rGetLocation());
        if (radius > max_mesh_radius)
        {
            max_mesh_radius = radius;
        }
    }

    // Translate centre of coarse nutrient mesh to the origin
    mpCoarseNutrientMesh->Translate(-centre_of_nutrient_mesh[0], -centre_of_nutrient_mesh[1]);

    // Scale nutrient mesh
    double scale_factor = (max_tissue_radius/max_mesh_radius)*coarseGrainScaleFactor;
    mpCoarseNutrientMesh->Scale(scale_factor, scale_factor);

    // Translate centre of coarse nutrient mesh to centre of the tissue
    mpCoarseNutrientMesh->Translate(centre_of_tissue[0], centre_of_tissue[1]);
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::InitialiseCoarseNutrientMesh()
{
    mCellNutrientElementMap.clear();

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        // Find the element of mpCoarseNutrientMesh that contains this cell
        const ChastePoint<DIM>& r_position_of_cell = this->mrTissue.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = mpCoarseNutrientMesh->GetContainingElementIndex(r_position_of_cell);
        mCellNutrientElementMap[&(*cell_iter)] = elem_index;
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::AfterSolve()
{
    if (this->mrTissue.Begin() != this->mrTissue.End() // if there are any cells
    && PetscTools::AmMaster())
    {
        mpNutrientResultsFile->close();

        if (mWriteAverageRadialNutrientResults)
        {
            WriteAverageRadialNutrientDistribution(SimulationTime::Instance()->GetTime(), mNumRadialIntervals);
            mpAverageRadialNutrientResultsFile->close();
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
//                             PostSolve methods                            //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SolveNutrientPde()
{
    if (mpCoarseNutrientMesh!=NULL)
    {
        SolveNutrientPdeUsingCoarseMesh();
        return;
    }

    assert(mpAveragedSinksPde == NULL);
    assert(mpPde);

    // Note: If not using a coarse nutrient mesh, we MUST be using a MeshBasedTissue
    // Make sure the mesh is in a nice state
    this->mrTissue.Update();

    TetrahedralMesh<DIM,DIM>& r_mesh = static_cast<MeshBasedTissue<DIM>*>(&(this->mrTissue))->rGetMesh();
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_boundary_condition = new ConstBoundaryCondition<DIM>(1.0);
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = r_mesh.GetBoundaryNodeIteratorBegin();
         node_iter != r_mesh.GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
        bcc.AddDirichletBoundaryCondition(*node_iter, p_boundary_condition);
    }

    /*
     * Set up assembler. This is a purpose-made elliptic assembler which must
     * interpolate contributions to source terms from nodes onto Gauss points,
     * because the nutrient concentration is only stored at the cells (nodes).
     */
    TissueSimulationWithNutrientsAssembler<DIM> assembler(&r_mesh, mpPde, &bcc);

    PetscInt size_of_soln_previous_step = 0;

    if (mNutrientSolution)
    {
        VecGetSize(mNutrientSolution, &size_of_soln_previous_step);
    }
    if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
    {
        // We make an initial guess which gets copied by the Solve method of
        // SimpleLinearSolver, so we need to delete it too.
        Vec initial_guess;
        VecDuplicate(mNutrientSolution, &initial_guess);
        VecCopy(mNutrientSolution, initial_guess);

        // Use current solution as the initial guess
        VecDestroy(mNutrientSolution);    // Solve method makes its own mNutrientSolution
        mNutrientSolution = assembler.Solve(initial_guess);
        VecDestroy(initial_guess);
    }
    else
    {
        if (mNutrientSolution)
        {
            assert(size_of_soln_previous_step != 0);
            VecDestroy(mNutrientSolution);
        }
        mNutrientSolution = assembler.Solve();
    }

    ReplicatableVector result_repl(mNutrientSolution);

    // Update cellwise data
    for (unsigned i=0; i<r_mesh.GetNumNodes(); i++)
    {
        double oxygen_conc = result_repl[i];
        CellwiseData<DIM>::Instance()->SetValue(oxygen_conc, r_mesh.GetNode(i));
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SolveNutrientPdeUsingCoarseMesh()
{
    assert(mpPde==NULL);
    assert(mpAveragedSinksPde);

    TetrahedralMesh<DIM,DIM>& r_mesh = *mpCoarseNutrientMesh;
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // Loop over cells and calculate centre of distribution
    c_vector<double, DIM> centre = zero_vector<double>(DIM);
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        centre += this->mrTissue.GetLocationOfCellCentre(*cell_iter);
    }
    centre /= this->mrTissue.GetNumRealCells();

    // Find max radius
    double max_radius = 0.0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        double radius = norm_2(centre - this->mrTissue.GetLocationOfCellCentre(*cell_iter));
        if (radius > max_radius)
        {
            max_radius = radius;
        }
    }

    // Set up boundary conditions
    BoundaryConditionsContainer<DIM,DIM,1> bcc;
    ConstBoundaryCondition<DIM>* p_boundary_condition = new ConstBoundaryCondition<DIM>(1.0);

    // Get the set of coarse element indices that contain tissue cells
    std::set<unsigned> coarse_element_indices_in_map;
    for (typename MeshBasedTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        coarse_element_indices_in_map.insert(mCellNutrientElementMap[&(*cell_iter)]);
    }

    // Find the node indices that associated with elements whose
    // indices are NOT in the set coarse_element_indices_in_map
    std::set<unsigned> coarse_mesh_boundary_node_indices;

    for (unsigned i=0; i<r_mesh.GetNumElements(); i++)
    {
        // If the element index is NOT in the set...
        if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
        {
            // ... then get the element...
            Element<DIM,DIM>* p_element = r_mesh.GetElement(i);

            // ... and add its associated nodes to coarse_mesh_boundary_node_indices
            for (unsigned local_index=0; local_index<DIM+1; local_index++)
            {
                unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
                coarse_mesh_boundary_node_indices.insert(node_index);
            }
        }
    }

    // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
    for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
         iter != coarse_mesh_boundary_node_indices.end();
         ++iter)
    {
        bcc.AddDirichletBoundaryCondition(r_mesh.GetNode(*iter), p_boundary_condition, 0, false);
    }

    PetscInt size_of_soln_previous_step = 0;

    if (mNutrientSolution)
    {
        VecGetSize(mNutrientSolution, &size_of_soln_previous_step);
    }

    mpAveragedSinksPde->SetupSourceTerms(*mpCoarseNutrientMesh);

    SimpleLinearEllipticAssembler<DIM,DIM> assembler(mpCoarseNutrientMesh, mpAveragedSinksPde, &bcc);

    if (size_of_soln_previous_step == (int)r_mesh.GetNumNodes())
    {
        // We make an initial guess which gets copied by the Solve method of
        // SimpleLinearSolver, so we need to delete it too.
        Vec initial_guess;
        VecDuplicate(mNutrientSolution, &initial_guess);
        VecCopy(mNutrientSolution, initial_guess);

        // Use current solution as the initial guess
        VecDestroy(mNutrientSolution);    // Solve method makes its own mNutrientSolution
        mNutrientSolution = assembler.Solve(initial_guess);
        VecDestroy(initial_guess);
    }
    else
    {
        assert(mNutrientSolution == NULL);
        /**
         * Eventually we will enable the coarse nutrient mesh to change size, for example
         * in the case of a spheroid that grows a lot (see #630). In this case we should
         * uncomment the following code.
         *
        if (mNutrientSolution)
        {
            assert(0);
            VecDestroy(mNutrientSolution);
        }
        *
        */
        mNutrientSolution = assembler.Solve();
    }

    // Update cellwise data - since the cells are not nodes on the coarse
    // mesh, we have to interpolate from the nodes of the coarse mesh onto
    // the cell locations
    ReplicatableVector nutrient_repl(mNutrientSolution);

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
        cell_iter != this->mrTissue.End();
        ++cell_iter)
    {
        // Find coarse mesh element containing cell
        unsigned elem_index = FindElementContainingCell(*cell_iter);

        Element<DIM,DIM>* p_element = mpCoarseNutrientMesh->GetElement(elem_index);

        const ChastePoint<DIM>& r_position_of_cell = this->mrTissue.GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(r_position_of_cell);

        double interpolated_nutrient = 0.0;
        for (unsigned i=0; i<DIM+1/*num_nodes*/; i++)
        {
            double nodal_value = nutrient_repl[ p_element->GetNodeGlobalIndex(i) ];
            interpolated_nutrient += nodal_value*weights(i);
        }

        CellwiseData<DIM>::Instance()->SetValue(interpolated_nutrient, (static_cast<AbstractCellCentreBasedTissue<DIM>*>(&(this->mrTissue)))->GetNodeCorrespondingToCell(*cell_iter));
    }
}

template<unsigned DIM>
unsigned TissueSimulationWithNutrients<DIM>::FindElementContainingCell(TissueCell& rCell)
{
    // Get containing element at last timestep from mCellNutrientElementMap
    unsigned old_element_index = mCellNutrientElementMap[&rCell];

    // Create a std::set of guesses for the current containing element
    std::set<unsigned> test_elements;
    test_elements.insert(old_element_index);

    Element<DIM,DIM>* p_element = mpCoarseNutrientMesh->GetElement(old_element_index);

    for (unsigned local_index=0; local_index<DIM+1; local_index++)
    {
        std::set<unsigned> element_indices = p_element->GetNode(local_index)->rGetContainingElementIndices();

        for (std::set<unsigned>::iterator iter = element_indices.begin();
             iter != element_indices.end();
             ++iter)
        {
            test_elements.insert(*iter);
        }
    }

    // Find new element, using the previous one as a guess
    const ChastePoint<DIM>& r_cell_position = this->mrTissue.GetLocationOfCellCentre(rCell);
    unsigned new_element_index = mpCoarseNutrientMesh->GetContainingElementIndex(r_cell_position, false, test_elements);

    // Update mCellNutrientElementMap
    mCellNutrientElementMap[&rCell] = new_element_index;

    return new_element_index;
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::PostSolve()
{
    SolveNutrientPde();

    // Save results to file
    SimulationTime* p_time = SimulationTime::Instance();

    double time_next_step = p_time->GetTime() + p_time->GetTimeStep();

    if ((p_time->GetTimeStepsElapsed()+1)%this->mSamplingTimestepMultiple == 0)
    {
        WriteNutrient(time_next_step);
    }

#define COVERAGE_IGNORE
    if (mWriteDailyAverageRadialNutrientResults)
    {
        ///\todo Worry about round-off errors
        unsigned num_timesteps_per_day = (unsigned) (DBL_EPSILON + 24/SimulationTime::Instance()->GetTimeStep());

        if ((p_time->GetTimeStepsElapsed()+1) % num_timesteps_per_day == 0)
        {
            WriteAverageRadialNutrientDistribution(time_next_step, mNumRadialIntervals);
        }
    }
#undef COVERAGE_IGNORE

}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::WriteNutrient(double time)
{
    if (PetscTools::AmMaster())
    {
        // Since there are no ghost nodes, the number of nodes must equal the number of real cells
        assert(this->mrTissue.GetNumNodes()==this->mrTissue.GetNumRealCells());

        (*mpNutrientResultsFile) << time << "\t";

        for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
             cell_iter != this->mrTissue.End();
             ++cell_iter)
        {
            unsigned global_index = this->mrTissue.GetLocationIndexUsingCell(*cell_iter);
            (*mpNutrientResultsFile) << global_index << " ";

            const c_vector<double,DIM>& position = this->mrTissue.GetLocationOfCellCentre(*cell_iter);
            for (unsigned i=0; i<DIM; i++)
            {
                (*mpNutrientResultsFile) << position[i] << " ";
            }

            double nutrient = CellwiseData<DIM>::Instance()->GetValue(*cell_iter);
            (*mpNutrientResultsFile) << nutrient << " ";
        }
        (*mpNutrientResultsFile) << "\n";
    }
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::SetWriteAverageRadialNutrientResults(unsigned numRadialIntervals, bool writeDailyResults)
{
    mWriteAverageRadialNutrientResults = true;
    mNumRadialIntervals = numRadialIntervals;
    mWriteDailyAverageRadialNutrientResults = writeDailyResults;
}

template<unsigned DIM>
void TissueSimulationWithNutrients<DIM>::WriteAverageRadialNutrientDistribution(double time, unsigned numRadialIntervals)
{
    (*mpAverageRadialNutrientResultsFile) << time << " ";

    // Calculate the centre of the tissue
    c_vector<double,DIM> centre = zero_vector<double>(DIM);
    double num_nodes_as_double = (double) this->mrTissue.GetNumNodes();

    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
         cell_iter != this->mrTissue.End();
         ++cell_iter)
    {
       centre += (this->mrTissue.GetLocationOfCellCentre(*cell_iter)) / num_nodes_as_double; 
    }

    // Calculate the distance between each node and the centre of the tissue, as well as the maximum of these
    std::map<double, TissueCell*> distance_cell_map;

    double max_distance_from_centre = 0.0;
    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
         cell_iter != this->mrTissue.End();
         ++cell_iter)
    {
        double distance = norm_2(this->mrTissue.GetLocationOfCellCentre(*cell_iter) - centre);
        distance_cell_map[distance] = &(*cell_iter);

        if (distance > max_distance_from_centre)
        {
            max_distance_from_centre = distance;
        }
    }

    // Create vector of radius intervals
    std::vector<double> radius_intervals;
    for (unsigned i=0; i<numRadialIntervals; i++)
    {
        double upper_radius = max_distance_from_centre*((double) i+1)/((double) numRadialIntervals);
        radius_intervals.push_back(upper_radius);
    }

    // Calculate nutrient concentration in each radial interval
    double lower_radius = 0.0;
    for (unsigned i=0; i<numRadialIntervals; i++)
    {
        unsigned counter = 0;
        double average_conc = 0.0;

        for (std::map<double, TissueCell*>::iterator iter = distance_cell_map.begin();
             iter != distance_cell_map.end();
             ++iter)
        {
            if (iter->first > lower_radius && iter->first <= radius_intervals[i])
            {
                average_conc += CellwiseData<DIM>::Instance()->GetValue(*(iter->second));
                counter++;
            }
        }
        if (counter > 0)
        {
            average_conc /= (double) counter;
        }

        // Write results to file
        (*mpAverageRadialNutrientResultsFile) << radius_intervals[i] << " " << average_conc << " ";
        lower_radius = radius_intervals[i];
    }
    (*mpAverageRadialNutrientResultsFile) << "\n";
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////


template class TissueSimulationWithNutrients<1>;
template class TissueSimulationWithNutrients<2>;
template class TissueSimulationWithNutrients<3>;
