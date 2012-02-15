/*

Copyright (C) University of Oxford, 2005-2012

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

#include "CellBasedPdeHandler.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include "CellBasedPdeSolver.hpp"
#include "CellwiseData.hpp"
#include "Exception.hpp"

template<unsigned DIM>
CellBasedPdeHandler<DIM>::CellBasedPdeHandler(AbstractCellPopulation<DIM>* pCellPopulation,
                                              bool deleteMemberPointersInDestructor)
    : mpCellPopulation(pCellPopulation),
      mWriteAverageRadialPdeSolution(false),
      mWriteDailyAverageRadialPdeSolution(false),
      mSetBcsOnCoarseBoundary(true),
      mNumRadialIntervals(UNSIGNED_UNSET),
      mpCoarsePdeMesh(NULL),
      mDeleteMemberPointersInDestructor(deleteMemberPointersInDestructor)
{
    // We must be using a NodeBasedCellPopulation or MeshBasedCellPopulation, with at least one cell
    ///\todo change to exceptions (#1891)
    assert(IsNodeBasedCellPopulation() || IsMeshBasedCellPopulation() || IsPottsBasedCellPopulation() );
    assert(!IsMeshBasedCellPopulationWithGhostNodes());
    assert(mpCellPopulation->GetNumRealCells() != 0);
}

template<unsigned DIM>
CellBasedPdeHandler<DIM>::~CellBasedPdeHandler()
{
    /*
     * Avoid memory leaks. Note that we do not take responsibility for
     * deleting mpCellPopulation, as this object is usually owned by a
     * subclass of AbstractCellBasedSimulation, which deletes the cell
     * population upon destruction if restored from an archive.
     */
    if (mDeleteMemberPointersInDestructor)
    {
        for (unsigned i=0; i<mPdeAndBcCollection.size(); i++)
        {
            delete mPdeAndBcCollection[i];
        }
    }
    if (mpCoarsePdeMesh)
    {
        delete mpCoarsePdeMesh;
    }
}

template<unsigned DIM>
const AbstractCellPopulation<DIM>* CellBasedPdeHandler<DIM>::GetCellPopulation() const
{
    return mpCellPopulation;
}

template<unsigned DIM>
TetrahedralMesh<DIM,DIM>* CellBasedPdeHandler<DIM>::GetCoarsePdeMesh()
{
    return mpCoarsePdeMesh;
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::AddPdeAndBc(PdeAndBoundaryConditions<DIM>* pPdeAndBc)
{
    mPdeAndBcCollection.push_back(pPdeAndBc);
}

template<unsigned DIM>
Vec CellBasedPdeHandler<DIM>::GetPdeSolution(unsigned pdeIndex)
{
    assert(pdeIndex<mPdeAndBcCollection.size());
    return mPdeAndBcCollection[pdeIndex]->GetSolution();
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::InitialiseCellPdeElementMap()
{
    if (mpCoarsePdeMesh == NULL)
    {
        EXCEPTION("InitialiseCellPdeElementMap() should only be called if mpCoarsePdeMesh is set up.");
    }

    mCellPdeElementMap.clear();

    // Find the element of mpCoarsePdeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
         cell_iter != mpCellPopulation->End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = mpCoarsePdeMesh->GetContainingElementIndex(r_position_of_cell);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::UpdateCellPdeElementMap()
{
    // Find the element of mpCoarsePdeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
         cell_iter != mpCellPopulation->End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = mpCoarsePdeMesh->GetContainingElementIndexWithInitialGuess(r_position_of_cell, mCellPdeElementMap[*cell_iter]);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::OpenResultsFiles(std::string outputDirectory)
{
    // If using a NodeBasedCellPopulation or a PottsBasedCellPopulation, mpCoarsePdeMesh must be set up
    if (IsNodeBasedCellPopulation() && mpCoarsePdeMesh==NULL)
    {
        EXCEPTION("Trying to solve a PDE on a NodeBasedCellPopulation without setting up a coarse mesh. Try calling UseCoarsePdeMesh().");
    }
    if (IsPottsBasedCellPopulation() && mpCoarsePdeMesh==NULL)
    {
        EXCEPTION("Trying to solve a PDE on a PottsBasedCellPopulation without setting up a coarse mesh. Try calling UseCoarsePdeMesh().");
    }

    if (mpCoarsePdeMesh != NULL)
    {
        InitialiseCellPdeElementMap();

        // Write mesh to file
        TrianglesMeshWriter<DIM,DIM> mesh_writer(outputDirectory+"/coarse_mesh_output", "coarse_mesh",false);
        mesh_writer.WriteFilesUsingMesh(*mpCoarsePdeMesh);
    }

    if (PetscTools::AmMaster())
    {
        OutputFileHandler output_file_handler(outputDirectory+"/", false);

        if (mpCoarsePdeMesh != NULL)
        {
            mpVizPdeSolutionResultsFile = output_file_handler.OpenOutputFile("results.vizcoarsepdesolution");
        }
        else
        {
            mpVizPdeSolutionResultsFile = output_file_handler.OpenOutputFile("results.vizpdesolution");
        }

        if (mWriteAverageRadialPdeSolution)
        {
            mpAverageRadialPdeSolutionResultsFile = output_file_handler.OpenOutputFile("radial_dist.dat");
        }
    }

    double current_time = SimulationTime::Instance()->GetTime();
    WritePdeSolution(current_time);
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::CloseResultsFiles()
{
    // Close results files
    if (PetscTools::AmMaster())
    {
        mpVizPdeSolutionResultsFile->close();
        if (mWriteAverageRadialPdeSolution)
        {
            WriteAverageRadialPdeSolution(SimulationTime::Instance()->GetTime());
            mpAverageRadialPdeSolutionResultsFile->close();
        }
    }
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::UseCoarsePdeMesh(double stepSize, double meshWidth)
{
    // If solving PDEs on a coarse mesh, each PDE must have an averaged source term
    if (mPdeAndBcCollection.empty())
    {
        EXCEPTION("mPdeAndBcCollection should be populated prior to calling UseCoarsePdeMesh().");
    }
    for (unsigned pde_index=0; pde_index<mPdeAndBcCollection.size(); pde_index++)
    {
        if (mPdeAndBcCollection[pde_index]->HasAveragedSourcePde() == false)
        {
            EXCEPTION("UseCoarsePdeMesh() should only be called if averaged-source PDEs are specified.");
        }
    }

    // Create a regular coarse tetrahedral mesh
    mpCoarsePdeMesh = new TetrahedralMesh<DIM,DIM>;
    switch (DIM)
    {
        case 1:
            mpCoarsePdeMesh->ConstructRegularSlabMesh(stepSize, meshWidth);
            break;
        case 2:
            mpCoarsePdeMesh->ConstructRegularSlabMesh(stepSize, meshWidth, meshWidth);
            break;
        case 3:
            mpCoarsePdeMesh->ConstructRegularSlabMesh(stepSize, meshWidth, meshWidth, meshWidth);
            break;
        default:
            NEVER_REACHED;
    }

    // Find the centre of the coarse PDE mesh
    c_vector<double,DIM> centre_of_coarse_mesh = zero_vector<double>(DIM);
    for (unsigned i=0; i<mpCoarsePdeMesh->GetNumNodes(); i++)
    {
        centre_of_coarse_mesh += mpCoarsePdeMesh->GetNode(i)->rGetLocation();
    }
    centre_of_coarse_mesh /= mpCoarsePdeMesh->GetNumNodes();

    // Translate the centre of coarse PDE mesh to the centre of the cell population
    c_vector<double,DIM> centre_of_cell_population = mpCellPopulation->GetCentroidOfCellPopulation();
    mpCoarsePdeMesh->Translate(centre_of_cell_population - centre_of_coarse_mesh);
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::SolvePdeAndWriteResultsToFile(unsigned samplingTimestepMultiple)
{
    // Record whether we are solving PDEs on a coarse mesh
    bool using_coarse_pde_mesh = (mpCoarsePdeMesh != NULL);

    // If solving PDEs on a coarse mesh, each PDE should have an averaged source term; otherwise none should
    assert(!mPdeAndBcCollection.empty());
    for (unsigned pde_index=0; pde_index<mPdeAndBcCollection.size(); pde_index++)
    {
        assert(mPdeAndBcCollection[pde_index]);
        assert(mPdeAndBcCollection[pde_index]->HasAveragedSourcePde() == using_coarse_pde_mesh);
    }

    // Clear CellwiseData
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    // Make sure the cell population is in a nice state
    mpCellPopulation->Update();

    // Store a pointer to the (population-level or coarse) mesh
    TetrahedralMesh<DIM,DIM>* p_mesh;
    if (using_coarse_pde_mesh)
    {
        p_mesh = mpCoarsePdeMesh;
    }
    else
    {
        // If not using a coarse PDE mesh, we must be using a MeshBasedCellPopulation
        p_mesh = &(static_cast<MeshBasedCellPopulation<DIM>*>(mpCellPopulation)->rGetMesh());
    }

    // Loop over elements of mPdeAndBcCollection
    for (unsigned pde_index=0; pde_index<mPdeAndBcCollection.size(); pde_index++)
    {
        // Get pointer to this PdeAndBoundaryConditions object
        PdeAndBoundaryConditions<DIM>* p_pde_and_bc = mPdeAndBcCollection[pde_index];

        // Set up boundary conditions
        AbstractBoundaryCondition<DIM>* p_bc = p_pde_and_bc->GetBoundaryCondition();
        BoundaryConditionsContainer<DIM,DIM,1> bcc(false);

        if (p_pde_and_bc->IsNeumannBoundaryCondition()) // this BC is of Neumann type
        {
            if (using_coarse_pde_mesh)
            {
                ///\todo enable this (#1891)
                EXCEPTION("Neumann BCs not yet implemented when using a coarse PDE mesh");
            }
            else
            {
                for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = p_mesh->GetBoundaryElementIteratorBegin();
                     elem_iter != p_mesh->GetBoundaryElementIteratorEnd();
                     ++elem_iter)
                {
                    bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc);
                }
            }
        }
        else // assume that if the BC is of Neumann type, then it is Dirichlet
        {
            if (using_coarse_pde_mesh && !mSetBcsOnCoarseBoundary)
            {
                // Get the set of coarse element indices that contain cells
                std::set<unsigned> coarse_element_indices_in_map;
                for (typename AbstractCentreBasedCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
                     cell_iter != mpCellPopulation->End();
                     ++cell_iter)
                {
                    coarse_element_indices_in_map.insert(mCellPdeElementMap[*cell_iter]);
                }

                // Find the node indices associated with elements whose indices are NOT in the set coarse_element_indices_in_map
                std::set<unsigned> coarse_mesh_boundary_node_indices;
                for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
                {
                    if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
                    {
                        Element<DIM,DIM>* p_element = p_mesh->GetElement(i);
                        for (unsigned j=0; j<DIM+1; j++)
                        {
                            unsigned node_index = p_element->GetNodeGlobalIndex(j);
                            coarse_mesh_boundary_node_indices.insert(node_index);
                        }
                    }
                }

                // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
                for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
                     iter != coarse_mesh_boundary_node_indices.end();
                     ++iter)
                {
                    bcc.AddDirichletBoundaryCondition(p_mesh->GetNode(*iter), p_bc, 0, false);
                }
            }
            else // apply BC at boundary nodes of (population-level or coarse) mesh
            {
                for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = p_mesh->GetBoundaryNodeIteratorBegin();
                     node_iter != p_mesh->GetBoundaryNodeIteratorEnd();
                     ++node_iter)
                {
                    bcc.AddDirichletBoundaryCondition(*node_iter, p_bc);
                }
            }
        }

        // If the solution at the previous timestep exists...
        PetscInt previous_solution_size = 0;
        if (p_pde_and_bc->GetSolution())
        {
            VecGetSize(p_pde_and_bc->GetSolution(), &previous_solution_size);
        }

        // ...then record whether it is the correct size...
        bool is_previous_solution_size_correct = (previous_solution_size == (int)p_mesh->GetNumNodes());

        // ...and if it is, store it as an initial guess for the PDE solver
        Vec initial_guess;
        if (is_previous_solution_size_correct)
        {
            // This Vec is copied by the solver's Solve() method, so must be deleted here too
            VecDuplicate(p_pde_and_bc->GetSolution(), &initial_guess);
            VecCopy(p_pde_and_bc->GetSolution(), initial_guess);
            p_pde_and_bc->DestroySolution();
        }
        else
        {
            ///\todo enable the coarse PDE mesh to change size, e.g. for a growing domain (#630/#1891)
            if (!using_coarse_pde_mesh && p_pde_and_bc->GetSolution())
            {
                assert(previous_solution_size != 0);
                p_pde_and_bc->DestroySolution();
            }
        }

        // Create a PDE solver and solve the PDE on the (population-level or coarse) mesh
        if (using_coarse_pde_mesh)
        {
            // When using a coarse PDE mesh, we must set up the source terms before solving the PDE.
            // pass in mCellPdeElementMap to speed up finding cells.
            this->UpdateCellPdeElementMap();
            p_pde_and_bc->SetUpSourceTermsForAveragedSourcePde(p_mesh, &mCellPdeElementMap);

            SimpleLinearEllipticSolver<DIM,DIM> solver(p_mesh, p_pde_and_bc->GetPde(), &bcc);

            // If we have an initial guess, use this when solving the system...
            if (is_previous_solution_size_correct)
            {
                p_pde_and_bc->SetSolution(solver.Solve(initial_guess));
                PetscTools::Destroy(initial_guess);
            }
            else // ...otherwise do not supply one
            {
                p_pde_and_bc->SetSolution(solver.Solve());
            }
        }
        else
        {
            CellBasedPdeSolver<DIM> solver(p_mesh, p_pde_and_bc->GetPde(), &bcc);

            // If we have an initial guess, use this...
            if (is_previous_solution_size_correct)
            {
                p_pde_and_bc->SetSolution(solver.Solve(initial_guess));
                PetscTools::Destroy(initial_guess);
            }
            else // ...otherwise do not supply one
            {
                p_pde_and_bc->SetSolution(solver.Solve());
            }
        }

        // Store the PDE solution in an accessible form
        ReplicatableVector solution_repl(p_pde_and_bc->GetSolution());

        // Having solved the PDE, now update CellwiseData
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
             cell_iter != mpCellPopulation->End();
             ++cell_iter)
        {
            unsigned node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            double solution_at_node = 0.0;

            if (using_coarse_pde_mesh)
            {
                // When using a coarse PDE mesh, the cells are not nodes of the mesh, so we must interpolate

                // Find the element in the coarse mesh that contains this cell. CellElementMap has been updated so use this.
                unsigned elem_index = mCellPdeElementMap[*cell_iter];
                Element<DIM,DIM>* p_element = mpCoarsePdeMesh->GetElement(elem_index);

                const ChastePoint<DIM>& node_location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

                c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(node_location);
                for (unsigned i=0; i<DIM+1; i++)
                {
                    double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(i)];
                    solution_at_node += nodal_value * weights(i);
                }
            }
            else
            {
                solution_at_node = solution_repl[node_index];
            }

            CellwiseData<DIM>::Instance()->SetValue(solution_at_node, node_index, pde_index);
        }
    }

    // Write results to file if required
    SimulationTime* p_time = SimulationTime::Instance();
    double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
    if ((p_time->GetTimeStepsElapsed()+1)%samplingTimestepMultiple == 0)
    {
        WritePdeSolution(time_next_step);
    }

#define COVERAGE_IGNORE
    ///\todo enable this in the case where a coarse PDE mesh is used
    if (!using_coarse_pde_mesh)
    {
        if (mWriteDailyAverageRadialPdeSolution)
        {
            ///\todo Worry about round-off errors (#1891)
            SimulationTime* p_time = SimulationTime::Instance();
            double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
            unsigned num_timesteps_per_day = (unsigned) (DBL_EPSILON + 24/SimulationTime::Instance()->GetTimeStep());
            if ((p_time->GetTimeStepsElapsed()+1) % num_timesteps_per_day == 0)
            {
                WriteAverageRadialPdeSolution(time_next_step);
            }
        }
    }
#undef COVERAGE_IGNORE
}

template<unsigned DIM>
unsigned CellBasedPdeHandler<DIM>::FindCoarseElementContainingCell(CellPtr pCell)
{
    // Get containing element at last timestep from mCellPdeElementMap
    unsigned old_element_index = mCellPdeElementMap[pCell];

    // Create a std::set of guesses for the current containing element
    std::set<unsigned> test_elements;
    test_elements.insert(old_element_index);

    Element<DIM,DIM>* p_element = mpCoarsePdeMesh->GetElement(old_element_index);
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
    const ChastePoint<DIM>& r_cell_position = mpCellPopulation->GetLocationOfCellCentre(pCell);
    unsigned new_element_index = mpCoarsePdeMesh->GetContainingElementIndex(r_cell_position, false, test_elements);

    // Update mCellPdeElementMap
    mCellPdeElementMap[pCell] = new_element_index;

    return new_element_index;
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::WritePdeSolution(double time)
{
    if (PetscTools::AmMaster())
    {
        (*mpVizPdeSolutionResultsFile) << time << "\t";

        for (unsigned pde_index=0; pde_index<mPdeAndBcCollection.size(); pde_index++)
        {
            if (mpCoarsePdeMesh != NULL)
            {
                PdeAndBoundaryConditions<DIM>* p_pde_and_bc = mPdeAndBcCollection[pde_index];

                for (unsigned i=0; i<mpCoarsePdeMesh->GetNumNodes(); i++)
                {
                    (*mpVizPdeSolutionResultsFile) << i << " ";
                    c_vector<double,DIM> location = mpCoarsePdeMesh->GetNode(i)->rGetLocation();
                    for (unsigned k=0; k<DIM; k++)
                    {
                        (*mpVizPdeSolutionResultsFile) << location[k] << " ";
                    }

                    if (p_pde_and_bc->GetSolution())
                    {
                        ReplicatableVector solution_repl(p_pde_and_bc->GetSolution());
                        (*mpVizPdeSolutionResultsFile) << solution_repl[i] << " ";
                    }
                    else
                    {
                        ///\todo consider whether a different initial condition is more appropriate (#1891)

                        // Find the nearest cell to this coarse mesh node
                        unsigned nearest_node_index = 0;
                        double nearest_node_distance = DBL_MAX;
                        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
                             cell_iter != mpCellPopulation->End();
                             ++cell_iter)
                        {
                            c_vector<double, DIM> node_location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
                            if (norm_2(node_location - location) < nearest_node_distance)
                            {
                                nearest_node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                            }
                        }

                        CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(nearest_node_index);
                        double solution = CellwiseData<DIM>::Instance()->GetValue(p_cell, pde_index);
                        (*mpVizPdeSolutionResultsFile) << solution << " ";
                    }
                }
            }
            else
            {
                for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
                     cell_iter != mpCellPopulation->End();
                     ++cell_iter)
                {
                    unsigned node_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                    (*mpVizPdeSolutionResultsFile) << node_index << " ";
                    const c_vector<double,DIM>& position = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
                    for (unsigned i=0; i<DIM; i++)
                    {
                        (*mpVizPdeSolutionResultsFile) << position[i] << " ";
                    }
                    double solution = CellwiseData<DIM>::Instance()->GetValue(*cell_iter, pde_index);
                    (*mpVizPdeSolutionResultsFile) << solution << " ";
                }
            }
        }
        (*mpVizPdeSolutionResultsFile) << "\n";
    }
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::SetWriteAverageRadialPdeSolution(unsigned numRadialIntervals, bool writeDailyResults)
{
    mWriteAverageRadialPdeSolution = true;
    mNumRadialIntervals = numRadialIntervals;
    mWriteDailyAverageRadialPdeSolution = writeDailyResults;
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::SetImposeBcsOnCoarseBoundary(bool setBcsOnCoarseBoundary)
{
    mSetBcsOnCoarseBoundary = setBcsOnCoarseBoundary;
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::WriteAverageRadialPdeSolution(double time)
{
    (*mpAverageRadialPdeSolutionResultsFile) << time << " ";

    // Calculate the centre of the cell population
    c_vector<double,DIM> centre = mpCellPopulation->GetCentroidOfCellPopulation();

    // Calculate the distance between each node and the centre of the cell population, as well as the maximum of these
    std::map<double, CellPtr> radius_cell_map;
    double max_distance_from_centre = 0.0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
         cell_iter != mpCellPopulation->End();
         ++cell_iter)
    {
        double distance = norm_2(mpCellPopulation->GetLocationOfCellCentre(*cell_iter) - centre);
        radius_cell_map[distance] = *cell_iter;

        if (distance > max_distance_from_centre)
        {
            max_distance_from_centre = distance;
        }
    }

    // Create vector of radius intervals
    std::vector<double> radius_intervals;
    for (unsigned i=0; i<mNumRadialIntervals; i++)
    {
        double upper_radius = max_distance_from_centre*((double) i+1)/((double) mNumRadialIntervals);
        radius_intervals.push_back(upper_radius);
    }

    // Calculate PDE solution in each radial interval
    double lower_radius = 0.0;
    for (unsigned i=0; i<mNumRadialIntervals; i++)
    {
        unsigned counter = 0;
        double average_solution = 0.0;

        for (std::map<double, CellPtr>::iterator iter = radius_cell_map.begin(); iter != radius_cell_map.end(); ++iter)
        {
            if (iter->first > lower_radius && iter->first <= radius_intervals[i])
            {
                average_solution += CellwiseData<DIM>::Instance()->GetValue(iter->second);
                counter++;
            }
        }
        if (counter != 0)
        {
            average_solution /= (double) counter;
        }

        // Write results to file
        (*mpAverageRadialPdeSolutionResultsFile) << radius_intervals[i] << " " << average_solution << " ";
        lower_radius = radius_intervals[i];
    }
    (*mpAverageRadialPdeSolutionResultsFile) << "\n";
}

template<unsigned DIM>
bool CellBasedPdeHandler<DIM>::GetWriteAverageRadialPdeSolution()
{
    return mWriteAverageRadialPdeSolution;
}

template<unsigned DIM>
bool CellBasedPdeHandler<DIM>::GetWriteDailyAverageRadialPdeSolution()
{
    return mWriteDailyAverageRadialPdeSolution;
}

template<unsigned DIM>
bool CellBasedPdeHandler<DIM>::GetImposeBcsOnCoarseBoundary()
{
    return mSetBcsOnCoarseBoundary;
}

template<unsigned DIM>
unsigned CellBasedPdeHandler<DIM>::GetNumRadialIntervals()
{
    return mNumRadialIntervals;
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::OutputParameters(out_stream& rParamsFile)
{
    std::string type = GetIdentifier();

    *rParamsFile << "\t\t<" << type << ">\n";
    *rParamsFile << "\t\t<WriteAverageRadialPdeSolution>" << mWriteAverageRadialPdeSolution << "</WriteAverageRadialPdeSolution>\n";
    *rParamsFile << "\t\t<WriteDailyAverageRadialPdeSolution>" << mWriteDailyAverageRadialPdeSolution << "</WriteDailyAverageRadialPdeSolution>\n";
    *rParamsFile << "\t\t<SetBcsOnCoarseBoundary>" << mSetBcsOnCoarseBoundary << "</SetBcsOnCoarseBoundary>\n";
    *rParamsFile << "\t\t<NumRadialIntervals>" << mNumRadialIntervals << "</NumRadialIntervals>\n";
    *rParamsFile << "\t\t</" << type << ">\n";
}

template<unsigned DIM>
bool CellBasedPdeHandler<DIM>::IsNodeBasedCellPopulation()
{
	return (dynamic_cast<NodeBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL);
}

template<unsigned DIM>
bool CellBasedPdeHandler<DIM>::IsMeshBasedCellPopulation()
{
	return (dynamic_cast<MeshBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL);
}

template<unsigned DIM>
bool CellBasedPdeHandler<DIM>::IsPottsBasedCellPopulation()
{
	return (dynamic_cast<PottsBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL);
}

template<unsigned DIM>
bool CellBasedPdeHandler<DIM>::IsMeshBasedCellPopulationWithGhostNodes()
{
	return (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(mpCellPopulation) != NULL);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellBasedPdeHandler)

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellBasedPdeHandler<1>;
template class CellBasedPdeHandler<2>;
template class CellBasedPdeHandler<3>;
