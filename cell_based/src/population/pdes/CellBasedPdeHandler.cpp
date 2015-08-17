/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "CellBasedPdeHandler.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include "CellBasedPdeSolver.hpp"
#include "Exception.hpp"
#include "VtkMeshWriter.hpp"

template<unsigned DIM>
CellBasedPdeHandler<DIM>::CellBasedPdeHandler(AbstractCellPopulation<DIM>* pCellPopulation,
                                              bool deleteMemberPointersInDestructor)
    : mpCellPopulation(pCellPopulation),
      mWriteAverageRadialPdeSolution(false),
      mWriteDailyAverageRadialPdeSolution(false),
      mAverageRadialSolutionVariableName(""),
      mSetBcsOnCoarseBoundary(true),
      mNumRadialIntervals(UNSIGNED_UNSET),
      mpCoarsePdeMesh(NULL),
      mDeleteMemberPointersInDestructor(deleteMemberPointersInDestructor)
{
    // We must be using a CellPopulation with at least one cell
    ///\todo change to exceptions (#1891)
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
    if (!mPdeAndBcCollection.empty() && (pPdeAndBc->rGetDependentVariableName()==""))
    {
        EXCEPTION("When adding more than one PDE to CellBasedPdeHandler set the dependent variable name using SetDependentVariableName(name).");
    }
    for (unsigned i=0; i< mPdeAndBcCollection.size(); i++)
    {
        if (pPdeAndBc->rGetDependentVariableName() == mPdeAndBcCollection[i]->rGetDependentVariableName())
        {
            EXCEPTION("The name "+ pPdeAndBc->rGetDependentVariableName()+ " has already been used in the PDE collection");
        }
    }
    mPdeAndBcCollection.push_back(pPdeAndBc);
}

template<unsigned DIM>
Vec CellBasedPdeHandler<DIM>::GetPdeSolution(const std::string& rName)
{
    if (rName != "")
    {
        for (unsigned i=0; i<mPdeAndBcCollection.size(); i++)
        {
            if (mPdeAndBcCollection[i]->rGetDependentVariableName() == rName)
            {
                return mPdeAndBcCollection[i]->GetSolution();
            }
        }

        EXCEPTION("The PDE collection does not contain a PDE named " + rName);
    }
    else
    {
        assert(mPdeAndBcCollection.size()==1);
        return mPdeAndBcCollection[0]->GetSolution();
    }
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
    // If appropriate, make a coarse mesh which exactly overlays the lattice sites of a PottsMesh (used for all OnLattice simulations)
    if ((dynamic_cast<CaBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL) && mpCoarsePdeMesh==NULL)
    {
        assert(DIM ==2);
        ChasteCuboid<DIM> cuboid = mpCellPopulation->rGetMesh().CalculateBoundingBox();

        // Currently only works with square meshes
        assert(cuboid.GetWidth(0) == cuboid.GetWidth(1));

        UseCoarsePdeMesh(1, cuboid, false);
    }

    // If using a NodeBasedCellPopulation a VertexBasedCellPopulation, a CaBasedCellPopulation or a PottsBasedCellPopulation, mpCoarsePdeMesh must be set up
    if (PdeSolveNeedsCoarseMesh() && mpCoarsePdeMesh==NULL)
    {
        EXCEPTION("Trying to solve a PDE on a cell population that doesn't have a mesh. Try calling UseCoarsePdeMesh().");
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

    mDirPath = outputDirectory; // caching the path to the output directory for VTK
    //double current_time = SimulationTime::Instance()->GetTime();
    //WritePdeSolution(current_time);
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
void CellBasedPdeHandler<DIM>::UseCoarsePdeMesh(double stepSize, ChasteCuboid<DIM> meshCuboid, bool centreOnCellPopulation)
{
    if (mPdeAndBcCollection.empty())
    {
        EXCEPTION("mPdeAndBcCollection should be populated prior to calling UseCoarsePdeMesh().");
    }
    // If solving PDEs on a coarse mesh, each PDE must have an averaged source term
    for (unsigned pde_index=0; pde_index<mPdeAndBcCollection.size(); pde_index++)
    {
        if (mPdeAndBcCollection[pde_index]->HasAveragedSourcePde() == false && !dynamic_cast<CaBasedCellPopulation<DIM>*>(mpCellPopulation))
        {
            EXCEPTION("UseCoarsePdeMesh() should only be called if averaged-source PDEs are specified.");
        }
    }

    // Create a regular coarse tetrahedral mesh
    mpCoarsePdeMesh = new TetrahedralMesh<DIM,DIM>;
    switch (DIM)
    {
        case 1:
            mpCoarsePdeMesh->ConstructRegularSlabMesh(stepSize, meshCuboid.GetWidth(0));
            break;
        case 2:
            mpCoarsePdeMesh->ConstructRegularSlabMesh(stepSize, meshCuboid.GetWidth(0), meshCuboid.GetWidth(1));
            break;
        case 3:
            mpCoarsePdeMesh->ConstructRegularSlabMesh(stepSize, meshCuboid.GetWidth(0), meshCuboid.GetWidth(1), meshCuboid.GetWidth(2));
            break;
        default:
            NEVER_REACHED;
    }

    if (centreOnCellPopulation)
    {
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
    else
    {
        // Get centroid of meshCuboid
        ChastePoint<DIM> upper = meshCuboid.rGetUpperCorner();
        ChastePoint<DIM> lower = meshCuboid.rGetLowerCorner();
        c_vector<double,DIM> centre_of_cuboid = 0.5*(upper.rGetLocation() + lower.rGetLocation());

        // Find the centre of the coarse PDE mesh
        c_vector<double,DIM> centre_of_coarse_mesh = zero_vector<double>(DIM);
        for (unsigned i=0; i<mpCoarsePdeMesh->GetNumNodes(); i++)
        {
            centre_of_coarse_mesh += mpCoarsePdeMesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_mesh /= mpCoarsePdeMesh->GetNumNodes();

        mpCoarsePdeMesh->Translate(centre_of_cuboid - centre_of_coarse_mesh);
    }
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
        assert(mPdeAndBcCollection[pde_index]->HasAveragedSourcePde() == using_coarse_pde_mesh || dynamic_cast<CaBasedCellPopulation<DIM>*>(mpCellPopulation));
    }

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
        std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer(p_pde_and_bc, p_mesh);

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

        if (using_coarse_pde_mesh)
        {
            // We update mCellPdeElementMap before setting up source terms to speed up
            // finding cells in the case of AveragedSourcePdes. This is also used with
            // non-AveragedSourcePdes when using a CaBasedCellPopulation.
            this->UpdateCellPdeElementMap();
        }

        // Create a PDE solver and solve the PDE on the (population-level or coarse) mesh
        if (p_pde_and_bc->HasAveragedSourcePde())
        {
            // When using a coarse PDE mesh, we must set up the source terms before solving the PDE.
            // Pass in already updated CellPdeElementMap to speed up finding cells.
            p_pde_and_bc->SetUpSourceTermsForAveragedSourcePde(p_mesh, &mCellPdeElementMap);

            SimpleLinearEllipticSolver<DIM,DIM> solver(p_mesh, p_pde_and_bc->GetPde(), p_bcc.get());

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
            CellBasedPdeSolver<DIM> solver(p_mesh, p_pde_and_bc->GetPde(), p_bcc.get());

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

        // Having solved the PDE, now update CellData
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
            cell_iter->GetCellData()->SetItem(mPdeAndBcCollection[pde_index]->rGetDependentVariableName(), solution_at_node);
        }
    }

    // Write results to file if required
    SimulationTime* p_time = SimulationTime::Instance();
    if ((p_time->GetTimeStepsElapsed())%samplingTimestepMultiple == 0)
    {
        WritePdeSolution(p_time->GetTime());
    }
#define COVERAGE_IGNORE
    ///\todo enable this in the case where a coarse PDE mesh is used
    if (!using_coarse_pde_mesh)
    {
        if (mWriteDailyAverageRadialPdeSolution)
        {
            ///\todo Worry about round-off errors (#1891)
            p_time = SimulationTime::Instance();
            unsigned num_timesteps_per_day = (unsigned) (DBL_EPSILON + 24/SimulationTime::Instance()->GetTimeStep());
            if ((p_time->GetTimeStepsElapsed()) % num_timesteps_per_day == 0)
            {
                WriteAverageRadialPdeSolution(p_time->GetTime());
            }
        }
    }
#undef COVERAGE_IGNORE
}

template<unsigned DIM>
std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > CellBasedPdeHandler<DIM>::ConstructBoundaryConditionsContainer(
        PdeAndBoundaryConditions<DIM>* pPdeAndBc,
        TetrahedralMesh<DIM,DIM>* pMesh)
{
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    AbstractBoundaryCondition<DIM>* p_bc = pPdeAndBc->GetBoundaryCondition();

    if (pPdeAndBc->IsNeumannBoundaryCondition()) // this BC is of Neumann type
    {
        // Note p_mesh is the coarse mesh or the natural mesh as appropriate
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = pMesh->GetBoundaryElementIteratorBegin();
             elem_iter != pMesh->GetBoundaryElementIteratorEnd();
             ++elem_iter)
        {
            p_bcc->AddNeumannBoundaryCondition(*elem_iter, p_bc);
        }
    }
    else // assume that if the BC is not of Neumann type, then it is Dirichlet
    {
        bool using_coarse_pde_mesh = (mpCoarsePdeMesh != NULL);

        if (using_coarse_pde_mesh && !mSetBcsOnCoarseBoundary)
        {
            // Get the set of coarse element indices that contain cells
            std::set<unsigned> coarse_element_indices_in_map;
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
                 cell_iter != mpCellPopulation->End();
                 ++cell_iter)
            {
                coarse_element_indices_in_map.insert(mCellPdeElementMap[*cell_iter]);
            }

            // Find the node indices associated with elements whose indices are NOT in the set coarse_element_indices_in_map
            std::set<unsigned> coarse_mesh_boundary_node_indices;
            for (unsigned i=0; i<pMesh->GetNumElements(); i++)
            {
                if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
                {
                    Element<DIM,DIM>* p_element = pMesh->GetElement(i);
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
                p_bcc->AddDirichletBoundaryCondition(pMesh->GetNode(*iter), p_bc, 0, false);
            }
        }
        else // apply BC at boundary nodes of (population-level or coarse) mesh
        {
            for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = pMesh->GetBoundaryNodeIteratorBegin();
                 node_iter != pMesh->GetBoundaryNodeIteratorEnd();
                 ++node_iter)
            {
                p_bcc->AddDirichletBoundaryCondition(*node_iter, p_bc);
            }
        }
    }

    return p_bcc;
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

#ifdef CHASTE_VTK
        // Note that this mesh writer is only constructed and used if mpCoarsePdeMesh exists
        VtkMeshWriter<DIM,DIM>* p_vtk_mesh_writer = NULL;
        if (DIM>1 && mpCoarsePdeMesh != NULL )
        {
            std::ostringstream time_string;
            time_string << SimulationTime::Instance()->GetTimeStepsElapsed()+1;
            std::string results_file = "pde_results_"+time_string.str();
            // Note that this mesh writer is always constructed, but is only used if mpCoarsePdeMesh exists
            p_vtk_mesh_writer = new VtkMeshWriter<DIM,DIM>(mDirPath, results_file, false);
        }
#endif //CHASTE_VTK
        for (unsigned pde_index=0; pde_index<mPdeAndBcCollection.size(); pde_index++)
        {
            if (mpCoarsePdeMesh != NULL)
            {
                PdeAndBoundaryConditions<DIM>* p_pde_and_bc = mPdeAndBcCollection[pde_index];
                assert(p_pde_and_bc->rGetDependentVariableName() != "");

#ifdef CHASTE_VTK
                if (p_pde_and_bc->GetSolution())
                {
                    if (DIM>1)
                    {
                        ReplicatableVector solution_repl(p_pde_and_bc->GetSolution());
                        std::vector<double> pde_solution;
                        for (unsigned i=0; i<mpCoarsePdeMesh->GetNumNodes(); i++)
                        {
                            pde_solution.push_back(solution_repl[i]);
                        }

                        p_vtk_mesh_writer->AddPointData(p_pde_and_bc->rGetDependentVariableName(),pde_solution);
                    }
                }

#endif //CHASTE_VTK

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
                        // should only come into this method AFTER solving the PDE
                        NEVER_REACHED;
                    }
                }
            }
            else // Not coarse mesh
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
                    double solution = cell_iter->GetCellData()->GetItem(mPdeAndBcCollection[pde_index]->rGetDependentVariableName());
                    (*mpVizPdeSolutionResultsFile) << solution << " ";
                }
            }
        }
        (*mpVizPdeSolutionResultsFile) << "\n";
#ifdef CHASTE_VTK
        if (p_vtk_mesh_writer != NULL)
        {
            p_vtk_mesh_writer->WriteFilesUsingMesh(*mpCoarsePdeMesh);
            delete p_vtk_mesh_writer;
        }
#endif //CHASTE_VTK
    }
}

template<unsigned DIM>
void CellBasedPdeHandler<DIM>::SetWriteAverageRadialPdeSolution(const std::string& rName, unsigned numRadialIntervals, bool writeDailyResults)
{
    mAverageRadialSolutionVariableName = rName;
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
                average_solution += (iter->second)->GetCellData()->GetItem(mAverageRadialSolutionVariableName);
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
double CellBasedPdeHandler<DIM>::GetPdeSolutionAtPoint(const c_vector<double,DIM>& rPoint, const std::string& rVariable)
{
    double solution_at_point = 0.0;

    unsigned pde_index = UINT_MAX;

    // Loop over elements of mPdeAndBcCollection to find correct PDE
    for (unsigned i=0; i<mPdeAndBcCollection.size(); i++)
    {
        if (mPdeAndBcCollection[i]->rGetDependentVariableName() == rVariable)
        {
            pde_index = i;
            break;
        }
    }
    if (pde_index == UINT_MAX)
    {
        EXCEPTION("Tried to get the solution of a variable name: " + rVariable + ". There is no PDE with that variable.");
    }
    PdeAndBoundaryConditions<DIM>* p_pde_and_bc = mPdeAndBcCollection[pde_index];

    Element<DIM,DIM>* p_containing_element;

    if (mpCoarsePdeMesh != NULL)
    {
        // Find PDE element containing point
        unsigned elem_index = mpCoarsePdeMesh->GetContainingElementIndex(ChastePoint<DIM>(rPoint));
        p_containing_element = mpCoarsePdeMesh->GetElement(elem_index);
    }
    else // Tetrahedral mesh
    {
        // If not using a coarse PDE mesh, we must be using a MeshBasedCellPopulation
        TetrahedralMesh<DIM,DIM>* p_tetrahedral_mesh = &(static_cast<MeshBasedCellPopulation<DIM>*>(mpCellPopulation)->rGetMesh());

        unsigned elem_index = p_tetrahedral_mesh->GetContainingElementIndex(ChastePoint<DIM>(rPoint));
        p_containing_element = p_tetrahedral_mesh->GetElement(elem_index);
    }

    // Interpolate solution
    if (p_pde_and_bc->GetSolution())
    {
        ReplicatableVector solution_repl(p_pde_and_bc->GetSolution());
        c_vector<double,DIM+1> weights = p_containing_element->CalculateInterpolationWeights(rPoint);
        for (unsigned i=0; i<DIM+1; i++)
        {
            double nodal_value = solution_repl[p_containing_element->GetNodeGlobalIndex(i)];
            solution_at_point += nodal_value * weights(i);
        }
    }
    else
    {
        // should only come into this method AFTER solving the PDE
        NEVER_REACHED;
    }

    return solution_at_point;
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
bool CellBasedPdeHandler<DIM>::PdeSolveNeedsCoarseMesh()
{
    return ((dynamic_cast<NodeBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL)
            || (dynamic_cast<PottsBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL)
            || (dynamic_cast<CaBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL)
            || (dynamic_cast<VertexBasedCellPopulation<DIM>*>(mpCellPopulation) != NULL));
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellBasedPdeHandler)

///////// Explicit instantiation
template class CellBasedPdeHandler<1>;
template class CellBasedPdeHandler<2>;
template class CellBasedPdeHandler<3>;
