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

#include "AbstractBoxDomainPdeModifier.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"

template<unsigned DIM>
AbstractBoxDomainPdeModifier<DIM>::AbstractBoxDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
                                                                boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
                                                                bool isNeumannBoundaryCondition,
                                                                boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
                                                                double stepSize,
                                                                Vec solution)
    : AbstractPdeModifier<DIM>(pPde,
                               pBoundaryCondition,
                               isNeumannBoundaryCondition,
                               solution),
      mpMeshCuboid(pMeshCuboid),
      mStepSize(stepSize),
      mSetBcsOnBoxBoundary(true),
      mSetBcsOnBoundingSphere(false),
      mUseVoronoiCellsForInterpolation(false),
      mTypicalCellRadius(0.5) // defaults to 0.5
{
    if (pMeshCuboid)
    {
        // We only need to generate mpFeMesh once, as it does not vary with time
        this->GenerateFeMesh(mpMeshCuboid, mStepSize);
        this->mDeleteFeMesh = true;
        //initialise the boundary nodes
        this->mIsDirichletBoundaryNode = std::vector<double>(this->mpFeMesh->GetNumNodes(), 0.0);
    }
}

template<unsigned DIM>
AbstractBoxDomainPdeModifier<DIM>::~AbstractBoxDomainPdeModifier()
{
}

template<unsigned DIM>
double AbstractBoxDomainPdeModifier<DIM>::GetStepSize()
{
     return mStepSize;
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary)
{
    mSetBcsOnBoxBoundary = setBcsOnBoxBoundary;
}

template<unsigned DIM>
bool AbstractBoxDomainPdeModifier<DIM>::AreBcsSetOnBoxBoundary()
{
    return mSetBcsOnBoxBoundary;
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::SetBcsOnBoundingSphere(bool setBcsOnBoundingSphere)
{
    mSetBcsOnBoundingSphere = setBcsOnBoundingSphere;
}

template<unsigned DIM>
bool AbstractBoxDomainPdeModifier<DIM>::AreBcsSetOnBoundingSphere()
{
    return mSetBcsOnBoundingSphere;
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::SetUseVoronoiCellsForInterpolation(bool useVoronoiCellsForInterpolation)
{
    mUseVoronoiCellsForInterpolation = useVoronoiCellsForInterpolation;
}

template<unsigned DIM>
bool AbstractBoxDomainPdeModifier<DIM>::GetUseVoronoiCellsForInterpolation()
{
    return mUseVoronoiCellsForInterpolation;
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::SetTypicalCellRadius(double typicalCellRadius)
{
    assert(mTypicalCellRadius>=0.0);
    mTypicalCellRadius = typicalCellRadius;
}

template<unsigned DIM>
double AbstractBoxDomainPdeModifier<DIM>::GetTypicalCellRadius()
{
    return mTypicalCellRadius;
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainerHelper(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                                                   std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > pBcc)
{
    if (!this->mSetBcsOnBoxBoundary)
    {
        // Here we approximate the cell population by the bounding spehere and apply the boundary conditions outside the sphere.
        if (this->mSetBcsOnBoundingSphere)
        {
            // First find the centre of the tissue by choosing the mid point of the extrema.
            c_vector<double, DIM> tissue_maxima = zero_vector<double>(DIM);
            c_vector<double, DIM> tissue_minima = zero_vector<double>(DIM);
            for (unsigned i = 0; i < DIM; i++)
            {
                tissue_maxima[i] = -DBL_MAX;
                tissue_minima[i] = DBL_MAX;
            }
            
            
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
            {
                const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                
                for (unsigned i = 0; i < DIM; i++)
                {
                    if (r_position_of_cell[i] > tissue_maxima[i])
                    {
                        tissue_maxima[i] = r_position_of_cell[i];
                    }
                    if (r_position_of_cell[i] < tissue_minima[i])
                    {
                        tissue_minima[i] = r_position_of_cell[i];
                    }
                }               
            }

            c_vector<double, DIM> tissue_centre = 0.5*(tissue_maxima + tissue_minima);


            double tissue_radius = 0.0;
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
            {
                const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                double radius = norm_2(tissue_centre - r_position_of_cell.rGetLocation());
                
                if (tissue_radius < radius)
                {
                    tissue_radius = radius;
                }                
            }

            // Apply boundary condition to the nodes outside the tissue_radius
            if (this->IsNeumannBoundaryCondition())
            {
                NEVER_REACHED;
            }
            else
            {
                // Impose any Dirichlet boundary conditions
                for (unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++)
                {
                    double radius = norm_2(tissue_centre - this->mpFeMesh->GetNode(i)->rGetLocation());
                    
                    if (radius > tissue_radius)
                    {
                        pBcc->AddDirichletBoundaryCondition(this->mpFeMesh->GetNode(i), this->mpBoundaryCondition.get(), 0, false);
                        this->mIsDirichletBoundaryNode[i] = 1.0;
                    }
                }
            }
        }
        else // Set pde nodes as boundary node if elements dont contain cells or nodes aren't within 0.5CD of a cell centre
        {
            // Get the set of coarse element indices that contain cells
            std::set<unsigned> coarse_element_indices_in_map;
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
            {
                coarse_element_indices_in_map.insert(this->mCellPdeElementMap[*cell_iter]);
            }

            // Find the node indices associated with elements whose indices are NOT in the set coarse_element_indices_in_map
            std::set<unsigned> coarse_mesh_boundary_node_indices;
            for (unsigned i=0; i<this->mpFeMesh->GetNumElements(); i++)
            {
                if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
                {
                    Element<DIM,DIM>* p_element = this->mpFeMesh->GetElement(i);
                    for (unsigned j=0; j<DIM+1; j++)
                    {
                        unsigned node_index = p_element->GetNodeGlobalIndex(j);
                        coarse_mesh_boundary_node_indices.insert(node_index);
                    }
                }
            }

            // Also remove nodes that are within the typical cell radius from the centre of a cell.
            std::set<unsigned> nearby_node_indices;
            for (std::set<unsigned>::iterator node_iter = coarse_mesh_boundary_node_indices.begin();
                node_iter != coarse_mesh_boundary_node_indices.end();
                ++node_iter)
            {
                bool remove_node = false;

                c_vector<double,DIM> node_location = this->mpFeMesh->GetNode(*node_iter)->rGetLocation();

                for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
                {
                    const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                    double separation = norm_2(node_location - r_position_of_cell.rGetLocation());

                    if (separation <= mTypicalCellRadius)
                    {
                        remove_node = true;
                        break;
                    }                
                }

                if (remove_node)
                {
                    // Node near cell so set it to be removed from boundary set
                    nearby_node_indices.insert(*node_iter);
                }
            }

            // Remove nodes that are near cells from boundary set
            for (std::set<unsigned>::iterator node_iter = nearby_node_indices.begin();
                node_iter != nearby_node_indices.end();
                ++node_iter)
            {
                coarse_mesh_boundary_node_indices.erase(*node_iter);
            }


            // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
            if (this->IsNeumannBoundaryCondition())
            {
                NEVER_REACHED;
            }
            else
            {
                // Impose any Dirichlet boundary conditions
                for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
                    iter != coarse_mesh_boundary_node_indices.end();
                    ++iter)
                {
                    pBcc->AddDirichletBoundaryCondition(this->mpFeMesh->GetNode(*iter), this->mpBoundaryCondition.get(), 0, false);
                    this->mIsDirichletBoundaryNode[*iter] = 1.0;
                }
            }
        }
    }
    else // Apply BC at boundary of box domain FE mesh
    {
        if (this->IsNeumannBoundaryCondition())
        {
            // Impose any Neumann boundary conditions
            for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
                 elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
                 ++elem_iter)
            {
                pBcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryCondition.get());
            }
        }
        else
        {
            // Impose any Dirichlet boundary conditions
            for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                 node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                 ++node_iter)
            {
                pBcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryCondition.get());
                this->mIsDirichletBoundaryNode[(*node_iter)->GetIndex()] = 1.0;
            }
        }
    }
}



template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractPdeModifier<DIM>::SetupSolve(rCellPopulation, outputDirectory);

    InitialiseCellPdeElementMap(rCellPopulation);
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::GenerateFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize)
{
    // Create a regular coarse tetrahedral mesh
    this->mpFeMesh = new TetrahedralMesh<DIM,DIM>();
    
    GenerateAndReturnFeMesh(pMeshCuboid, stepSize, this->mpFeMesh);
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::GenerateAndReturnFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize, TetrahedralMesh<DIM,DIM>* pMesh)
{
    // Create a regular coarse tetrahedral mesh
    switch (DIM)
    {
        case 1:
            pMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0));
            break;
        case 2:
            pMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1));
            break;
        case 3:
            pMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1), pMeshCuboid->GetWidth(2));
            break;
        default:
            NEVER_REACHED;
    }

    // Get centroid of meshCuboid
    ChastePoint<DIM> upper = pMeshCuboid->rGetUpperCorner();
    ChastePoint<DIM> lower = pMeshCuboid->rGetLowerCorner();
    c_vector<double,DIM> centre_of_cuboid = 0.5*(upper.rGetLocation() + lower.rGetLocation());

    // Find the centre of the PDE mesh
    c_vector<double,DIM> centre_of_coarse_mesh = zero_vector<double>(DIM);
    for (unsigned i=0; i<pMesh->GetNumNodes(); i++)
    {
        centre_of_coarse_mesh += pMesh->GetNode(i)->rGetLocation();
    }
    centre_of_coarse_mesh /= pMesh->GetNumNodes();

    // Now move the mesh to the correct location
    pMesh->Translate(centre_of_cuboid - centre_of_coarse_mesh);
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution);

    if (mUseVoronoiCellsForInterpolation) 
    {
        unsigned num_nodes = rCellPopulation.GetNumNodes();

        std::vector<double> cell_data(num_nodes, -1);
        std::vector<unsigned> num_cells(num_nodes, -1);
        
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            unsigned cell_location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            cell_data[cell_location_index]=0.0;
            num_cells[cell_location_index]=0;    
        }

        // Loop over nodes of the finite element mesh and work out which voronoi region the node is in.
        for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
                node_iter != this->mpFeMesh->GetNodeIteratorEnd();
                ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();

            c_vector<double,DIM> node_location = node_iter->rGetLocation();

            double closest_separation = DBL_MAX;
            unsigned nearest_cell = UNSIGNED_UNSET;

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
            {
                c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                double separation = norm_2(node_location - cell_location);

                if (separation < closest_separation)
                {
                    closest_separation = separation;
                    nearest_cell = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                }                
            }
            assert(closest_separation<DBL_MAX);

            cell_data[nearest_cell] = cell_data[nearest_cell] + solution_repl[node_index];
            num_cells[nearest_cell] = num_cells[nearest_cell] + 1;
        }   
        
        // Now calculate the solution in the cell by averaging over all nodes in the voronoi region.
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
        {
            unsigned cell_location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);   

            
            if (num_cells[cell_location_index]==0)
            {
                EXCEPTION("One or more of the cells doesnt contain any pde nodes so cant use voroni CellData calculation in the ");
            }

            double  solution_at_cell = cell_data[cell_location_index]/num_cells[cell_location_index];
            
            cell_iter->GetCellData()->SetItem(this->mDependentVariableName, solution_at_cell);
        }

        if (this->mOutputGradient)
        {
            // This isnt implemented yet
            NEVER_REACHED;
        }
    }  
    else // Interpolate solutions 
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            // The cells are not nodes of the mesh, so we must interpolate
            double solution_at_cell = 0.0;

            // Find the element in the FE mesh that contains this cell. CellElementMap has been updated so use this.
            unsigned elem_index = mCellPdeElementMap[*cell_iter];
            Element<DIM,DIM>* p_element = this->mpFeMesh->GetElement(elem_index);

            const ChastePoint<DIM>& node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(node_location);

            for (unsigned i=0; i<DIM+1; i++)
            {
                double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(i)];
                solution_at_cell += nodal_value * weights(i);
            }

            cell_iter->GetCellData()->SetItem(this->mDependentVariableName, solution_at_cell);

            if (this->mOutputGradient)
            {
                // Now calculate the gradient of the solution and store this in CellVecData
                c_vector<double, DIM> solution_gradient = zero_vector<double>(DIM);

                // Calculate the basis functions at any point (e.g. zero) in the element
                c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
                double jacobian_det;
                this->mpFeMesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
                const ChastePoint<DIM> zero_point;
                c_matrix<double, DIM, DIM+1> grad_phi;
                LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

                for (unsigned node_index=0; node_index<DIM+1; node_index++)
                {
                    double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(node_index)];

                    for (unsigned j=0; j<DIM; j++)
                    {
                        solution_gradient(j) += nodal_value* grad_phi(j, node_index);
                    }
                }

                switch (DIM)
                {
                    case 1:
                        cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_x", solution_gradient(0));
                        break;
                    case 2:
                        cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_x", solution_gradient(0));
                        cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_y", solution_gradient(1));
                        break;
                    case 3:
                        cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_x", solution_gradient(0));
                        cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_y", solution_gradient(1));
                        cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_z", solution_gradient(2));
                        break;
                    default:
                        NEVER_REACHED;
                }
            }
        }
    }     
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::InitialiseCellPdeElementMap(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    mCellPdeElementMap.clear();

    // Find the element of mpFeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndex(r_position_of_cell);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::UpdateCellPdeElementMap(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Find the element of mpCoarsePdeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndexWithInitialGuess(r_position_of_cell, mCellPdeElementMap[*cell_iter]);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractBoxDomainPdeModifier<1>;
template class AbstractBoxDomainPdeModifier<2>;
template class AbstractBoxDomainPdeModifier<3>;
