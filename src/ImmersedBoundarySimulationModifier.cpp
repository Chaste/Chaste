/*

Copyright (c) 2005-2014, University of Oxford.
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

#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "Exception.hpp"
#include <complex>
#include <fftw3.h>
#include "Debug.hpp"

template<unsigned DIM>
ImmersedBoundarySimulationModifier<DIM>::ImmersedBoundarySimulationModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mpMesh(NULL),
      mpCellPopulation(NULL),
      mCellNeighbourUpdateFrequency(10u),
      mNumGridPtsX(0u),
      mNumGridPtsY(0u),
      mGridSpacingX(0.0),
      mGridSpacingY(0.0),
      mFftNorm(0.0),
      mDiffusionCoef(1.0),
      mFluidDensity(1.0),
      mI(0.0 + 1.0 * 1i)
{
}

template<unsigned DIM>
ImmersedBoundarySimulationModifier<DIM>::~ImmersedBoundarySimulationModifier()
{
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // We need to update cell neighbours occasionally, but not at every timestep
    if(SimulationTime::Instance()->GetTimeStepsElapsed() % mCellNeighbourUpdateFrequency == 0)
    {
        this->UpdateCellNeighbours();
    }

    // This will solve the fluid problem for all timesteps after the first, which is handled in SetupSolve()
    this->UpdateFluidVelocityGrids(rCellPopulation);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // We can set up some helper variables here which need only be set up once for the entire simulation
    this->SetupConstantMemberVariables(rCellPopulation);
    this->SetupTransmembraneProteinLevels();
    this->UpdateCellNeighbours();

    // This will solve the fluid problem based on the initial mesh setup
    this->UpdateFluidVelocityGrids(rCellPopulation);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateFluidVelocityGrids(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    this->ClearForces();
    this->CalculatePerimeterElasticForces();
    this->CalculateCellCellInteractionElasticForces();
    this->PropagateForcesToFluidGrid();
    this->SolveNavierStokesSpectral();
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupConstantMemberVariables(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<ImmersedBoundaryCellPopulation<DIM> *>(&rCellPopulation) == NULL)
    {
        EXCEPTION("Cell population must be Immersed Boundary");
    }

    mpCellPopulation = static_cast<ImmersedBoundaryCellPopulation<DIM> *>(&rCellPopulation);
    mpMesh = &(mpCellPopulation->rGetMesh());

    // Get the size of the mesh
    mNumGridPtsX = mpMesh->GetNumGridPtsX();
    mNumGridPtsY = mpMesh->GetNumGridPtsY();

    // Get the grid spacing
    mGridSpacingX = 1.0 / (double) mNumGridPtsX;
    mGridSpacingY = 1.0 / (double) mNumGridPtsY;

    // Set up force grids
    SetupGrid(mFluidForceGridX);
    SetupGrid(mFluidForceGridY);

    // FFT norm
    mFftNorm = sqrt((double) mNumGridPtsX * (double) mNumGridPtsY);

    // Create sine variables
    mSinX.resize(mNumGridPtsX);
    mSin2X.resize(mNumGridPtsX);
    mSinY.resize(mNumGridPtsY);
    mSin2Y.resize(mNumGridPtsY);
    for (unsigned x = 0; x < mNumGridPtsX; x++)
    {
        mSinX[x] = sin(M_PI * (double) x * mGridSpacingX);
        mSin2X[x] = sin(2 * M_PI * (double) x * mGridSpacingX);
    }
    for (unsigned y = 0; y < mNumGridPtsY; y++)
    {
        mSinY[y] = sin(M_PI * (double) y * mGridSpacingY);
        mSin2Y[y] = sin(2 * M_PI * (double) y * mGridSpacingY);
    }

    mInteractionDistance = mpCellPopulation->GetInteractionDistance();
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupTransmembraneProteinLevels()
{
    unsigned num_nodes = mpMesh->GetNumNodes();

    std::vector<double>& r_ecad_levels = mpCellPopulation->rGetModifiableECadherinLevels();
    std::vector<double>& r_pcad_levels = mpCellPopulation->rGetModifiablePCadherinLevels();
    std::vector<double>& r_integrin_levels = mpCellPopulation->rGetModifiableIntegrinLevels();

    r_ecad_levels.resize(num_nodes);
    r_pcad_levels.resize(num_nodes);
    r_integrin_levels.resize(num_nodes);

    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_iter = mpMesh->GetElementIteratorBegin();
         elem_iter != mpMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_idx = elem_iter->GetIndex();
        unsigned num_nodes_this_elem = elem_iter->GetNumNodes();

        for (unsigned node_idx = 0 ; node_idx < num_nodes_this_elem ; node_idx++)
        {
            // Current node
            Node<DIM>* current_node = elem_iter->GetNode(node_idx);
            unsigned global_idx = current_node->GetIndex();

            if (elem_idx == mpMesh->GetMembraneIndex())
            {
                r_integrin_levels[global_idx] = 10.0;
                r_ecad_levels[global_idx] = 0.0;
                r_pcad_levels[global_idx] = 0.0;
            }
            else // element is not the membrane element
            {
                r_integrin_levels[global_idx] = 0.0;
                r_ecad_levels[global_idx] = 1.0;
                r_pcad_levels[global_idx] = 0.0;
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateTransmembraneProteinLevels()
{
    //\todo decide how protein levels should be updated each timestep.  This can be overridden in derived classes
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateCellNeighbours()
{
    // This method needs to be called first, as centroids are not calculated each timestep
    mpCellPopulation->UpdateCellCentroids();

    // Get the neighbouring indices for each cell
    for (typename ImmersedBoundaryCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
         cell_iter != mpCellPopulation->End();
         ++cell_iter)
    {
        //\todo: when cells divide etc, may need to call this method directly
        mCellNeighbours[*cell_iter] = mpCellPopulation->GetNeighbouringLocationIndices(*cell_iter);
    }

    /**
     * Next, we can prune the set so we're not duplicating work calculating cell-cell interactions.
     *
     * First, loop over each cell.  For each cell, see which cells are neighbours to it.  For each of these cells,
     * check their entry in mCellNeighbours for the key of the current cell in the outer loop.  If it's there,
     * we can remove it.
     */
    for (typename ImmersedBoundaryCellPopulation<DIM>::Iterator cell_iter = mpCellPopulation->Begin();
         cell_iter != mpCellPopulation->End();
         ++cell_iter)
    {
        // Get the location index of our reference cell
        unsigned location_index = mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        // Loop over each neighbour
        for (typename std::set<unsigned>::iterator set_iter = mCellNeighbours[*cell_iter].begin();
             set_iter != mCellNeighbours[*cell_iter].end();
             ++set_iter)
        {
            // Grab a pointer to the cell, which can be used as a key in the map mCellNeighbours
            CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(*set_iter);

            // Erase the element of the set which coincides with location index of our reference cell
            // If no such index is found, it should still be fine.
            mCellNeighbours[p_cell].erase(mCellNeighbours[p_cell].find(location_index));
        }
    }
    MARK;
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::ClearForces()
{
    for (typename ImmersedBoundaryMesh<DIM, DIM>::NodeIterator node_iter = mpMesh->GetNodeIteratorBegin(false);
            node_iter != mpMesh->GetNodeIteratorEnd();
            ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

    for (unsigned y = 0; y < mNumGridPtsY; y++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            mFluidForceGridX[y][x] = 0.0;
            mFluidForceGridY[y][x] = 0.0;
        }
    }
}


template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::CalculatePerimeterElasticForces()
{
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_iter = mpMesh->GetElementIteratorBegin();
            elem_iter != mpMesh->GetElementIteratorEnd();
            ++elem_iter)
    {
        // Get number of nodes in current element
        unsigned num_nodes = elem_iter->GetNumNodes();
        assert(num_nodes > 0);

        // Get spring parameters (owned by the elements as they may differ within the element population)
        double spring_constant = elem_iter->GetMembraneSpringConstant();
        double rest_length = elem_iter->GetMembraneRestLength();

        // Helper variables
        double normed_dist;
        c_vector<double, DIM> aggregate_force;

        // Make a vector to store the force on node i+1 from node i
        std::vector<c_vector<double, DIM> > elastic_force_to_next_node(num_nodes);

        // Loop over nodes and calculate the force exerted on node i+1 by node i
        for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
        {
            // Index of the next node, calculated modulo number of nodes in this element
            unsigned next_idx = (node_idx + 1) % num_nodes;

            // Hooke's law linear spring force
            elastic_force_to_next_node[node_idx] = mpMesh->GetVectorFromAtoB(elem_iter->GetNodeLocation(next_idx), elem_iter->GetNodeLocation(node_idx));
            normed_dist = norm_2(elastic_force_to_next_node[node_idx]);
            elastic_force_to_next_node[node_idx] *= spring_constant * (normed_dist - rest_length) / normed_dist;
        }

        // Add the contributions of springs springs adjacent to each node
        for (unsigned node_idx = 0 ; node_idx < num_nodes ; node_idx++)
        {
            // Index of previous node, but -1%n doesn't work, so we add num_nodes when calculating
            unsigned prev_idx = (node_idx + num_nodes - 1) % num_nodes;

            aggregate_force = elastic_force_to_next_node[prev_idx] - elastic_force_to_next_node[node_idx];

            // Add the aggregate force contribution to the node
            elem_iter->GetNode(node_idx)->AddAppliedForceContribution(aggregate_force);
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::CalculateCellCellInteractionElasticForces()
{
    /**
     * We iterate over each cell, and then over each neighbour, by using the map mCellNeighbours.
     * For each neighbour, we iterate over each node in the current cell and then each node in the
     * neighbouring cell.
     */

    // Helper variables for loops
    double normed_dist;
    c_vector<double, DIM> vector_between_nodes;
    c_vector<double, DIM> force_a_b;
    c_vector<double, DIM> force_b_a;
    double protein_mult;
    Node<DIM>* p_curr_node;
    Node<DIM>* p_comp_node;

    // Transmembrane protein levels
    const std::vector<double>& r_ecad_levels = mpCellPopulation->rGetECadherinLevels();
    const std::vector<double>& r_pcad_levels = mpCellPopulation->rGetPCadherinLevels();
    const std::vector<double>& r_integrin_levels = mpCellPopulation->rGetIntegrinLevels();

    for (typename std::map<CellPtr, std::set<unsigned> >::iterator cell_iter = mCellNeighbours.begin();
         cell_iter != mCellNeighbours.end();
         ++cell_iter)
    {
        // cell_iter->first  (key value) is a CellPtr
        // cell_iter->second (mapped value) is a std::set<unsigned>; indices of cells which neighbour the CellPtr

        // Get the current element and associated information
        ImmersedBoundaryElement<DIM,DIM>* p_curr_elem = mpCellPopulation->GetElementCorrespondingToCell(cell_iter->first);
        double spring_constant = p_curr_elem->GetCellCellSpringConstant();
        double rest_length = p_curr_elem->GetCellCellRestLength();
        unsigned curr_num_nodes = p_curr_elem->GetNumNodes();

        // Iterate over neighbours, by using the mapped value from mCellNeighbours
        for (typename std::set<unsigned>::iterator set_iter = (cell_iter->second).begin();
             set_iter != (cell_iter->second).end();
             ++set_iter)
        {
            // *set_iter is an unsigned int; index of current cell neighbouring CellPtr

            // Get the comparison-element and its number of nodes
            ImmersedBoundaryElement<DIM,DIM>* p_comp_elem = mpCellPopulation->GetElementCorrespondingToCell(mpCellPopulation->GetCellUsingLocationIndex(*set_iter));
            unsigned comp_num_nodes = p_comp_elem->GetNumNodes();
            double comp_s_c = p_comp_elem->GetCellCellSpringConstant();
            double comp_r_l = p_comp_elem->GetCellCellRestLength();

            // Loop over nodes in current element, then over nodes in comparison element
            for (unsigned curr_idx = 0 ; curr_idx < curr_num_nodes ; curr_idx++)
            {
                p_curr_node = p_curr_elem->GetNode(curr_idx);

                // Pre-define current node location so we're not grabbing it more times than necessary
                c_vector<double, DIM> curr_location = p_curr_node->rGetLocation();

                for (unsigned comp_idx = 0 ; comp_idx < comp_num_nodes ; comp_idx++)
                {
                    p_comp_node = p_comp_elem->GetNode(comp_idx);

                    vector_between_nodes = mpMesh->GetVectorFromAtoB(curr_location, p_comp_node->rGetLocation());
                    normed_dist = norm_2(vector_between_nodes);

                    // Apply linear spring cell-cell interaction if dist is smaller than the pre-defined interaction distance
                    if (normed_dist < mInteractionDistance)
                    {
                        // The protein multiplier is a function of the levels of each protein in the current and comparison nodes
                        protein_mult = std::min(r_ecad_levels[p_curr_node->GetIndex()], r_ecad_levels[p_comp_node->GetIndex()]) +
                                       std::min(r_pcad_levels[p_curr_node->GetIndex()], r_pcad_levels[p_comp_node->GetIndex()]) +
                                       std::max(r_integrin_levels[p_curr_node->GetIndex()], r_integrin_levels[p_comp_node->GetIndex()]);

                        force_a_b = vector_between_nodes * spring_constant * protein_mult * (normed_dist - rest_length) / normed_dist;
                        p_curr_node->AddAppliedForceContribution(force_a_b);

                        force_b_a = -1.0 * vector_between_nodes * comp_s_c * protein_mult * (normed_dist - comp_r_l) / normed_dist;
                        p_comp_node->AddAppliedForceContribution(force_b_a);
                    }
                }
            }
        }
    }

    // Now add force contributions from each AbstractForce
    for (typename std::vector<boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > >::iterator iter = mForceCollection.begin();
            iter != mForceCollection.end();
            ++iter)
    {
        (*iter)->AddForceContribution(*(this->mpCellPopulation));
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::PropagateForcesToFluidGrid()
{
    // Helper variables
    double dl = mpMesh->GetCharacteristicNodeSpacing();
    double dist_x;
    double dist_y;
    double force_x;
    double force_y;
    int first_idx_x;
    int first_idx_y;

    // Iterate over all nodes and grab their position
    for (typename ImmersedBoundaryMesh<DIM, DIM>::NodeIterator node_iter = mpMesh->GetNodeIteratorBegin(false);
            node_iter != mpMesh->GetNodeIteratorEnd();
            ++node_iter)
    {

        // Get location and applied force contribution of current node
        c_vector<double, DIM> node_location = node_iter->rGetLocation();
        c_vector<double, DIM> applied_force = node_iter->rGetAppliedForce();


        // Get index of grid positions ignoring possible wrap-around
        first_idx_x = (int)floor(node_location[0] / mGridSpacingX) - 1;
        first_idx_y = (int)floor(node_location[1] / mGridSpacingY) - 1;

        // Loop over the 4x4 grid used to spread the force on the nodes to the fluid grid
        for (unsigned x_idx = 0 ; x_idx < 4 ; x_idx ++)
        {
            // Calculate distance between current x index and node, then account for possible wrap-around
            dist_x = fabs((double)(first_idx_x + x_idx) * mGridSpacingX - node_location[0]);
            if(first_idx_x == -1)
            {
                first_idx_x += mNumGridPtsX;
            }

            for (unsigned y_idx = 0 ; y_idx < 4 ; y_idx ++)
            {
                // Calculate distance between current x index and node, then account for possible wrap-around
                dist_y = fabs((double)(first_idx_y + y_idx) * mGridSpacingY - node_location[1]);
                if(first_idx_y == -1)
                {
                    first_idx_y += mNumGridPtsY;
                }

                // The applied force is weighted by the delta function
                force_x = applied_force[0] * Delta1D(dist_x, mGridSpacingX) * Delta1D(dist_y, mGridSpacingY) * dl;
                force_y = applied_force[1] * Delta1D(dist_x, mGridSpacingX) * Delta1D(dist_y, mGridSpacingY) * dl;

                mFluidForceGridX[(first_idx_y + y_idx) % mNumGridPtsY][(first_idx_x + x_idx) % mNumGridPtsX] += force_x;
                mFluidForceGridY[(first_idx_y + y_idx) % mNumGridPtsY][(first_idx_x + x_idx) % mNumGridPtsX] += force_y;
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SolveNavierStokesSpectral()
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    // Get non modifiable fluid grids
    const std::vector<std::vector<double> >& VelX = mpMesh->rGetFluidVelocityGridX();
    const std::vector<std::vector<double> >& VelY = mpMesh->rGetFluidVelocityGridY();

    // Create upwind variables
    std::vector<std::vector<double> > Upwind_x;
    std::vector<std::vector<double> > Upwind_y;
    UpwindScheme(VelX, VelY, Upwind_x, Upwind_y);

    // Create RHS of linear system
    std::vector<std::vector<double> > rhsX; SetupGrid(rhsX);
    std::vector<std::vector<double> > rhsY; SetupGrid(rhsY);

    for(unsigned y = 0 ; y < mNumGridPtsY ; y++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            rhsX[y][x] = VelX[y][x] + dt * (mFluidForceGridX[y][x] - Upwind_x[y][x]);
            rhsY[y][x] = VelY[y][x] + dt * (mFluidForceGridY[y][x] - Upwind_y[y][x]);
        }
    }

    // Perform fft on rhsX
    std::vector<std::vector<std::complex<double> > > VelX_hat;
    Fft2DForwardRealToComplex(rhsX, VelX_hat);

    // Perform fft on rhsY
    std::vector<std::vector<std::complex<double> > > VelY_hat;
    Fft2DForwardRealToComplex(rhsY, VelY_hat);

    std::vector<std::vector<std::complex<double> > > p_hat; SetupGrid(p_hat);

    // Calculate pressure
    for( unsigned y = 0 ; y < mNumGridPtsY ; y++ )
    {
        for( unsigned x = 0 ; x < mNumGridPtsX ; x++ )
        {
            std::complex<double> numerator = -mI * (mSin2X[x] * VelX_hat[y][x] / mGridSpacingX + mSin2Y[y] * VelY_hat[y][x] / mGridSpacingY);

            double denominator = (mSin2X[x] * mSin2X[x] / (mGridSpacingX * mGridSpacingX)) + (mSin2Y[y] * mSin2Y[y] / (mGridSpacingY * mGridSpacingY));
            denominator *= (dt / mFluidDensity);

            p_hat[y][x] = numerator / denominator;
        }
    }

    // Set some values to zero
    p_hat[0][0] = 0.0;
    p_hat[0][mNumGridPtsX/2] = 0.0;
    p_hat[mNumGridPtsY/2][mNumGridPtsX/2] = 0.0;
    p_hat[mNumGridPtsY/2][0] = 0.0;

    // Do final stage of computation before inverse FFT
    std::vector<std::vector<std::complex<double> > > pre_inverse_X; SetupGrid(pre_inverse_X);
    std::vector<std::vector<std::complex<double> > > pre_inverse_Y; SetupGrid(pre_inverse_Y);

    for( unsigned y = 0 ; y < mNumGridPtsY ; y++ )
    {
        for( unsigned x = 0 ; x < mNumGridPtsX ; x++ )
        {
            double op = 1 + (4 * mDiffusionCoef * dt / mFluidDensity) * ( (mSinX[x] * mSinX[x] / (mGridSpacingX * mGridSpacingX)) + (mSinY[y] * mSinY[y] / (mGridSpacingY * mGridSpacingY)));
            pre_inverse_X[y][x] = (VelX_hat[y][x] - (mI * dt / (mFluidDensity * mGridSpacingX)) * mSin2X[x] * p_hat[y][x]) / op;
            pre_inverse_Y[y][x] = (VelY_hat[y][x] - (mI * dt / (mFluidDensity * mGridSpacingY)) * mSin2Y[y] * p_hat[y][x]) / op;
        }
    }

    // Perform inverse FFT on x data
    std::vector<std::vector<double> > new_velocity_x;
    Fft2DInverseComplexToReal(pre_inverse_X, new_velocity_x);

    // Perform inverse FFT on y data
    std::vector<std::vector<double> > new_velocity_y;
    Fft2DInverseComplexToReal(pre_inverse_Y, new_velocity_y);

    // Get modifiable fluid grids
    std::vector<std::vector<double> >& modifiable_vel_x = mpMesh->rGetModifiableFluidVelocityGridX();
    std::vector<std::vector<double> >& modifiable_vel_y = mpMesh->rGetModifiableFluidVelocityGridY();

    for(unsigned y = 0 ; y < mNumGridPtsY ; y++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            modifiable_vel_x[y][x] = new_velocity_x[y][x];
            modifiable_vel_y[y][x] = new_velocity_y[y][x];
        }
    }
}

template<unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::Delta1D(double dist, double spacing)
{
    return (0.25 * (1.0 + cos(M_PI * dist / (2 * spacing)))) / spacing;
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpwindScheme(const std::vector<std::vector<double> >& in_x, const std::vector<std::vector<double> >& in_y, std::vector<std::vector<double> >& out_x, std::vector<std::vector<double> >& out_y)
{
    // Ensure the derivative grid is set up to the correct size
    SetupGrid(out_x);
    SetupGrid(out_y);

    unsigned prev_x = mNumGridPtsX - 1;
    unsigned prev_y = mNumGridPtsY - 1;

    unsigned next_x = 1;
    unsigned next_y = 1;

    for(unsigned y = 0 ; y < mNumGridPtsY ; y++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            // Set values for output from conditional on x grid
            if(in_x[y][x] > 0)
            {
                out_x[y][x] = in_x[y][x] * (in_x[y][x] - in_x[y][prev_x]) / mGridSpacingX;
                out_y[y][x] = in_x[y][x] * (in_y[y][x] - in_y[y][prev_x]) / mGridSpacingX;
            }
            else
            {
                out_x[y][x] = in_x[y][x] * (in_x[y][next_x] - in_x[y][x]) / mGridSpacingX;
                out_y[y][x] = in_x[y][x] * (in_y[y][next_x] - in_y[y][x]) / mGridSpacingX;
            }

            // Then add values from conditional on y grid
            if(in_y[y][x] > 0)
            {
                out_x[y][x] += in_y[y][x] * (in_x[y][x] - in_x[prev_y][x]) / mGridSpacingY;
                out_y[y][x] += in_y[y][x] * (in_y[y][x] - in_y[prev_y][x]) / mGridSpacingY;
            }
            else
            {
                out_x[y][x] += in_y[y][x] * (in_x[next_y][x] - in_x[y][x]) / mGridSpacingY;
                out_y[y][x] += in_y[y][x] * (in_y[next_y][x] - in_y[y][x]) / mGridSpacingY;
            }

            prev_x = (prev_x + 1) % mNumGridPtsX;
            next_x = (next_x + 1) % mNumGridPtsX;
        }

        prev_y = (prev_y + 1) % mNumGridPtsY;
        next_y = (next_y + 1) % mNumGridPtsY;
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::Fft2DForwardRealToComplex(std::vector<std::vector<double> >& input, std::vector<std::vector<std::complex<double> > >& output)
{
    // Ensure output grid is set up correctly
    SetupGrid(output);

    // Set up some variables needed for fftw
    fftw_complex *complex_in_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * mNumGridPtsY * mNumGridPtsX);
    fftw_plan plan;

    // Rearrange input array into 1D vector of fftw_complexs
    unsigned count = 0;
    for(unsigned y = 0 ; y < mNumGridPtsY ; y++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            complex_in_out[count][0] = input[y][x];
            complex_in_out[count][1] = 0.0; // maybe don't need to do this
            count++;
        }
    }

    // Plan and perform in-place fft using fftw
    plan = fftw_plan_dft_2d(mNumGridPtsY, mNumGridPtsX, complex_in_out, complex_in_out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Rearrange 1D output into 2D vector of complex numbers, and normalise
    for( unsigned y = 0 ; y < mNumGridPtsY ; y++ )
    {
        for( unsigned x = 0 ; x < mNumGridPtsX ; x++ )
        {
            output[y][x] = (complex_in_out[y * mNumGridPtsX + x][0] + mI * complex_in_out[y * mNumGridPtsX + x][1]);
        }
    }

    // Free memory used for fft array
    fftw_free(complex_in_out);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::Fft2DInverseComplexToReal(std::vector<std::vector<std::complex<double> > >& input, std::vector<std::vector<double> >& output)
{
    // Ensure output grid is set up correctly
    SetupGrid(output);

    // Set up some variables needed for fftw
    fftw_complex *complex_in_out;
    complex_in_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * mNumGridPtsY * mNumGridPtsX);
    fftw_plan plan;

    // Rearrange input array into 1D vector of fftw_complexs
    unsigned count = 0;
    for(unsigned y = 0 ; y < mNumGridPtsY ; y++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            complex_in_out[count][0] = input[y][x].real();
            complex_in_out[count][1] = input[y][x].imag();
            count++;
        }
    }

    // Plan and perform in-place inverse FFT using FFTW
    plan = fftw_plan_dft_2d(mNumGridPtsY, mNumGridPtsX, complex_in_out, complex_in_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Rearrange 1D output into 2D vector of real doubles
    for( unsigned y = 0 ; y < mNumGridPtsY ; y++ )
    {
        for( unsigned x = 0 ; x < mNumGridPtsX ; x++ )
        {
            output[y][x] = complex_in_out[y * mNumGridPtsX + x][0] /(mFftNorm * mFftNorm);
        }
    }

    fftw_free(complex_in_out);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupGrid(std::vector<std::vector<std::complex<double> > >& grid)
{
    grid.resize(mNumGridPtsY);
    for (unsigned grid_y = 0; grid_y < mNumGridPtsY; grid_y++)
    {
        grid[grid_y].resize(mNumGridPtsX);
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupGrid(std::vector<std::vector<double> >& grid)
{
    grid.resize(mNumGridPtsY);
    for (unsigned grid_y = 0; grid_y < mNumGridPtsY; grid_y++)
    {
        grid[grid_y].resize(mNumGridPtsX);
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::PrintGrid(const std::vector<std::vector<double> >& grid)
{
    for (unsigned y = 0; y < mNumGridPtsY; y++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            std::cout << std::setprecision(5) << grid[y][x] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::PrintGrid(const std::vector<std::vector<std::complex<double>  > >& grid)
{
    for (unsigned y = 0; y < mNumGridPtsY; y++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            std::cout << std::setprecision(5) << grid[y][x].real() << "+" << grid[y][x].imag() << "i ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetMemberVariablesForTesting(unsigned numGridPtsY, unsigned numGridPtsX)
{
    mNumGridPtsY = numGridPtsY;
    mNumGridPtsX = numGridPtsX;
    mFftNorm = sqrt(numGridPtsX * numGridPtsY);

    mGridSpacingY = 1.0 / (double)mNumGridPtsY;
    mGridSpacingX = 1.0 / (double)mNumGridPtsX;
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetCellNeighbourUpdateFrequency(unsigned newFrequency)
{
    assert(newFrequency > 0);
    mCellNeighbourUpdateFrequency = newFrequency;
}

template<unsigned DIM>
unsigned ImmersedBoundarySimulationModifier<DIM>::GetCellNeighbourUpdateFrequency()
{
    return mCellNeighbourUpdateFrequency;
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::AddImmersedBoundaryForce(boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}


// Explicit instantiation
template class ImmersedBoundarySimulationModifier<1>;
template class ImmersedBoundarySimulationModifier<2>;
template class ImmersedBoundarySimulationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundarySimulationModifier)

