/*

Copyright (c) 2005-2017, University of Oxford.
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
#include "FluidSource.hpp"

template<unsigned DIM>
ImmersedBoundarySimulationModifier<DIM>::ImmersedBoundarySimulationModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mpMesh(NULL),
      mpCellPopulation(NULL),
      mNodeNeighbourUpdateFrequency(1u),
      mNumGridPtsX(0u),
      mNumGridPtsY(0u),
      mGridSpacingX(0.0),
      mGridSpacingY(0.0),
      mFftNorm(0.0),
      mpBoxCollection(NULL),
      mReynoldsNumber(1e-4),
      mI(0.0, 1.0),
      mpArrays(NULL),
      mpFftInterface(NULL)
{
}

template<unsigned DIM>
ImmersedBoundarySimulationModifier<DIM>::~ImmersedBoundarySimulationModifier()
{
    if (mpBoxCollection)
    {
        delete(mpBoxCollection);
    }
    if (mpArrays)
    {
        delete(mpArrays);
    }
    if (mpFftInterface)
    {
        delete(mpFftInterface);
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // We need to update node neighbours occasionally, but not necessarily each timestep
    if (SimulationTime::Instance()->GetTimeStepsElapsed() % mNodeNeighbourUpdateFrequency == 0)
    {
        mpBoxCollection->CalculateNodePairs(mpMesh->rGetNodes(), mNodePairs);
    }

    // This will solve the fluid problem for all timesteps after the first, which is handled in SetupSolve()
    this->UpdateFluidVelocityGrids(rCellPopulation);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // We can set up some helper variables here which need only be set up once for the entire simulation
    this->SetupConstantMemberVariables(rCellPopulation);

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
    this->ClearForcesAndSources();
    this->AddImmersedBoundaryForceContributions();
    this->PropagateForcesToFluidGrid();

    // If sources are active, we must propagate them from their nodes to the grid
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        this->PropagateFluidSourcesToGrid();
    }

    this->SolveNavierStokesSpectral();
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupConstantMemberVariables(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<ImmersedBoundaryCellPopulation<DIM> *>(&rCellPopulation) == NULL)
    {
        EXCEPTION("Cell population must be immersed boundary");
    }

    mpCellPopulation = static_cast<ImmersedBoundaryCellPopulation<DIM> *>(&rCellPopulation);
    mpMesh = &(mpCellPopulation->rGetMesh());

    // Get the size of the mesh
    mNumGridPtsX = mpMesh->GetNumGridPtsX();
    mNumGridPtsY = mpMesh->GetNumGridPtsY();

    // Get the grid spacing
    mGridSpacingX = 1.0 / (double) mNumGridPtsX;
    mGridSpacingY = 1.0 / (double) mNumGridPtsY;

    // Set up the box collection
    c_vector<double, 2 * 2> domain_size;
    domain_size(0) = 0.0;
    domain_size(1) = 1.0;
    domain_size(2) = 0.0;
    domain_size(3) = 1.0;
    mpBoxCollection = new ObsoleteBoxCollection<DIM>(mpCellPopulation->GetInteractionDistance(), domain_size, true, true);
    mpBoxCollection->SetupLocalBoxesHalfOnly();
    mpBoxCollection->CalculateNodePairs(mpMesh->rGetNodes(), mNodePairs);

    // Set up dimension-dependent variables
    switch (DIM)
    {
        case 2:
        {
            mpArrays = new ImmersedBoundary2dArrays<DIM>(mpMesh, SimulationTime::Instance()->GetTimeStep(), mReynoldsNumber, mpCellPopulation->DoesPopulationHaveActiveSources());
            mpFftInterface = new ImmersedBoundaryFftInterface<DIM>(mpMesh,
                                                                   &(mpArrays->rGetModifiableRightHandSideGrids()[0][0][0]),
                                                                   &(mpArrays->rGetModifiableFourierGrids()[0][0][0]),
                                                                   &(mpMesh->rGetModifiable2dVelocityGrids()[0][0][0]),
                                                                   mpCellPopulation->DoesPopulationHaveActiveSources());

            mpFftInterface_correction = new ImmersedBoundaryFftInterface<DIM>(mpMesh,
                                                                              &(mpArrays->rGetMeshGrid()[0][0][0]),
                                                                              &(mpArrays->rGetModifiablePressureCorrectionGrid()[0][0][0]),
                                                                              &(mpMesh->rGetModifiable2dVelocityGrids()[0][0][0]),
                                                                              mpCellPopulation->DoesPopulationHaveActiveSources());

            mpFftInterface_correction->FftExecuteCorrection();

            mFftNorm = (double) mNumGridPtsX * (double) mNumGridPtsY;
            break;
        }
        default:
            NEVER_REACHED; // Not yet implemented in 3D
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::ClearForcesAndSources()
{
    // Clear applied forces on each node
    for (typename ImmersedBoundaryMesh<DIM, DIM>::NodeIterator node_iter = mpMesh->GetNodeIteratorBegin(false);
         node_iter != mpMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

    // Reset force grids to 0 everywhere
    multi_array<double, 3>& r_force_grids = mpArrays->rGetModifiableForceGrids();

    for (unsigned dim = 0; dim < 2; dim++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                r_force_grids[dim][x][y] = 0.0;
            }
        }
    }

    // If there are active sources, the relevant grid needs to be reset to zero everywhere
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        multi_array<double, 3> &r_rhs_grid = mpArrays->rGetModifiableRightHandSideGrids();

        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                r_rhs_grid[2][x][y] = 0.0;
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::AddImmersedBoundaryForceContributions()
{
    // Add contributions from each immersed boundary force
    for (typename std::vector<boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        (*iter)->AddImmersedBoundaryForceContribution(mNodePairs, *mpCellPopulation);
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::PropagateForcesToFluidGrid()
{
    // Helper variables, pre-defined for efficiency
    double dl;
    double weight;

    unsigned first_idx_x;
    unsigned first_idx_y;

    std::vector<unsigned> x_indices(4);
    std::vector<unsigned> y_indices(4);

    std::vector<double> x_deltas(4);
    std::vector<double> y_deltas(4);

    c_vector<double, DIM> node_location;
    c_vector<double, DIM> applied_force;

    // Get a reference to the force grids which we spread the applied forces to
    multi_array<double, 3>& force_grids = mpArrays->rGetModifiableForceGrids();

    // Here, we loop over elements and then nodes as the length scale dl varies per element
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_iter = mpMesh->GetElementIteratorBegin(false);
         elem_iter != mpMesh->GetElementIteratorEnd();
         ++elem_iter)
    {
        dl = mpMesh->GetAverageNodeSpacingOfElement(elem_iter->GetIndex(), false);

        for (unsigned node_idx = 0; node_idx < elem_iter->GetNumNodes(); node_idx++)
        {
            Node<DIM> *p_node = elem_iter->GetNode(node_idx);

            // Get location and applied force contribution of current node
            node_location = p_node->rGetLocation();
            applied_force = p_node->rGetAppliedForce();

            // Get first grid index in each dimension, taking account of possible wrap-around
            first_idx_x = unsigned(floor(node_location[0] / mGridSpacingX)) + mNumGridPtsX - 1;
            first_idx_y = unsigned(floor(node_location[1] / mGridSpacingY)) + mNumGridPtsY - 1;

            // Calculate all four indices and deltas in each dimension
            for (unsigned i = 0; i < 4; i++)
            {
                x_indices[i] = (first_idx_x + i) % mNumGridPtsX;
                y_indices[i] = (first_idx_y + i) % mNumGridPtsY;

                x_deltas[i] = Delta1D(fabs(x_indices[i] * mGridSpacingX - node_location[0]), mGridSpacingX);
                y_deltas[i] = Delta1D(fabs(y_indices[i] * mGridSpacingY - node_location[1]), mGridSpacingY);
            }

            // Loop over the 4x4 grid used to spread the force on the nodes to the fluid grid
            for (unsigned x_idx = 0; x_idx < 4; x_idx++)
            {
                for (unsigned y_idx = 0; y_idx < 4; y_idx++)
                {
                    // The applied force is weighted by the delta function
                    weight = x_deltas[x_idx] * y_deltas[y_idx] * dl / (mGridSpacingX * mGridSpacingY);

                    force_grids[0][x_indices[x_idx]][y_indices[y_idx]] += applied_force[0] * weight;
                    force_grids[1][x_indices[x_idx]][y_indices[y_idx]] += applied_force[1] * weight;
                }
            }
        }
    }

    // Here, we loop over laminas and then nodes as the length scale dl varies per element
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryLaminaIterator lam_iter = mpMesh->GetLaminaIteratorBegin(false);
         lam_iter != mpMesh->GetLaminaIteratorEnd();
         ++lam_iter)
    {
        dl = mpMesh->GetAverageNodeSpacingOfLamina(lam_iter->GetIndex(), false);

        for (unsigned node_idx = 0; node_idx < lam_iter->GetNumNodes(); node_idx++)
        {
            Node<DIM> *p_node = lam_iter->GetNode(node_idx);

            // Get location and applied force contribution of current node
            node_location = p_node->rGetLocation();
            applied_force = p_node->rGetAppliedForce();

            // Get first grid index in each dimension, taking account of possible wrap-around
            first_idx_x = unsigned(floor(node_location[0] / mGridSpacingX)) + mNumGridPtsX - 1;
            first_idx_y = unsigned(floor(node_location[1] / mGridSpacingY)) + mNumGridPtsY - 1;

            // Calculate all four indices and deltas in each dimension
            for (unsigned i = 0; i < 4; i++)
            {
                x_indices[i] = (first_idx_x + i) % mNumGridPtsX;
                y_indices[i] = (first_idx_y + i) % mNumGridPtsY;

                x_deltas[i] = Delta1D(fabs(x_indices[i] * mGridSpacingX - node_location[0]), mGridSpacingX);
                y_deltas[i] = Delta1D(fabs(y_indices[i] * mGridSpacingY - node_location[1]), mGridSpacingY);
            }

            // Loop over the 4x4 grid used to spread the force on the nodes to the fluid grid
            for (unsigned x_idx = 0; x_idx < 4; x_idx++)
            {
                for (unsigned y_idx = 0; y_idx < 4; y_idx++)
                {
                    // The applied force is weighted by the delta function
                    weight = x_deltas[x_idx] * y_deltas[y_idx] * dl / (mGridSpacingX * mGridSpacingY);

                    force_grids[0][x_indices[x_idx]][y_indices[y_idx]] += applied_force[0] * weight;
                    force_grids[1][x_indices[x_idx]][y_indices[y_idx]] += applied_force[1] * weight;
                }
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::PropagateFluidSourcesToGrid()
{
    // Helper variables, pre-defined for efficiency
    double weight;

    unsigned first_idx_x;
    unsigned first_idx_y;

    std::vector<unsigned> x_indices(4);
    std::vector<unsigned> y_indices(4);

    std::vector<double> x_deltas(4);
    std::vector<double> y_deltas(4);

    // Currently the fluid source grid is the final part of the right hand side grid, as having all three grids
    // contiguous helps improve the Fourier transform performance.
    //\todo could make this nicer by using boost multiarray 'slice'?
    multi_array<double, 3>& rhs_grids = mpArrays->rGetModifiableRightHandSideGrids();

    std::vector<FluidSource<DIM>*>& r_element_sources = mpMesh->rGetElementFluidSources();
    std::vector<FluidSource<DIM>*>& r_balance_sources = mpMesh->rGetBalancingFluidSources();

    // Construct a vector of all sources combined
    std::vector<FluidSource<DIM>*> combined_sources;
    combined_sources.insert(combined_sources.end(), r_element_sources.begin(), r_element_sources.end());
    combined_sources.insert(combined_sources.end(), r_balance_sources.begin(), r_balance_sources.end());

    // Find the combined element source strength
    double cumulative_strength = 0.0;
    for (unsigned source_idx = 0; source_idx < r_element_sources.size(); source_idx++)
    {
        cumulative_strength += r_element_sources[source_idx]->GetStrength();
    }

    // Calculate the required balancing strength, and apply it to all balancing sources
    double balance_strength = -1.0 * cumulative_strength / (double)r_balance_sources.size();
    for (unsigned source_idx = 0; source_idx < r_balance_sources.size(); source_idx++)
    {
        r_balance_sources[source_idx]->SetStrength(balance_strength);
    }

    // Iterate over all sources and propagate their effects to the source grid
    for (unsigned source_idx = 0; source_idx < combined_sources.size(); source_idx++)
    {
        FluidSource<DIM>* this_source = combined_sources[source_idx];

        // Get location and strength of this source
        c_vector<double, DIM> source_location = this_source->rGetLocation();
        double source_strength = this_source->GetStrength();

        // Get first grid index in each dimension, taking account of possible wrap-around
        first_idx_x = unsigned(floor(source_location[0] / mGridSpacingX)) + mNumGridPtsX - 1;
        first_idx_y = unsigned(floor(source_location[1] / mGridSpacingY)) + mNumGridPtsY - 1;

        // Calculate all four indices and deltas in each dimension
        for (unsigned i = 0; i < 4; i ++)
        {
            x_indices[i] = (first_idx_x + i) % mNumGridPtsX;
            y_indices[i] = (first_idx_y + i) % mNumGridPtsY;

            x_deltas[i] = Delta1D(fabs(x_indices[i] * mGridSpacingX - source_location[0]), mGridSpacingX);
            y_deltas[i] = Delta1D(fabs(y_indices[i] * mGridSpacingY - source_location[1]), mGridSpacingY);
        }

        // Loop over the 4x4 grid needed to spread the source strength to the source grid
        for (unsigned x_idx = 0; x_idx < 4; x_idx ++)
        {
            for (unsigned y_idx = 0; y_idx < 4; y_idx ++)
            {
                // The strength is weighted by the delta function
                weight = x_deltas[x_idx] * y_deltas[y_idx] / (mGridSpacingX * mGridSpacingY);

                rhs_grids[2][x_indices[x_idx]][y_indices[y_idx]] += source_strength * weight;
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SolveNavierStokesSpectral()
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned reduced_size = 1 + (mNumGridPtsY/2);

    // Get references to all the necessary grids
    multi_array<double, 3>& vel_grids   = mpMesh->rGetModifiable2dVelocityGrids();
    multi_array<double, 3>& force_grids = mpArrays->rGetModifiableForceGrids();
    multi_array<double, 3>& rhs_grids   = mpArrays->rGetModifiableRightHandSideGrids();
    multi_array<double, 3>& source_gradient_grids   = mpArrays->rGetModifiableSourceGradientGrids();

    const multi_array<double, 2>& op_1  = mpArrays->rGetOperator1();
    const multi_array<double, 2>& op_2  = mpArrays->rGetOperator2();
    const std::vector<double>& sin_2x   = mpArrays->rGetSin2x();
    const std::vector<double>& sin_2y   = mpArrays->rGetSin2y();

    multi_array<std::complex<double>, 3>& fourier_grids = mpArrays->rGetModifiableFourierGrids();
    multi_array<std::complex<double>, 2>& pressure_grid = mpArrays->rGetModifiablePressureGrid();
    multi_array<std::complex<double>, 3>& correction_term_grid = mpArrays->rGetModifiablePressureCorrectionGrid(); // Correction term
    multi_array<double, 3>& acceleration_grid = mpArrays->rGetModifiableAccelerationGrid();

    // Perform upwind differencing and create RHS of linear system
    Upwind2d(vel_grids, rhs_grids);

    // If the population has active sources, the Right Hand Side grids are calculated differently.
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        double factor = 1.0 / (3.0 * mReynoldsNumber);

        CalculateSourceGradients(rhs_grids, source_gradient_grids);

        for (unsigned dim = 0; dim < 2; dim++)
        {
            for (unsigned x = 0; x < mNumGridPtsX; x++)
            {
                for (unsigned y = 0; y < mNumGridPtsY; y++)
                {
                    rhs_grids[dim][x][y] = vel_grids[dim][x][y] + dt * (force_grids[dim][x][y] + factor * source_gradient_grids[dim][x][y] - rhs_grids[dim][x][y]);
                }
            }
        }
    }
    else
    {
        for (unsigned dim = 0; dim < 2; dim++)
        {
            for (unsigned x = 0; x < mNumGridPtsX; x++)
            {
                for (unsigned y = 0; y < mNumGridPtsY; y++)
                {
                    rhs_grids[dim][x][y] = (vel_grids[dim][x][y] + dt * (force_grids[dim][x][y] - rhs_grids[dim][x][y]));
                }
            }
        }
    }

    CalculateCorrectionTerm(force_grids, acceleration_grid);

    // Perform fft on rhs_grids; results go to fourier_grids
    mpFftInterface->FftExecuteForward();

    /*
     * The result of a DFT of n real datapoints is n/2 + 1 complex values, due to redundancy: element n-1 is conj(2),
     * etc.  A similar pattern of redundancy occurs in a 2D transform.  Calculations in the Fourier domain preserve this
     * redundancy, and so all calculations need only be done on reduced-size arrays, saving memory and computation.
     */

    // If the population has active fluid sources, the computation is slightly more complicated
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < reduced_size; y++)
            {
                pressure_grid[x][y] = (op_2[x][y] * fourier_grids[2][x][y] - mI * (sin_2x[x] * fourier_grids[0][x][y] / mGridSpacingX +
                                                                                   sin_2y[y] * fourier_grids[1][x][y] / mGridSpacingY)) / op_1[x][y]
                                      - mDeltaPx * correction_term_grid[0][x][y] - mDeltaPy * correction_term_grid[1][x][y];
            }
        }
    }
    else // no active fluid sources
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < reduced_size; y++)
            {
                pressure_grid[x][y] = -mI * (sin_2x[x] * fourier_grids[0][x][y] / mGridSpacingX +
                                             sin_2y[y] * fourier_grids[1][x][y] / mGridSpacingY) / op_1[x][y]
                                      - mDeltaPx * correction_term_grid[0][x][y] - mDeltaPy * correction_term_grid[1][x][y];
            }
        }
    }


    // Set some values to zero
    pressure_grid[0][0] = 0.0;
    pressure_grid[mNumGridPtsX/2][0] = 0.0;
    pressure_grid[mNumGridPtsX/2][mNumGridPtsY/2] = 0.0;
    pressure_grid[0][mNumGridPtsY/2] = 0.0;

    /*
     * Do final stage of computation before inverse FFT.  We do the necessary DFT scaling at this stage so the output
     * from the inverse DFT is correct.
     */
    for (unsigned x = 0; x < mNumGridPtsX; x++)
    {
        for (unsigned y = 0; y < reduced_size; y++)
        {
            fourier_grids[0][x][y] = (fourier_grids[0][x][y] - (mI * dt / (mReynoldsNumber * mGridSpacingX)) * sin_2x[x] * pressure_grid[x][y]) / (op_2[x][y]);
            fourier_grids[1][x][y] = (fourier_grids[1][x][y] - (mI * dt / (mReynoldsNumber * mGridSpacingY)) * sin_2y[y] * pressure_grid[x][y]) / (op_2[x][y]);
        }
    }

    acceleration_grid = vel_grids;
    // Perform inverse fft on fourier_grids; results are in vel_grids.  Then, normalise the DFT.
    mpFftInterface->FftExecuteInverse();
    for (unsigned dim = 0; dim < 2; dim++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                vel_grids[dim][x][y] /= mFftNorm;
            }
        }
    }

    for (unsigned dim = 0; dim < 2; dim++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                acceleration_grid[dim][x][y] = (vel_grids[dim][x][y] - acceleration_grid[dim][x][y])/dt;
            }
        }
    }

}

template<unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::Delta1D(double dist, double spacing)
{
    return (0.25 * (1.0 + cos(M_PI * dist / (2 * spacing))));
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::Upwind2d(const multi_array<double, 3>& input, multi_array<double, 3>& output)
{
    unsigned prev_x = mNumGridPtsX - 1;
    unsigned prev_y = mNumGridPtsY - 1;

    unsigned next_x = 1;
    unsigned next_y = 1;

    for (unsigned x = 0; x < mNumGridPtsX; x++)
    {
        for (unsigned y = 0; y < mNumGridPtsY; y++)
        {
            // Set values for output from conditional on x grid
            if (input[0][x][y] > 0)
            {
                output[0][x][y] = input[0][x][y] * (input[0][x][y] - input[0][prev_x][y]) / mGridSpacingX;
                output[1][x][y] = input[0][x][y] * (input[1][x][y] - input[1][prev_x][y]) / mGridSpacingX;
            }
            else
            {
                output[0][x][y] = input[0][x][y] * (input[0][next_x][y] - input[0][x][y]) / mGridSpacingX;
                output[1][x][y] = input[0][x][y] * (input[1][next_x][y] - input[1][x][y]) / mGridSpacingX;
            }

            // Then add values from conditional on y grid
            if (input[1][x][y] > 0)
            {
                output[0][x][y] += input[1][x][y] * (input[0][x][y] - input[0][x][prev_y]) / mGridSpacingY;
                output[1][x][y] += input[1][x][y] * (input[1][x][y] - input[1][x][prev_y]) / mGridSpacingY;
            }
            else
            {
                output[0][x][y] += input[1][x][y] * (input[0][x][next_y] - input[0][x][y]) / mGridSpacingY;
                output[1][x][y] += input[1][x][y] * (input[1][x][next_y] - input[1][x][y]) / mGridSpacingY;
            }

            prev_y = (prev_y + 1) % mNumGridPtsY;
            next_y = (next_y + 1) % mNumGridPtsY;
        }

        prev_x = (prev_x + 1) % mNumGridPtsX;
        next_x = (next_x + 1) % mNumGridPtsX;
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::CalculateSourceGradients(const multi_array<double, 3>& rhs, multi_array<double, 3>& gradients)
{
    // \todo: this assumes mGridSpacingX = mGridSpacingY (fine for the foreseeable future)
    double factor = 1.0 / (2.0 * mGridSpacingX);

    // The fluid sources are stored in the third slice of the rhs grids
    for (unsigned x = 0; x < mNumGridPtsX; x++)
    {
        unsigned next_x = (x + 1) % mNumGridPtsX;
        unsigned prev_x = (x + mNumGridPtsX - 1) % mNumGridPtsX;

        for (unsigned y = 0; y < mNumGridPtsY; y++)
        {
            unsigned next_y = (y + 1) % mNumGridPtsY;
            unsigned prev_y = (y + mNumGridPtsY - 1) % mNumGridPtsY;

            gradients[0][x][y] = factor * (rhs[2][next_x][y] - rhs[2][prev_x][y]);
            gradients[1][x][y] = factor * (rhs[2][x][next_y] - rhs[2][x][prev_y]);
        }
    }
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
void ImmersedBoundarySimulationModifier<DIM>::SetNodeNeighbourUpdateFrequency(unsigned newFrequency)
{
    assert(newFrequency > 0);
    mNodeNeighbourUpdateFrequency = newFrequency;
}

template<unsigned DIM>
unsigned ImmersedBoundarySimulationModifier<DIM>::GetNodeNeighbourUpdateFrequency()
{
    return mNodeNeighbourUpdateFrequency;
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::AddImmersedBoundaryForce(boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetReynoldsNumber(double reynoldsNumber)
{
    assert(reynoldsNumber > 0.0);
    mReynoldsNumber = reynoldsNumber;
}

template<unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::GetReynoldsNumber()
{
    return mReynoldsNumber;
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::CalculateCorrectionTerm(const multi_array<double, 3>& force_grids, const multi_array<double, 3>& acceleration_grids)
{
    // extra input -> , multi_array<double, 3> &output
    // input array should to be the the force grid

    // Trapezium rule is used to approximate the integral of the force on the fluid over the whole domain


    multi_array<double, 3> input;
    input.resize(extents[2][mNumGridPtsX][mNumGridPtsY]);

    for (unsigned dim = 0; dim < 2; dim++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                input[dim][x][y] = acceleration_grids[dim][x][y] - force_grids[dim][x][y];
            }
        }
    }

    double delta_p_x = 1.0 * (input[0][0][0] + input[0][mNumGridPtsX-1][0] + input[0][0][mNumGridPtsY-1] + input[0][mNumGridPtsX-1][mNumGridPtsY-1]);
    for (unsigned i=1; i<mNumGridPtsX-2; i++)
    {
        delta_p_x += 2.0 * input[0][i][0];
        delta_p_x += 2.0 * input[0][i][mNumGridPtsX-1];
    }
    for (unsigned i=1; i<mNumGridPtsY-2; i++)
    {
        delta_p_x += 2.0 * input[0][0][i];
        delta_p_x += 2.0 * input[0][mNumGridPtsY-1][i];
    }
    for (unsigned i=1; i<mNumGridPtsX-2; i++)
    {
        for (unsigned j=1; j<mNumGridPtsY-2; j++)
        {
            delta_p_x += 4.0 * input[0][i][j];
        }
    }

    double delta_p_y = 1.0 * (input[1][0][0] + input[1][mNumGridPtsX-1][0] + input[1][0][mNumGridPtsY-1] + input[1][mNumGridPtsX-1][mNumGridPtsY-1]);
    for (unsigned i=1; i<mNumGridPtsX-2; i++)
    {
        delta_p_x += 2.0 * input[1][i][0];
        delta_p_x += 2.0 * input[1][i][mNumGridPtsX-1];
    }
    for (unsigned i=1; i<mNumGridPtsY-2; i++)
    {
        delta_p_x += 2.0 * input[1][0][i];
        delta_p_x += 2.0 * input[1][mNumGridPtsY-1][i];
    }
    for (unsigned i=1; i<mNumGridPtsX-2; i++)
    {
        for (unsigned j=1; j<mNumGridPtsY-2; j++)
        {
            delta_p_x += 4.0 * input[1][i][j];
        }
    }

    mDeltaPx = -0.25 * mReynoldsNumber * mGridSpacingX * mGridSpacingY * delta_p_x;
    mDeltaPy = -0.25 * mReynoldsNumber * mGridSpacingX * mGridSpacingY * delta_p_y;
}

// Explicit instantiation
template class ImmersedBoundarySimulationModifier<1>;
template class ImmersedBoundarySimulationModifier<2>;
template class ImmersedBoundarySimulationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundarySimulationModifier)

