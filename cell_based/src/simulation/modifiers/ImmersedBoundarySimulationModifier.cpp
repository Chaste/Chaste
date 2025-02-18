/*

Copyright (c) 2005-2025, University of Oxford.
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

#include <memory>

#include "FluidSource.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

template<unsigned DIM>
ImmersedBoundarySimulationModifier<DIM>::ImmersedBoundarySimulationModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mpMesh(nullptr),
      mpCellPopulation(nullptr),
      mNodeNeighbourUpdateFrequency(1u),
      mGridSpacingX(0.0),
      mGridSpacingY(0.0),
      mFftNorm(0.0),
      mAdditiveNormalNoise(false),
      mNoiseStrength(0.0),
      mNoiseSkip(1u),
      mNoiseLengthScale(0.1),
      mZeroFieldSums(false),
      mpBoxCollection(nullptr),
      mReynoldsNumber(1e-4),
      mI(0.0, 1.0),
      mpArrays(nullptr),
      mpFftInterface(nullptr),
      mpRandomField(nullptr)
{
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateAtEndOfTimeStep(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
void ImmersedBoundarySimulationModifier<DIM>::SetupSolve(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation,
    std::string outputDirectory)
{
    /*
     * We can set up some helper variables here which need only be set up once
     * for the entire simulation.
     */
    this->SetupConstantMemberVariables(rCellPopulation);

    // This will solve the fluid problem based on the initial mesh setup
    this->UpdateFluidVelocityGrids(rCellPopulation);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::OutputSimulationModifierParameters(
    out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateFluidVelocityGrids(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    this->ClearForcesAndSources();
    this->RecalculateAverageNodeSpacings();
    this->AddImmersedBoundaryForceContributions();
    this->PropagateForcesToFluidGrid();

    // If random noise is required, add it to the force grids
    if (mAdditiveNormalNoise)
    {
        this->AddNormalNoise();
    }

    // If sources are active, we must propagate them from their nodes to the grid
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        this->PropagateFluidSourcesToGrid();
    }

    this->SolveNavierStokesSpectral();
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetupConstantMemberVariables(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (dynamic_cast<ImmersedBoundaryCellPopulation<DIM> *>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("Cell population must be immersed boundary");
    }

    mpCellPopulation = static_cast<ImmersedBoundaryCellPopulation<DIM> *>(&rCellPopulation);
    mpMesh = &(mpCellPopulation->rGetMesh());

    // Get the size of the mesh
    unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
    unsigned num_grid_ptsY = mpMesh->GetNumGridPtsY();

    // Get the grid spacing
    mGridSpacingX = 1.0 / (double) num_grid_ptsX;
    mGridSpacingY = 1.0 / (double) num_grid_ptsY;

    // Set up the box collection
    c_vector<double, 2 * 2> domain_size;
    domain_size(0) = 0.0;
    domain_size(1) = 1.0;
    domain_size(2) = 0.0;
    domain_size(3) = 1.0;

    mpBoxCollection = std::make_unique<ObsoleteBoxCollection<DIM>>(
            mpCellPopulation->GetInteractionDistance(),
            domain_size,
            true, // LCOV_EXCL_LINE
            true
    );
    mpBoxCollection->SetupLocalBoxesHalfOnly();
    mpBoxCollection->CalculateNodePairs(mpMesh->rGetNodes(), mNodePairs);

    // Set up dimension-dependent variables
    switch (DIM)
    {
        case 2:
        {
            mpArrays = std::make_unique<ImmersedBoundary2dArrays<DIM>>(
                    mpMesh,
                    SimulationTime::Instance()->GetTimeStep(),
                    mReynoldsNumber,
                    mpCellPopulation->DoesPopulationHaveActiveSources()
            );

            mpFftInterface = std::make_unique<ImmersedBoundaryFftInterface<DIM>>(
                    mpMesh,
                    &(mpArrays->rGetModifiableRightHandSideGrids()[0][0][0]),
                    &(mpArrays->rGetModifiableFourierGrids()[0][0][0]),
                    &(mpMesh->rGetModifiable2dVelocityGrids()[0][0][0]),
                    mpCellPopulation->DoesPopulationHaveActiveSources()
            );

            mFftNorm = (double) num_grid_ptsX * (double) num_grid_ptsY;
            break;
        }
        default:
            NEVER_REACHED; // Not yet implemented in 3D
    }

    // Set up random noise, if required
    if(mAdditiveNormalNoise)
    {
        std::array<double, DIM> lower_corner;
        lower_corner.fill(0.0);

        std::array<double, DIM> upper_corner;
        upper_corner.fill(1.0);

        std::array<unsigned, DIM> num_grid_pts;
        num_grid_pts.fill(num_grid_ptsX / mNoiseSkip);

        std::array<bool, DIM> periodicity;
        periodicity.fill(true);

        const double total_gridpts = std::accumulate(num_grid_pts.begin(),
                                                     num_grid_pts.end(),
                                                     1.0,
                                                     std::multiplies<double>());

        // Warn at this point if parameters are not sensible
        if (num_grid_ptsX % mNoiseSkip != 0)
        {
            WARNING("mNoiseSkip should perfectly divide num_grid_ptsX or adding forces will likely not work as expected.");
        }
        if (total_gridpts > 100.0 * 100.0)  // This is about the maximum that can be handled
        {
            WARNING("You probably have too many grid points for creating a random field.  Increase mNoiseSkip");
        }

        // Set up the random field generator and save it to cache
        mpRandomField = std::make_unique<UniformGridRandomFieldGenerator<DIM>>(
                lower_corner,
                upper_corner,
                num_grid_pts,
                periodicity,
                mNoiseLengthScale
        );

    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::ClearForcesAndSources()
{
    unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
    unsigned num_grid_ptsY = mpMesh->GetNumGridPtsY();

    // Clear applied forces on each node
    for (auto node_iter = mpMesh->GetNodeIteratorBegin(false);
         node_iter != mpMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

    // Reset force grids to 0 everywhere
    multi_array<double, 3>& r_force_grids = mpArrays->rGetModifiableForceGrids();

    for (unsigned dim = 0; dim < 2; dim++)
    {
        for (unsigned x = 0; x < num_grid_ptsX; x++)
        {
            for (unsigned y = 0; y < num_grid_ptsY; y++)
            {
                r_force_grids[dim][x][y] = 0.0;
            }
        }
    }

    // If there are active sources, the relevant grid needs to be reset to zero everywhere
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        multi_array<double, 3> &r_rhs_grid = mpArrays->rGetModifiableRightHandSideGrids();

        for (unsigned x = 0; x < num_grid_ptsX; x++)
        {
            for (unsigned y = 0; y < num_grid_ptsY; y++)
            {
                r_rhs_grid[2][x][y] = 0.0;
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::RecalculateAverageNodeSpacings()
{
    for (auto elem_it = mpMesh->GetElementIteratorBegin(false);
         elem_it != mpMesh->GetElementIteratorEnd();
         ++elem_it)
    {
        mpMesh->GetAverageNodeSpacingOfElement(elem_it->GetIndex(), true);
    }

    for (auto lam_it = mpMesh->GetLaminaIteratorBegin(false);
         lam_it != mpMesh->GetLaminaIteratorEnd();
         ++lam_it)
    {
        mpMesh->GetAverageNodeSpacingOfLamina(lam_it->GetIndex(), true);
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::AddImmersedBoundaryForceContributions()
{
    // Add contributions from each immersed boundary force
    for (auto iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        (*iter)->AddImmersedBoundaryForceContribution(mNodePairs, *mpCellPopulation);
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::PropagateForcesToFluidGrid()
{
    if constexpr (DIM == 2)
    {
        unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
        unsigned num_grid_ptsY = mpMesh->GetNumGridPtsY();

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
        multi_array<double, 3>& r_force_grids = mpArrays->rGetModifiableForceGrids();

        // Here, we loop over elements and then nodes as the length scale dl varies per element
        for (auto elem_iter = mpMesh->GetElementIteratorBegin(false);
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
                first_idx_x = unsigned(floor(node_location[0] / mGridSpacingX)) + num_grid_ptsX - 1;
                first_idx_y = unsigned(floor(node_location[1] / mGridSpacingY)) + num_grid_ptsY - 1;

                // Calculate all four indices and deltas in each dimension
                for (unsigned i = 0; i < 4; i++)
                {
                    x_indices[i] = (first_idx_x + i) % num_grid_ptsX;
                    y_indices[i] = (first_idx_y + i) % num_grid_ptsY;

                    x_deltas[i] = Delta1D(fabs(x_indices[i] * mGridSpacingX - node_location[0]), mGridSpacingX);
                    y_deltas[i] = Delta1D(fabs(y_indices[i] * mGridSpacingY - node_location[1]), mGridSpacingY);
                }

                // Loop over the 4x4 grid used to spread the force on the nodes to the fluid grid
                for (unsigned x_idx = 0; x_idx < 4; ++x_idx)
                {
                    for (unsigned y_idx = 0; y_idx < 4; ++y_idx)
                    {
                        // The applied force is weighted by the delta function
                        weight = x_deltas[x_idx] * y_deltas[y_idx] * dl / (mGridSpacingX * mGridSpacingY);

                        r_force_grids[0][x_indices[x_idx]][y_indices[y_idx]] += applied_force[0] * weight;
                        r_force_grids[1][x_indices[x_idx]][y_indices[y_idx]] += applied_force[1] * weight;
                    }
                }
            }
        }

        // Here, we loop over laminas and then nodes as the length scale dl varies per element
        for (auto lam_iter = mpMesh->GetLaminaIteratorBegin(false);
             lam_iter != mpMesh->GetLaminaIteratorEnd();
             ++lam_iter)
        {
            dl = mpMesh->GetAverageNodeSpacingOfLamina(lam_iter->GetIndex(), false);

            for (unsigned node_idx = 0; node_idx < lam_iter->GetNumNodes(); ++node_idx)
            {
                Node<DIM> *p_node = lam_iter->GetNode(node_idx);

                // Get location and applied force contribution of current node
                node_location = p_node->rGetLocation();
                applied_force = p_node->rGetAppliedForce();

                // Get first grid index in each dimension, taking account of possible wrap-around
                first_idx_x = unsigned(floor(node_location[0] / mGridSpacingX)) + num_grid_ptsX - 1;
                first_idx_y = unsigned(floor(node_location[1] / mGridSpacingY)) + num_grid_ptsY - 1;

                // Calculate all four indices and deltas in each dimension
                for (unsigned i = 0; i < 4; i++)
                {
                    x_indices[i] = (first_idx_x + i) % num_grid_ptsX;
                    y_indices[i] = (first_idx_y + i) % num_grid_ptsY;

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

                        r_force_grids[0][x_indices[x_idx]][y_indices[y_idx]] += applied_force[0] * weight;
                        r_force_grids[1][x_indices[x_idx]][y_indices[y_idx]] += applied_force[1] * weight;
                    }
                }
            }
        }

        // Finally, zero out any systematic small errors on the force grids that may lead to a drift, if required
        if (mZeroFieldSums)
        {
            ZeroFieldSums(r_force_grids);
        }
    }
    else
    {
        NEVER_REACHED;
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::PropagateFluidSourcesToGrid()
{
    if constexpr (DIM == 2)
    {
        unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
        unsigned num_grid_ptsY = mpMesh->GetNumGridPtsY();

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
        multi_array<double, 3>& r_rhs_grids = mpArrays->rGetModifiableRightHandSideGrids();

        std::vector<std::shared_ptr<FluidSource<DIM>>>& r_element_sources = mpMesh->rGetElementFluidSources();
        std::vector<std::shared_ptr<FluidSource<DIM>>>& r_balance_sources = mpMesh->rGetBalancingFluidSources();

        // Construct a vector of all sources combined
        std::vector<std::shared_ptr<FluidSource<DIM>>> combined_sources;
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
            std::shared_ptr<FluidSource<DIM>> this_source = combined_sources[source_idx];

            // Get location and strength of this source
            c_vector<double, DIM> source_location = this_source->rGetLocation();
            double source_strength = this_source->GetStrength();

            // Get first grid index in each dimension, taking account of possible wrap-around
            first_idx_x = unsigned(floor(source_location[0] / mGridSpacingX)) + num_grid_ptsX - 1;
            first_idx_y = unsigned(floor(source_location[1] / mGridSpacingY)) + num_grid_ptsY - 1;

            // Calculate all four indices and deltas in each dimension
            for (unsigned i = 0; i < 4; i ++)
            {
                x_indices[i] = (first_idx_x + i) % num_grid_ptsX;
                y_indices[i] = (first_idx_y + i) % num_grid_ptsY;

                x_deltas[i] = Delta1D(fabs(x_indices[i] * mGridSpacingX - source_location[0]), mGridSpacingX);
                y_deltas[i] = Delta1D(fabs(y_indices[i] * mGridSpacingY - source_location[1]), mGridSpacingY);
            }

            // Loop over the 4x4 grid needed to spread the source strength to the source grid
            for (unsigned x_idx = 0; x_idx < 4; ++x_idx)
            {
                for (unsigned y_idx = 0; y_idx < 4; ++y_idx)
                {
                    // The strength is weighted by the delta function
                    weight = x_deltas[x_idx] * y_deltas[y_idx] / (mGridSpacingX * mGridSpacingY);

                    r_rhs_grids[2][x_indices[x_idx]][y_indices[y_idx]] += source_strength * weight;
                }
            }
        }
    }
    else
    {
        NEVER_REACHED;
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SolveNavierStokesSpectral()
{
    unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
    unsigned num_grid_ptsY = mpMesh->GetNumGridPtsY();

    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned reduced_size = 1 + (num_grid_ptsY/2);

    // Get references to all the necessary grids
    multi_array<double, 3>& r_vel_grids = mpMesh->rGetModifiable2dVelocityGrids();
    multi_array<double, 3>& r_force_grids = mpArrays->rGetModifiableForceGrids();
    multi_array<double, 3>& r_rhs_grids = mpArrays->rGetModifiableRightHandSideGrids();
    multi_array<double, 3>& r_source_gradient_grids = mpArrays->rGetModifiableSourceGradientGrids();

    const multi_array<double, 2>& r_op_1 = mpArrays->rGetOperator1();
    const multi_array<double, 2>& r_op_2 = mpArrays->rGetOperator2();
    const std::vector<double>& r_sin_2x = mpArrays->rGetSin2x();
    const std::vector<double>& r_sin_2y = mpArrays->rGetSin2y();

    multi_array<std::complex<double>, 3>& r_fourier_grids = mpArrays->rGetModifiableFourierGrids();
    multi_array<std::complex<double>, 2>& r_pressure_grid = mpArrays->rGetModifiablePressureGrid();

    // Perform upwind differencing and create RHS of linear system
    Upwind2d(r_vel_grids, r_rhs_grids);

    // If the population has active sources, the Right Hand Side grids are calculated differently.
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        double factor = 1.0 / (3.0 * mReynoldsNumber);

        CalculateSourceGradients(r_rhs_grids, r_source_gradient_grids);

        for (unsigned dim = 0; dim < 2; ++dim)
        {
            for (unsigned x = 0; x < num_grid_ptsX; ++x)
            {
                for (unsigned y = 0; y < num_grid_ptsY; ++y)
                {
                    r_rhs_grids[dim][x][y] = r_vel_grids[dim][x][y] + dt * (r_force_grids[dim][x][y] + factor * r_source_gradient_grids[dim][x][y] - r_rhs_grids[dim][x][y]);
                }
            }
        }
    }
    else
    {
        for (unsigned dim = 0; dim < 2; ++dim)
        {
            for (unsigned x = 0; x < num_grid_ptsX; ++x)
            {
                for (unsigned y = 0; y < num_grid_ptsY; ++y)
                {
                    r_rhs_grids[dim][x][y] = (r_vel_grids[dim][x][y] + dt * (r_force_grids[dim][x][y] - r_rhs_grids[dim][x][y]));
                }
            }
        }
    }

    // Perform fft on r_rhs_grids; results go to r_fourier_grids
    mpFftInterface->FftExecuteForward();

    /*
     * The result of a DFT of n real datapoints is n/2 + 1 complex values, due
     * to redundancy: element n-1 is conj(2), etc. A similar pattern of
     * redundancy occurs in a 2D transform. Calculations in the Fourier domain
     * preserve this redundancy, and so all calculations need only be done on
     * reduced-size arrays, saving memory and computation.
     */

    // If the population has active fluid sources, the computation is slightly more complicated
    if (mpCellPopulation->DoesPopulationHaveActiveSources())
    {
        for (unsigned x = 0; x < num_grid_ptsX; ++x)
        {
            for (unsigned y = 0; y < reduced_size; ++y)
            {
                r_pressure_grid[x][y] = (r_op_2[x][y] * r_fourier_grids[2][x][y] - mI * (r_sin_2x[x] * r_fourier_grids[0][x][y] / mGridSpacingX +
                                                                                   r_sin_2y[y] * r_fourier_grids[1][x][y] / mGridSpacingY)) / r_op_1[x][y];
            }
        }
    }
    else // no active fluid sources
    {
        for (unsigned x = 0; x < num_grid_ptsX; ++x)
        {
            for (unsigned y = 0; y < reduced_size; ++y)
            {
                r_pressure_grid[x][y] = -mI * (r_sin_2x[x] * r_fourier_grids[0][x][y] / mGridSpacingX +
                                             r_sin_2y[y] * r_fourier_grids[1][x][y] / mGridSpacingY) / r_op_1[x][y];
            }
        }
    }

    // Set some values to zero
    r_pressure_grid[0][0] = 0.0;
    r_pressure_grid[num_grid_ptsX/2][0] = 0.0;
    r_pressure_grid[num_grid_ptsX/2][num_grid_ptsY/2] = 0.0;
    r_pressure_grid[0][num_grid_ptsY/2] = 0.0;

    /*
     * Do final stage of computation before inverse FFT. We do the necessary DFT
     * scaling at this stage so the output from the inverse DFT is correct.
     */
    for (unsigned x = 0; x < num_grid_ptsX; ++x)
    {
        for (unsigned y = 0; y < reduced_size; ++y)
        {
            r_fourier_grids[0][x][y] = (r_fourier_grids[0][x][y] - (mI * dt / (mReynoldsNumber * mGridSpacingX)) * r_sin_2x[x] * r_pressure_grid[x][y]) / (r_op_2[x][y]);
            r_fourier_grids[1][x][y] = (r_fourier_grids[1][x][y] - (mI * dt / (mReynoldsNumber * mGridSpacingY)) * r_sin_2y[y] * r_pressure_grid[x][y]) / (r_op_2[x][y]);
        }
    }

    // Perform inverse fft on r_fourier_grids; results are in r_vel_grids. Then, normalise the DFT.
    mpFftInterface->FftExecuteInverse();
    for (unsigned dim = 0; dim < 2; ++dim)
    {
        for (unsigned x = 0; x < num_grid_ptsX; ++x)
        {
            for (unsigned y = 0; y < num_grid_ptsY; ++y)
            {
                r_vel_grids[dim][x][y] /= mFftNorm;
            }
        }
    }

    /*
     * Finally, zero out any systematic small errors on the velocity grids that
     * may lead to a drift, if required.
     */
    if (mZeroFieldSums)
    {
        ZeroFieldSums(r_vel_grids);
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::ZeroFieldSums(
    multi_array<double, 3>& rField)
{
    unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
    unsigned num_grid_ptsY = mpMesh->GetNumGridPtsY();

    for (unsigned dim = 0; dim < 2; ++dim)
    {
        double total_field_val = 0.0;
        long count = 0;

        for (unsigned x = 0; x < num_grid_ptsX; ++x)
        {
            for (unsigned y = 0; y < num_grid_ptsY; ++y)
            {
                if (rField[dim][x][y] != 0.0)
                {
                    total_field_val += rField[dim][x][y];
                    count++;
                }
            }
        }

        const double average_field_val = (total_field_val / count);

        for (unsigned x = 0; x < num_grid_ptsX; ++x)
        {
            for (unsigned y = 0; y < num_grid_ptsY; ++y)
            {
                if (rField[dim][x][y] != 0.0)
                {
                    rField[dim][x][y] -= average_field_val;
                }
            }
        }
    }
}

template<unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::Delta1D(
    double dist,
    double spacing)
{
    return (0.25 * (1.0 + cos(M_PI * dist / (2 * spacing))));
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::Upwind2d(
    const multi_array<double, 3>& rInput,
    multi_array<double, 3>& rOutput)
{
    unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
    unsigned num_grid_ptsY = mpMesh->GetNumGridPtsX();

    unsigned prev_x = num_grid_ptsX - 1;
    unsigned prev_y = num_grid_ptsY - 1;

    unsigned next_x = 1;
    unsigned next_y = 1;

    for (unsigned x = 0; x < num_grid_ptsX; x++)
    {
        for (unsigned y = 0; y < num_grid_ptsY; y++)
        {
            // Set values for rOutput from conditional on x grid
            if (rInput[0][x][y] > 0)
            {
                rOutput[0][x][y] = rInput[0][x][y] * (rInput[0][x][y] - rInput[0][prev_x][y]) / mGridSpacingX;
                rOutput[1][x][y] = rInput[0][x][y] * (rInput[1][x][y] - rInput[1][prev_x][y]) / mGridSpacingX;
            }
            else
            {
                rOutput[0][x][y] = rInput[0][x][y] * (rInput[0][next_x][y] - rInput[0][x][y]) / mGridSpacingX;
                rOutput[1][x][y] = rInput[0][x][y] * (rInput[1][next_x][y] - rInput[1][x][y]) / mGridSpacingX;
            }

            // Then add values from conditional on y grid
            if (rInput[1][x][y] > 0)
            {
                rOutput[0][x][y] += rInput[1][x][y] * (rInput[0][x][y] - rInput[0][x][prev_y]) / mGridSpacingY;
                rOutput[1][x][y] += rInput[1][x][y] * (rInput[1][x][y] - rInput[1][x][prev_y]) / mGridSpacingY;
            }
            else
            {
                rOutput[0][x][y] += rInput[1][x][y] * (rInput[0][x][next_y] - rInput[0][x][y]) / mGridSpacingY;
                rOutput[1][x][y] += rInput[1][x][y] * (rInput[1][x][next_y] - rInput[1][x][y]) / mGridSpacingY;
            }

            prev_y = (prev_y + 1) % num_grid_ptsY;
            next_y = (next_y + 1) % num_grid_ptsY;
        }

        prev_x = (prev_x + 1) % num_grid_ptsX;
        next_x = (next_x + 1) % num_grid_ptsX;
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::CalculateSourceGradients(
    const multi_array<double, 3>& rRhs,
    multi_array<double, 3>& rGradients)
{
    unsigned num_grid_ptsX = mpMesh->GetNumGridPtsX();
    unsigned num_grid_ptsY = mpMesh->GetNumGridPtsY();

    // Assumes 1.0 x 1.0 square domain
    double factor = 1.0 / (2.0 * mGridSpacingX);

    // The fluid sources are stored in the third slice of the rRhs grids
    for (unsigned x = 0; x < num_grid_ptsX; ++x)
    {
        unsigned next_x = (x + 1) % num_grid_ptsX;
        unsigned prev_x = (x + num_grid_ptsX - 1) % num_grid_ptsX;

        for (unsigned y = 0; y < num_grid_ptsY; ++y)
        {
            unsigned next_y = (y + 1) % num_grid_ptsY;
            unsigned prev_y = (y + num_grid_ptsY - 1) % num_grid_ptsY;

            rGradients[0][x][y] = factor * (rRhs[2][next_x][y] - rRhs[2][prev_x][y]);
            rGradients[1][x][y] = factor * (rRhs[2][x][next_y] - rRhs[2][x][prev_y]);
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetNodeNeighbourUpdateFrequency(
    unsigned newFrequency)
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
void ImmersedBoundarySimulationModifier<DIM>::AddImmersedBoundaryForce(
    boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::AddNormalNoise() const
{
    auto& r_force_grids = mpArrays->rGetModifiableForceGrids();

    // This scaling is to ensure "correct" scaling with timestep
    const double force_scale_factor = std::sqrt(2.0 * mNoiseStrength / SimulationTime::Instance()->GetTimeStep());

    // Add the noise
    for (unsigned dim = 0; dim < r_force_grids.shape()[0]; dim++)
    {
        // Get an instance of the random field for the current dimension
        const std::vector<double> field = mpRandomField->SampleRandomField();

        // Calculate the sum of the random field, which we must adjust to have a net-zero impact
        const double adjustment = std::accumulate(field.begin(), field.end(), 0.0) / field.size();

        for (unsigned x = 0, idx = 0; x < r_force_grids.shape()[1]; x += mNoiseSkip)
        {
            for (unsigned y = 0; y < r_force_grids.shape()[2]; y += mNoiseSkip, idx++)
            {
                // Value from the random grid at this location
                const double val = force_scale_factor * (field[idx] - adjustment);

                // Fill in the mSkip x mSkip square of points
                for (unsigned i = 0; i < mNoiseSkip; ++i)
                {
                    for (unsigned j = 0; j < mNoiseSkip; ++j)
                    {
                        r_force_grids[dim][x+i][y+j] += val;
                    }
                }
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetReynoldsNumber(
    double reynoldsNumber)
{
    assert(reynoldsNumber > 0.0);
    mReynoldsNumber = reynoldsNumber;
}

template<unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::GetReynoldsNumber()
{
    return mReynoldsNumber;
}

template <unsigned DIM>
bool ImmersedBoundarySimulationModifier<DIM>::GetAdditiveNormalNoise() const
{
    return mAdditiveNormalNoise;
}

template <unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetAdditiveNormalNoise(
    bool additiveNormalNoise)
{
    mAdditiveNormalNoise = additiveNormalNoise;
}

template <unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::GetNoiseStrength() const
{
    return mNoiseStrength;
}

template <unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetNoiseStrength(
    double noiseStrength)
{
    mNoiseStrength = noiseStrength;
}

template <unsigned DIM>
unsigned ImmersedBoundarySimulationModifier<DIM>::GetNoiseSkip() const
{
    return mNoiseSkip;
}

template <unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetNoiseSkip(unsigned noiseSkip)
{
    mNoiseSkip = noiseSkip;
}

template <unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::GetNoiseLengthScale() const
{
    return mNoiseLengthScale;
}

template <unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetNoiseLengthScale(
    double noiseLengthScale)
{
    mNoiseLengthScale = noiseLengthScale;
}

template<unsigned int DIM>
bool ImmersedBoundarySimulationModifier<DIM>::GetZeroFieldSums() const
{
    return mZeroFieldSums;
}

template<unsigned int DIM>
void ImmersedBoundarySimulationModifier<DIM>::SetZeroFieldSums(
    bool zeroFieldSums)
{
    mZeroFieldSums = zeroFieldSums;
}

// Explicit instantiation
template class ImmersedBoundarySimulationModifier<1>;
template class ImmersedBoundarySimulationModifier<2>;
template class ImmersedBoundarySimulationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundarySimulationModifier)

