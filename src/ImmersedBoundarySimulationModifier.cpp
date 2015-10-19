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
#include "Warnings.hpp"
#include <complex>
#include <fftw3.h>
#include "Debug.hpp"
#include "Timer.hpp"

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
      mReynolds(1e-4),
      mI(0.0 + 1.0 * 1i),
      mpArrays(NULL)
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

    fftw_destroy_plan(mFftwForwardPlanX);
    fftw_destroy_plan(mFftwForwardPlanY);

    fftw_destroy_plan(mFftwInversePlanX);
    fftw_destroy_plan(mFftwInversePlanY);
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // We need to update node neighbours occasionally, but not necessarily each timestep
    if(SimulationTime::Instance()->GetTimeStepsElapsed() % mNodeNeighbourUpdateFrequency == 0)
    {
        mpBoxCollection->CalculateNodePairs(mpMesh->rGetNodes(), mNodePairs, mNodeNeighbours);
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
    this->ClearForces();
    this->AddForceContributions();
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

    // Set up dimension-dependent variables
    switch (DIM)
    {
        case 2:
            mpArrays = new ImmersedBoundary2dArrays(mNumGridPtsX, mNumGridPtsY);

            mFftNorm = sqrt((double) mNumGridPtsX * (double) mNumGridPtsY);
            break;

        case 3:
            EXCEPTION("Not implemented yet in 3D");
            break;

        default:
            NEVER_REACHED;
    }

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

    // Set up the box collection
    c_vector<double, 2 * 2> domain_size;
    domain_size(0) = 0.0;
    domain_size(1) = 1.0;
    domain_size(2) = 0.0;
    domain_size(3) = 1.0;
    mpBoxCollection = new BoxCollection<DIM>(mpCellPopulation->GetInteractionDistance(), domain_size, true, true);
    mpBoxCollection->SetupLocalBoxesHalfOnly();
    mpBoxCollection->CalculateNodePairs(mpMesh->rGetNodes(), mNodePairs, mNodeNeighbours);

    std::string filename = "./projects/ImmersedBoundary/src/fftw.wisdom";
    int wisdom_flag = fftw_import_wisdom_from_filename(filename.c_str());

    // 1 means it's read correctly, 0 indicates a failure
    if (wisdom_flag != 1)
    {
        WARNING("FFTW wisdom file not imported correctly; DFT may take much longer than usual.");
    }

    /*
     * Plan the discrete Fourier transforms:
     *
     *  * Forward are real-to-complex and out-of-place
     *  * Backward are complex-to-complex and-in-place
     */

    // We first set the pointers to arrays that have already been assigned space, as the planning procedure may
    // overwrite data
    mpFftwForwardInX = &(mpArrays->rGetModifiableRightHandSideGrids()[0][0][0]);
    mpFftwForwardInY = &(mpArrays->rGetModifiableRightHandSideGrids()[1][0][0]);

    mpFftwForwardOutX = reinterpret_cast<fftw_complex*>(&(mpArrays->rGetModifiableFourierGrids()[0][0][0]));
    mpFftwForwardOutY = reinterpret_cast<fftw_complex*>(&(mpArrays->rGetModifiableFourierGrids()[1][0][0]));

    mpFftwBackwardInOutX = reinterpret_cast<fftw_complex*>(&(mpArrays->rGetModifiableNewVelocityGrids()[0][0][0]));
    mpFftwBackwardInOutY = reinterpret_cast<fftw_complex*>(&(mpArrays->rGetModifiableNewVelocityGrids()[1][0][0]));

    mFftwForwardPlanX = fftw_plan_dft_r2c_2d(mNumGridPtsX, mNumGridPtsY, mpFftwForwardInX, mpFftwForwardOutX, FFTW_EXHAUSTIVE);
    mFftwForwardPlanY = fftw_plan_dft_r2c_2d(mNumGridPtsX, mNumGridPtsY, mpFftwForwardInY, mpFftwForwardOutY, FFTW_EXHAUSTIVE);

    mFftwInversePlanX = fftw_plan_dft_2d(mNumGridPtsX, mNumGridPtsY, mpFftwBackwardInOutX, mpFftwBackwardInOutX, FFTW_BACKWARD, FFTW_EXHAUSTIVE);
    mFftwInversePlanY = fftw_plan_dft_2d(mNumGridPtsX, mNumGridPtsY, mpFftwBackwardInOutY, mpFftwBackwardInOutY, FFTW_BACKWARD, FFTW_EXHAUSTIVE);


//    // Set up threads for fftw
//    int potential_thread_errors = fftw_init_threads();
//
//    if (potential_thread_errors == 0)
//    {
//        EXCEPTION("fftw thread error");
//    }
//
//    fftw_plan_with_nthreads(2);
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

    multi_array<double, 3>& r_force_grids = mpArrays->rGetModifiableForceGrids();

    for (unsigned x = 0; x < mNumGridPtsX; x++)
    {
        for (unsigned y = 0; y < mNumGridPtsY; y++)
        {
            r_force_grids[0][x][y] = 0.0;
        }
    }
    for (unsigned x = 0; x < mNumGridPtsX; x++)
    {
        for (unsigned y = 0; y < mNumGridPtsY; y++)
        {
            r_force_grids[1][x][y] = 0.0;
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::AddForceContributions()
{
    // Add contributions from each immersed boundary force
    for (typename std::vector<boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > >::iterator iter = mForceCollection.begin();
            iter != mForceCollection.end();
            ++iter)
    {
        (*iter)->AddForceContribution(mNodePairs);
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

    multi_array<double, 3>& force_grids = mpArrays->rGetModifiableForceGrids();

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

                force_grids[0][(first_idx_y + y_idx) % mNumGridPtsY][(first_idx_x + x_idx) % mNumGridPtsX] += force_x;
                force_grids[1][(first_idx_y + y_idx) % mNumGridPtsY][(first_idx_x + x_idx) % mNumGridPtsX] += force_y;
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SolveNavierStokesSpectral()
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned reduced_size = 1 + (mNumGridPtsY/2);


    std::string match;

    // Get non modifiable fluid grids
    const std::vector<std::vector<double> >& VelX = mpMesh->rGetFluidVelocityGridX();
    const std::vector<std::vector<double> >& VelY = mpMesh->rGetFluidVelocityGridY();

    // Get references to all the necessary grids
    multi_array<double, 3>& vel_grids   = mpMesh->rGetModifiable2dVelocityGrids();
    multi_array<double, 3>& force_grids = mpArrays->rGetModifiableForceGrids();
    multi_array<double, 3>& rhs_grids   = mpArrays->rGetModifiableRightHandSideGrids();

    multi_array<std::complex<double>, 3>& fourier_grids = mpArrays->rGetModifiableFourierGrids();
    multi_array<std::complex<double>, 2>& pressure_grid = mpArrays->rGetModifiablePressureGrid();
    multi_array<std::complex<double>, 3>& new_vel_grids  = mpArrays->rGetModifiableNewVelocityGrids();

    for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
    {
        for (unsigned y = 0; y < mNumGridPtsY; y++)
        {
            vel_grids[0][x][y] = VelX[y][x];
            vel_grids[1][x][y] = VelY[y][x];
        }
    }

//    {//DEBUG
//        for (unsigned x = 0; x < mNumGridPtsX; x++)
//        {
//            for (unsigned y = 0; y < mNumGridPtsY; y++)
//            {
//                std::cout << std::setprecision(4) << std::setw(12) << VelX[y][x];
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        for (unsigned x = 0; x < mNumGridPtsX; x++)
//        {
//            for (unsigned y = 0; y < mNumGridPtsY; y++)
//            {
//                std::cout << std::setprecision(4) << std::setw(12) << vel_grids[0][x][y];
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        match = "true";
//        for (unsigned x = 0; x < mNumGridPtsX; x++)
//            for (unsigned y = 0; y < mNumGridPtsY; y++)
//                if (fabs(vel_grids[0][x][y] - VelX[y][x]) > 1e-10)
//                    match = "false";
//        PRINT_VARIABLE(match);
//        std::cout << std::endl;
//    }

    // Create upwind variables
    std::vector<std::vector<double> > Upwind_x;
    std::vector<std::vector<double> > Upwind_y;

    UpwindScheme(VelX, VelY, Upwind_x, Upwind_y);
    Upwind2d(vel_grids, rhs_grids);

    double large_number = 4294967296.0;

//    {//DEBUG
//        for (unsigned x = 0; x < mNumGridPtsX; x++)
//        {
//            for (unsigned y = 0; y < mNumGridPtsY; y++)
//            {
//                std::cout << std::setprecision(4) << std::setw(12) << Upwind_x[y][x];
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        for (unsigned x = 0; x < mNumGridPtsX; x++)
//        {
//            for (unsigned y = 0; y < mNumGridPtsY; y++)
//            {
//                std::cout << std::setprecision(4) << std::setw(12) << rhs_grids[0][x][y];
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        match = "true";
//        for (unsigned x = 0; x < mNumGridPtsX; x++)
//            for (unsigned y = 0; y < mNumGridPtsY; y++)
//                if (fabs(rhs_grids[0][x][y] - Upwind_x[y][x]) > 1e-10)
//                    match = "false";
//        PRINT_VARIABLE(match);
//        std::cout << std::endl;
//    }


    // Create RHS of linear system
    std::vector<std::vector<double> > rhsX; SetupGrid(rhsX);
    std::vector<std::vector<double> > rhsY; SetupGrid(rhsY);

    for(unsigned y = 0 ; y < mNumGridPtsY ; y++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            rhsX[y][x] = large_number * (VelX[y][x] + dt * (mFluidForceGridX[y][x] - Upwind_x[y][x]));
            rhsY[y][x] = large_number * (VelY[y][x] + dt * (mFluidForceGridY[y][x] - Upwind_y[y][x]));
        }
    }

    for (unsigned dim = 0 ; dim < 2 ; dim++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                rhs_grids[dim][x][y] = large_number * (vel_grids[dim][x][y] + dt * (force_grids[dim][x][y] - rhs_grids[dim][x][y]));
            }
        }
    }

    {//DEBUG
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                std::cout << std::setprecision(4) << std::setw(12) << rhsX[y][x];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                std::cout << std::setprecision(4) << std::setw(12) << rhs_grids[0][x][y];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        match = "true";
        for (unsigned x = 0; x < mNumGridPtsX; x++)
            for (unsigned y = 0; y < mNumGridPtsY; y++)
                if (fabs(rhs_grids[0][x][y] - rhsX[y][x]) > 1e-10)
                    match = "false";
        PRINT_VARIABLE(match);
        std::cout << std::endl;
    }

    // Perform fft on rhsX
    std::vector<std::vector<std::complex<double> > > VelX_hat;
    Fft2DForwardRealToComplex(rhsX, VelX_hat);

    // Perform fft on rhsY
    std::vector<std::vector<std::complex<double> > > VelY_hat;
    Fft2DForwardRealToComplex(rhsY, VelY_hat);

    // Perform fft on rhs_grids; results go to fourier_grids
    fftw_execute(mFftwForwardPlanX);
    fftw_execute(mFftwForwardPlanY);

//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << VelX_hat[y][x].real();
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < reduced_size; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << fourier_grids[0][x][y].real();
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;

    std::vector<std::vector<std::complex<double> > > p_hat; SetupGrid(p_hat);

    // Calculate pressure
    for( unsigned y = 0 ; y < mNumGridPtsY ; y++ )
    {
        for( unsigned x = 0 ; x < mNumGridPtsX ; x++ )
        {
            std::complex<double> numerator = -mI * (mSin2X[x] * VelX_hat[y][x] / mGridSpacingX + mSin2Y[y] * VelY_hat[y][x] / mGridSpacingY);

            double denominator = (mSin2X[x] * mSin2X[x] / (mGridSpacingX * mGridSpacingX)) + (mSin2Y[y] * mSin2Y[y] / (mGridSpacingY * mGridSpacingY));
            denominator *= (dt / mReynolds);

            p_hat[y][x] = numerator / denominator;
        }
    }

    /*
     * The result of a DFT of n real datapoints is n/2 + 1 complex values, due to redundancy: element n-1 is conj(2),
     * etc.  Calculation of the pressure grid preserves this property, and so the final pressure grid will take up
     * the first half of the pressure_grids array.
     */

    for (unsigned x = 0 ; x < mNumGridPtsX ; x++)
    {
        for (unsigned y = 0 ; y < reduced_size ; y++)
        {
            std::complex<double> numerator = -mI * (mSin2X[y] * fourier_grids[0][x][y] / mGridSpacingX + mSin2Y[x] * fourier_grids[1][x][y] / mGridSpacingY);

            double denominator = (mSin2X[y] * mSin2X[y] / (mGridSpacingX * mGridSpacingX)) + (mSin2Y[x] * mSin2Y[x] / (mGridSpacingY * mGridSpacingY));
            denominator *= (dt / mReynolds);

            pressure_grid[x][y] = numerator * large_number / denominator;
        }

        for (unsigned y = reduced_size ; y < mNumGridPtsY ; y++)
        {
            unsigned z = mNumGridPtsY - y;
            std::complex<double> numerator = -mI * (mSin2X[y] * conj(fourier_grids[0][x][z]) / mGridSpacingX + mSin2Y[x] * conj(fourier_grids[1][x][z]) / mGridSpacingY);

            double denominator = (mSin2X[y] * mSin2X[y] / (mGridSpacingX * mGridSpacingX)) + (mSin2Y[x] * mSin2Y[x] / (mGridSpacingY * mGridSpacingY));
            denominator *= (dt / mReynolds);

            pressure_grid[x][y] = numerator / denominator;
        }
    }

//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << p_hat[y][x].imag();
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << pressure_grid[x][y].imag();
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;


    // Set some values to zero
    p_hat[0][0] = 0.0;
    p_hat[0][mNumGridPtsX/2] = 0.0;
    p_hat[mNumGridPtsY/2][mNumGridPtsX/2] = 0.0;
    p_hat[mNumGridPtsY/2][0] = 0.0;

    pressure_grid[0][0] = 0.0;
    pressure_grid[mNumGridPtsX/2][0] = 0.0;
    pressure_grid[mNumGridPtsX/2][mNumGridPtsY/2] = 0.0;
    pressure_grid[0][mNumGridPtsY/2] = 0.0;

    // Do final stage of computation before inverse FFT
    std::vector<std::vector<std::complex<double> > > pre_inverse_X; SetupGrid(pre_inverse_X);
    std::vector<std::vector<std::complex<double> > > pre_inverse_Y; SetupGrid(pre_inverse_Y);


    for( unsigned y = 0 ; y < mNumGridPtsY ; y++ )
    {
        for( unsigned x = 0 ; x < mNumGridPtsX ; x++ )
        {
            double op = 1 + (4 * dt / mReynolds) * ( (mSinX[x] * mSinX[x] / (mGridSpacingX * mGridSpacingX)) + (mSinY[y] * mSinY[y] / (mGridSpacingY * mGridSpacingY)));
            pre_inverse_X[y][x] = (VelX_hat[y][x] - (mI * dt / (mReynolds * mGridSpacingX)) * mSin2X[x] * p_hat[y][x]) / op;
            pre_inverse_Y[y][x] = (VelY_hat[y][x] - (mI * dt / (mReynolds * mGridSpacingY)) * mSin2Y[y] * p_hat[y][x]) / op;
        }
    }
    for (unsigned x = 0 ; x < mNumGridPtsX ; x++)
    {
        for (unsigned y = 0 ; y < reduced_size ; y++)
        {
            double op = 1 + (4 * dt / mReynolds) * ( (mSinX[y] * mSinX[y] / (mGridSpacingX * mGridSpacingX)) + (mSinY[x] * mSinY[x] / (mGridSpacingY * mGridSpacingY)));
            new_vel_grids[0][x][y] = (fourier_grids[0][x][y] - (mI * dt / (mReynolds * mGridSpacingX)) * mSin2X[y] * pressure_grid[x][y]) / op;
            new_vel_grids[1][x][y] = (fourier_grids[1][x][y] - (mI * dt / (mReynolds * mGridSpacingY)) * mSin2X[y] * pressure_grid[x][y]) / op;
        }

        for (unsigned y = reduced_size ; y < mNumGridPtsY ; y++)
        {
            unsigned z = mNumGridPtsY - y;
            double op = 1 + (4 * dt / mReynolds) * ( (mSinX[y] * mSinX[y] / (mGridSpacingX * mGridSpacingX)) + (mSinY[x] * mSinY[x] / (mGridSpacingY * mGridSpacingY)));

            new_vel_grids[0][x][y] = (conj(fourier_grids[0][x][z]) - (mI * dt / (mReynolds * mGridSpacingX)) * mSin2X[y] * conj(pressure_grid[x][z])) / op;
            new_vel_grids[1][x][y] = (conj(fourier_grids[1][x][z]) - (mI * dt / (mReynolds * mGridSpacingY)) * mSin2X[y] * conj(pressure_grid[x][z])) / op;
        }
    }

//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << pre_inverse_X[y][x].real();
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << new_vel_grids[0][x][y].real();
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;









    // Perform inverse FFT on x data
    std::vector<std::vector<double> > new_velocity_x;
    Fft2DInverseComplexToReal(pre_inverse_X, new_velocity_x);

    // Perform inverse FFT on y data
    std::vector<std::vector<double> > new_velocity_y;
    Fft2DInverseComplexToReal(pre_inverse_Y, new_velocity_y);

    // Perform inverse fft on new_vel_grids; results are in-place
    fftw_execute(mFftwInversePlanX);
    fftw_execute(mFftwInversePlanY);

//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << new_velocity_y[y][x];
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << new_vel_grids[1][x][y].real() / 64.0;
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;


    // Get modifiable fluid grids
    std::vector<std::vector<double> >& modifiable_vel_x = mpMesh->rGetModifiableFluidVelocityGridX();
    std::vector<std::vector<double> >& modifiable_vel_y = mpMesh->rGetModifiableFluidVelocityGridY();

    for(unsigned y = 0 ; y < mNumGridPtsY ; y++)
    {
        for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
        {
            modifiable_vel_x[y][x] = new_velocity_x[y][x] / large_number;
            modifiable_vel_y[y][x] = new_velocity_y[y][x] / large_number;
        }
    }

    for (unsigned dim = 0 ; dim < 2 ; dim++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < reduced_size; y++)
            {
                vel_grids[dim][x][y] = new_vel_grids[dim][x][y].real() / large_number;;// / (double(mNumGridPtsX) * double(mNumGridPtsY));
            }
        }
    }

//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << modifiable_vel_x[y][x];
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    for (unsigned x = 0; x < mNumGridPtsX; x++)
//    {
//        for (unsigned y = 0; y < mNumGridPtsY; y++)
//        {
//            std::cout << std::setprecision(4) << std::setw(12) << vel_grids[0][x][y];
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
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
void ImmersedBoundarySimulationModifier<DIM>::Upwind2d(const multi_array<double, 3>& input, multi_array<double, 3>& output)
{
    unsigned prev_x = mNumGridPtsX - 1;
    unsigned prev_y = mNumGridPtsY - 1;

    unsigned next_x = 1;
    unsigned next_y = 1;

    for(unsigned x = 0 ; x < mNumGridPtsX ; x++)
    {
        for (unsigned y = 0 ; y < mNumGridPtsY ; y++)
        {
            // Set values for output from conditional on x grid
            if(input[0][x][y] > 0)
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
            if(input[1][x][y] > 0)
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
void ImmersedBoundarySimulationModifier<DIM>::PrintGrid(const std::vector<std::vector<std::complex<double> > > & grid)
{
    for (unsigned y = 0; y < mNumGridPtsY; y++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            std::cout << std::setprecision(10) << grid[y][x].real() << "+" << grid[y][x].imag() << "i ";
        }
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
    mReynolds = reynoldsNumber;
}

template<unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::GetReynoldsNumber()
{
    return mReynolds;
}


// Explicit instantiation
template class ImmersedBoundarySimulationModifier<1>;
template class ImmersedBoundarySimulationModifier<2>;
template class ImmersedBoundarySimulationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundarySimulationModifier)

