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
#include "Debug.hpp"
#include "Timer.hpp"
#include "FileFinder.hpp"

#include <boost/thread.hpp>

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


    // Set up dimension-dependent variables
    switch (DIM)
    {
        case 2:
            mpArrays = new ImmersedBoundary2dArrays<DIM>(mpMesh, SimulationTime::Instance()->GetTimeStep(), mReynolds);
            mpFftInterface = new ImmersedBoundaryFftInterface<DIM>(mpMesh,
                                                                   &(mpArrays->rGetModifiableRightHandSideGrids()[0][0][0]),
                                                                   &(mpArrays->rGetModifiableFourierGrids()[0][0][0]),
                                                                   &(mpMesh->rGetModifiable2dVelocityGrids()[0][0][0]),
                                                                   2);

            mFftNorm = (double) mNumGridPtsX * (double) mNumGridPtsY;
            break;

        case 3:
            EXCEPTION("Not implemented yet in 3D");
            break;

        default:
            NEVER_REACHED;
    }
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

    multi_array<double, 3>& r_force_grids = mpArrays->rGetModifiableForceGrids();

    for (unsigned dim = 0 ; dim < 2 ; dim++)
    {
        for (unsigned x = 0; x < mNumGridPtsX; x++)
        {
            for (unsigned y = 0; y < mNumGridPtsY; y++)
            {
                r_force_grids[dim][x][y] = 0.0;
            }
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

                force_grids[0][(first_idx_x + x_idx) % mNumGridPtsX][(first_idx_y + y_idx) % mNumGridPtsY] += force_x;
                force_grids[1][(first_idx_x + x_idx) % mNumGridPtsX][(first_idx_y + y_idx) % mNumGridPtsY] += force_y;
            }
        }
    }
}

template<unsigned DIM>
void ImmersedBoundarySimulationModifier<DIM>::SolveNavierStokesSpectral()
{
    double dt = SimulationTime::Instance()->GetTimeStep();
    double large_number = 268435456.0;
    unsigned reduced_size = 1 + (mNumGridPtsY/2);

    // Get references to all the necessary grids
    multi_array<double, 3>& vel_grids   = mpMesh->rGetModifiable2dVelocityGrids();
    multi_array<double, 3>& force_grids = mpArrays->rGetModifiableForceGrids();
    multi_array<double, 3>& rhs_grids   = mpArrays->rGetModifiableRightHandSideGrids();

    const multi_array<double, 2>& op_1  = mpArrays->rGetOperator1();
    const multi_array<double, 2>& op_2  = mpArrays->rGetOperator2();
    const std::vector<double>& sin_2x   = mpArrays->rGetSin2x();
    const std::vector<double>& sin_2y   = mpArrays->rGetSin2y();

    multi_array<std::complex<double>, 3>& fourier_grids = mpArrays->rGetModifiableFourierGrids();
    multi_array<std::complex<double>, 2>& pressure_grid = mpArrays->rGetModifiablePressureGrid();

    // Perform upwind differencing and create RHS of linear system
    Upwind2d(vel_grids, rhs_grids);

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

    // Perform fft on rhs_grids; results go to fourier_grids
    mpFftInterface->FftExecuteForward();

    /*
     * The result of a DFT of n real datapoints is n/2 + 1 complex values, due to redundancy: element n-1 is conj(2),
     * etc.  A similar pattern of redundancy occurs in a 2D transform.  Calculations in the Fourier domain preserve this
     * redundancy, and so all calculations need only be done on reduced-size arrays, saving memory and computation.
     */

    // Calculate the pressure grid
    for (unsigned x = 0 ; x < mNumGridPtsX ; x++)
    {
        for (unsigned y = 0 ; y < reduced_size ; y++)
        {
            pressure_grid[x][y] = -mI * (sin_2x[x] * fourier_grids[0][x][y] / mGridSpacingX + sin_2y[y] * fourier_grids[1][x][y] / mGridSpacingY) / op_1[x][y];
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
    for (unsigned x = 0 ; x < mNumGridPtsX ; x++)
    {
        for (unsigned y = 0 ; y < reduced_size ; y++)
        {
            fourier_grids[0][x][y] = (fourier_grids[0][x][y] - (mI * dt / (mReynolds * mGridSpacingX)) * sin_2x[x] * pressure_grid[x][y]) / (op_2[x][y] * large_number * mFftNorm);
            fourier_grids[1][x][y] = (fourier_grids[1][x][y] - (mI * dt / (mReynolds * mGridSpacingY)) * sin_2y[y] * pressure_grid[x][y]) / (op_2[x][y] * large_number * mFftNorm);
        }
    }

    // Perform inverse fft on fourier_grids; results are in vel_grids
    mpFftInterface->FftExecuteInverse();
}

template<unsigned DIM>
double ImmersedBoundarySimulationModifier<DIM>::Delta1D(double dist, double spacing)
{
    return (0.25 * (1.0 + cos(M_PI * dist / (2 * spacing)))) / spacing;
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

