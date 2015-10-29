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

#ifndef IMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_
#define IMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_


// Chaste includes
#include "AbstractCellBasedSimulationModifier.hpp"
#include "BoxCollection.hpp"
#include "ChasteSerialization.hpp"

// Immersed boundary includes
#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundary2dArrays.hpp"

// Other includes
#include <complex>
#include <fftw3.h>
#include <boost/serialization/base_object.hpp>
#include <boost/multi_array.hpp>

//typedef boost::multi_array<std::complex<double>, 2> boost_complex_array;
//typedef boost::multi_array<double, 2> boost_real_array;

/**
 * A modifier class which at each simulation time step implements the immersed boundary algorithm similar to
 * Rejniak, K. A., Kliman, H. J., & Fauci, L. J. (2004). A computational model of the mechanics of growth of the villous trophoblast bilayer.
 * Bulletin of Mathematical Biology, 66(2), 199–232. doi:10.1016/j.bulm.2003.06.001
 */
template<unsigned DIM>
class ImmersedBoundarySimulationModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** To allow tests to directly access solver methods */
    friend class TestImmersedBoundaryPdeSolveMethods;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

private:
    /** Pointer to the immersed boundary mesh */
    ImmersedBoundaryMesh<DIM,DIM>* mpMesh;

    /** Pointer to the immersed boundary cell population */
    ImmersedBoundaryCellPopulation<DIM>* mpCellPopulation;

    /** How often we calculate which cells are neighbours */
    unsigned mNodeNeighbourUpdateFrequency;

    /** Number of grid points in the x direction */
    unsigned mNumGridPtsX;

    /** Number of grid points in the y direction */
    unsigned mNumGridPtsY;

    /** Number of grid points in the x direction */
    double mGridSpacingX;

    /** Number of grid points in the y direction */
    double mGridSpacingY;

    /** Normalising constant needed for FFT */
    double mFftNorm;

    /** Vector of sin(pi x / Nx) needed for FFT, but constant after grid size is known */
    std::vector<double> mSinX;

    /** Vector of sin(2pi x / Nx) needed for FFT, but constant after grid size is known */
    std::vector<double> mSin2X;

    /** Vector of sin(pi y / Ny) needed for FFT, but constant after grid size is known */
    std::vector<double> mSinY;

    /** Vector of sin(2pi y / Ny) needed for FFT, but constant after grid size is known */
    std::vector<double> mSin2Y;

    /** A vector of pairs of pointers to nodes, representing all possible node-node interactions */
    std::vector<std::pair<Node<DIM>*, Node<DIM>*> > mNodePairs;

    /** A box collection to efficiently keep track of node neighbours */
    BoxCollection<DIM>* mpBoxCollection;

    /** A map between node indices and a set of their possible neighbours, used calculating cell-cell interactions */
    std::map<unsigned, std::set<unsigned> > mNodeNeighbours;

    /** The fluid Reynolds number */
    double mReynolds;

    /** Imaginary unit */
    std::complex<double> mI;

    /** A list of force laws to determine the force applied to each node */
    std::vector<boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > > mForceCollection;

    /** The fftw plan for forward transforms in X */
    fftw_plan mFftwForwardPlanX;

    /** The fftw plan for forward transforms in Y */
    fftw_plan mFftwForwardPlanY;

    /** The fftw plan for inverse transforms in X */
    fftw_plan mFftwInversePlanX;

    /** The fftw plan for inverse transforms in Y */
    fftw_plan mFftwInversePlanY;

    /** Pointer to structure storing all necessary arrays */
    ImmersedBoundary2dArrays* mpArrays;

    /**
     * Helper method to calculate elastic forces, propagate these to the fluid grid
     * and solve Navier-Stokes to update the fluid velocity grids
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateFluidVelocityGrids(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method for SetupSolve()
     * Sets up all variables which need not change throughout the simulation
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetupConstantMemberVariables(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method for UpdateFluidVelocityGrids()
     * Ensures force applied to each node is reset to zero
     * Ensures fluid force grids are reset
     */
    void ClearForces();

    /**
     * Loops over each immersed boundary force and invokes AddForceContribution()
     */
    void AddForceContributions();

    /**
     * Helper method for UpdateFluidVelocityGrids()
     * Propagates elastic forces to fluid grid
     */
    void PropagateForcesToFluidGrid();

    /**
     * Helper method for UpdateFluidVelocityGrids()
     * Updates fluid velocity grids by solving Navier-Stokes
     */
    void SolveNavierStokesSpectral();

    /**
     * Helper method for SolveNavierStokesSpectral()
     * Performs forward fourier transforms.  In separate method only for profiling
     */
    void FftwForward();

    /**
     * Helper method for SolveNavierStokesSpectral()
     * Performs inverse fourier transforms.  In separate method only for profiling
     */
    void FftwInverse();

    /**
     * Helper method for PropagateForcesToFluidGrid()
     * Calculates the discrete delta approximation based on distance and grid spacing
     *
     * @param absolute 1-D distance between boundary-node and fluid-node
     * @param the grid spacing
     */
    double Delta1D(double dist, double spacing);

    /**
     * Calculates upwind difference
     *
     * @param const reference to input grids
     * @param reference to output grids
     */
    void Upwind2d(const multi_array<double, 3>& input, multi_array<double, 3>& output);

    /**
     * Helper method to set member variables for testing purposes
     *
     * @param numGridPtsY the number of grid points in the y direction
     * @param numGridPtsX the number of grid points in the x direction
     */
    void SetMemberVariablesForTesting(unsigned numGridPtsY, unsigned numGridPtsX);


public:

    /**
     * Default constructor.
     */
    ImmersedBoundarySimulationModifier();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundarySimulationModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    /**
     * @param newFrequency the new number of time steps after which node neighbours are re-calculated
     */
    void SetNodeNeighbourUpdateFrequency(unsigned newFrequency);

    /**
     * @return the current number of time steps after which node neighbours are re-calculated
     */
    unsigned GetNodeNeighbourUpdateFrequency();

    /**
     * Add an immersed boundary force to be used in this modifier.
     *
     * @param pForce pointer to a force law
     */
    void AddImmersedBoundaryForce(boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > pForce);

    /**
     * @param reynoldsNumber the new Reynolds number
     */
    void SetReynoldsNumber(double reynoldsNumber);

    /**
     * @return the current Reynolds number
     */
    double GetReynoldsNumber();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundarySimulationModifier)

#endif /*IMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_*/
