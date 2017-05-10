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

#ifndef IMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_
#define IMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_

// Chaste includes
#include "AbstractCellBasedSimulationModifier.hpp"
#include "ObsoleteBoxCollection.hpp"
#include "ChasteSerialization.hpp"

// Immersed boundary includes
#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundary2dArrays.hpp"
#include "ImmersedBoundaryFftInterface.hpp"

// Other includes
#include <complex>
#include <boost/serialization/base_object.hpp>
#include <boost/multi_array.hpp>

/**
 * A modifier class which at each simulation time step implements the immersed
 * boundary algorithm similar to Rejniak et al (2004). A computational model of
 * the mechanics of growth of the villous trophoblast bilayer. Bull. Math. Biol.
 * 66:199–232. doi:10.1016/j.bulm.2003.06.001.
 */
template<unsigned DIM>
class ImmersedBoundarySimulationModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

    /** To allow tests to directly access solver methods */
    friend class TestImmersedBoundaryPdeSolveMethods;
    friend class TestImmersedBoundarySimulationModifier;

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

    /** Pointer to the immersed boundary mesh */
    ImmersedBoundaryMesh<DIM,DIM>* mpMesh;

    /** Pointer to the immersed boundary cell population */
    ImmersedBoundaryCellPopulation<DIM>* mpCellPopulation;

    /** How often we calculate which cells are neighbours */
    unsigned mNodeNeighbourUpdateFrequency;

    /**
     * Number of grid points in the x direction.
     *
     * \todo Duplicate of ImmersedBoundaryMesh::mNumGridPtsX - remove?
     */
    unsigned mNumGridPtsX;

    /**
     * Number of grid points in the y direction.
     *
     * \todo Duplicate of ImmersedBoundaryMesh::mNumGridPtY - remove?
     */
    unsigned mNumGridPtsY;

    /** Number of grid points in the x direction */
    double mGridSpacingX;

    /** Number of grid points in the y direction */
    double mGridSpacingY;

    /** Normalising constant needed for FFT */
    double mFftNorm;

    /** Whether to apply multiplicative normal noise to the calculated force */
    bool mAdditiveNormalNoise;

    /** The mean of the Normal distribution from which random noise variations are drawn */
    double mNormalNoiseMean;

    /** The standard deviation of the Normal distribution from which random noise variations are drawn */
    double mNormalNoiseStdDev;

    /** A box collection to efficiently keep track of node neighbours */
    ObsoleteBoxCollection<DIM>* mpBoxCollection;

    /** A vector of pairs of pointers to nodes, representing all possible node-node interactions */
    std::vector<std::pair<Node<DIM>*, Node<DIM>*> > mNodePairs;

    /** A map between node indices and a set of their possible neighbours, used calculating cell-cell interactions */
    std::map<unsigned, std::set<unsigned> > mNodeNeighbours;

    /**
     * The fluid Reynolds number.
     *
     * Initialised to 1e-4 in the constructor.
     */
    double mReynoldsNumber;

    /** Imaginary unit. */
    std::complex<double> mI;

    /** A list of force laws to determine the force applied to each node */
    std::vector<boost::shared_ptr<AbstractImmersedBoundaryForce<DIM> > > mForceCollection;

    /** Pointer to structure storing all necessary arrays */
    ImmersedBoundary2dArrays<DIM>* mpArrays;

    ///\todo Document class member
    ImmersedBoundaryFftInterface<DIM>* mpFftInterface;

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
     * Ensures fluid force grids and source grid are reset
     */
    void ClearForcesAndSources();

    /**
     * Loops over each immersed boundary force and invokes AddImmersedBoundaryForceContribution()
     */
    void AddImmersedBoundaryForceContributions();

    /**
     * Helper method for UpdateFluidVelocityGrids()
     * Propagates elastic forces to fluid grid
     */
    void PropagateForcesToFluidGrid();

    /**
     * Helper method for UpdateFluidVelocityGrids()
     * Propagates fluid sources to grid
     */
    void PropagateFluidSourcesToGrid();

    /**
     * Helper method for UpdateFluidVelocityGrids()
     * Updates fluid velocity grids by solving Navier-Stokes
     */
    void SolveNavierStokesSpectral();

    /**
     * Helper method for PropagateForcesToFluidGrid()
     * Calculates the discrete delta approximation based on distance and grid spacing
     *
     * @param dist absolute 1-D distance between boundary-node and fluid-node
     * @param spacing the grid spacing
     */
    double Delta1D(double dist, double spacing);

    /**
     * Calculates upwind difference of fluid velocity grids
     *
     * @param input const reference to input grids
     * @param output reference to output grids
     */
    void Upwind2d(const multi_array<double, 3>& input, multi_array<double, 3>& output);

    /**
     * Calculates the vector of central differences of the fluid source grid
     *
     * @param rhs const reference to rhs grids which contain the fluid source strengths
     * @param gradients reference to grids storing the graidents
     */
    void CalculateSourceGradients(const multi_array<double, 3>& rhs, multi_array<double, 3>& gradients);

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
     * Set #mReynoldsNumber.
     *
     * @param reynoldsNumber the new Reynolds number
     */
    void SetReynoldsNumber(double reynoldsNumber);

    /**
     * @return #mReynoldsNumber
     */
    double GetReynoldsNumber();

    /** @return mAdditiveNormalNoise */
    bool GetAdditiveNormalNoise() const;

    /** @param additiveNormalNoise whether to include multiplicative normal noise */
    void SetAdditiveNormalNoise(bool additiveNormalNoise);

    /** @return mNormalNoiseMean */
    double GetNormalNoiseMean() const;

    /** @param normalNoiseMean the new value of mNormalNoiseMean */
    void SetNormalNoiseMean(double normalNoiseMean);

    /** @return mNormalNoiseStdDev */
    double GetNormalNoiseStdDev() const;

    /** @param normalNoiseStdDev the new value of mNormalNoiseStdDev */
    void SetNormalNoiseStdDev(double normalNoiseStdDev);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundarySimulationModifier)

#endif /*IMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_*/
