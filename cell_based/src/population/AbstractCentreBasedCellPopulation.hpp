/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef ABSTRACTCENTREBASEDCELLPOPULATION_HPP_
#define ABSTRACTCENTREBASEDCELLPOPULATION_HPP_

#include "AbstractOffLatticeCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"

/**
 * An abstract facade class encapsulating a centre-based cell population, in which
 * each cell corresponds to a Node.
 */
template<unsigned DIM>
class AbstractCentreBasedCellPopulation : public AbstractOffLatticeCellPopulation<DIM>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOffLatticeCellPopulation<DIM> >(*this);
        archive & mMeinekeDivisionSeparation;
    }

protected:

    /**
     * Initial separation placement of mother/daughter cells at birth.
     * Has units of cell size at equilibrium rest length
     */
    double mMeinekeDivisionSeparation;

    /**
     * Constructor that just takes in a mesh.
     */
    AbstractCentreBasedCellPopulation(AbstractMesh<DIM, DIM>& rMesh);

    /**
     * Write the current results to mpVtkMetaFile.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void WriteVtkResultsToFile()=0;

public:

    /**
     * Default constructor.
     *
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractCentreBasedCellPopulation( AbstractMesh<DIM, DIM>& rMesh,
										std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Overridden GetLocationOfCellCentre() method.
     * Find where a given cell is in space.
     *
     * @param pCell the cell
     *
     * @return the location of the node corresponding to this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell);

    /**
     * Get a pointer to the node corresponding to a given cell.
     *
     * @param pCell the cell
     *
     * @return address of the node
     */
    Node<DIM>* GetNodeCorrespondingToCell(CellPtr pCell);

    /**
     * Add a new cell to the cell population.
     *
     * @param pNewCell  the cell to add
     * @param rCellDivisionVector  the position in space at which to put it
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list
     */
    CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell=CellPtr());

    /**
     * Overridden IsCellAssociatedWithADeletedLocation() method.
     *
     * @param pCell the cell
     *
     * @return whether a given cell is associated with a deleted node.
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell);

    /**
     * Overridden UpdateNodeLocations() method.
     *
     * @param rNodeForces a vector containing the force on each node in the cell population
     * @param dt the time step
     */
    virtual void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt);

    /**
     * Overridden GetDampingConstant() method.
     *
     * Get the damping constant for the cell associated with this node,
     * i.e. d in drdt = F/d.
     *
     * If the cell is wild-type, then the normal damping constant is returned. If
     * the cell has a mutation, then the mutant damping constant is returned.
     *
     * Note that by default, the normal and mutant damping constants are the same.
     * To alter the damping constant for mutant cells, call the method
     * SetDampingConstantMutant() on the CellPopulation.
     *
     * @param nodeIndex the global index of this node
     * @return the damping constant at the Cell associated with this node
     */
    virtual double GetDampingConstant(unsigned nodeIndex);

    /**
     * Overridden GenerateCellResults() method.
     *  Generate results for a given cell in the current cell population state to output files.
     *
     * @param locationIndex location index of the cell
     * @param rCellProliferativeTypeCounter cell type counter
     * @param rCellCyclePhaseCounter cell cycle phase counter
     */
    void GenerateCellResults(unsigned locationIndex,
                             std::vector<unsigned>& rCellProliferativeTypeCounter,
                             std::vector<unsigned>& rCellCyclePhaseCounter);

    /**
     * Overridden GenerateCellResultsAndWriteToFiles() method.
     */
    virtual void GenerateCellResultsAndWriteToFiles();

    /**
     * Overridden WriteTimeAndNodeResultsToFiles() method.
     */
    virtual void WriteTimeAndNodeResultsToFiles();

    /**
     * Find if a given node is a ghost node. The method always returns false
     * but is overridden in MeshBasedCellPopulationWithGhostNodes.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is a ghost node
     */
    virtual bool IsGhostNode(unsigned index);

    /**
     * Find if a given node is a particle. The method always returns false
     * but is overridden in NodeBasedCellPopulationWithParticles.
     *
     * @param index the global index of a specified node
     *
     * @return whether the node is a particle
     */
    virtual bool IsParticle(unsigned index);

    /**
     * @return mMeinekeDivisionSeparation
     */
    double GetMeinekeDivisionSeparation();

    /**
     * Set mMeinekeDivisionSeparation.
     *
     * @param divisionSeparation the new value of mMeinekeDivisionSeparation
     */
    void SetMeinekeDivisionSeparation(double divisionSeparation);

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#endif /*ABSTRACTCENTREBASEDCELLPOPULATION_HPP_*/
