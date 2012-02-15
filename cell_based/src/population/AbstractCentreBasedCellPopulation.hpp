/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
     * Constructor for use by archiving only.
     */
    AbstractCentreBasedCellPopulation();

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
    AbstractCentreBasedCellPopulation(std::vector<CellPtr>& rCells,
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
