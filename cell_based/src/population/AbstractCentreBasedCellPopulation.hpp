/*

Copyright (c) 2005-2019, University of Oxford.
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
#include "AbstractCentreBasedDivisionRule.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule; // Circular definition thing.

/**
 * An abstract facade class encapsulating a centre-based cell population, in which
 * each cell corresponds to a Node.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractCentreBasedCellPopulation : public AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * This test uses the private constructor to simplify testing.
     */
    friend class TestCentreBasedDivisionRules;

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
        archive & boost::serialization::base_object<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mMeinekeDivisionSeparation;
        archive & mMarkedSprings;
        archive & mpCentreBasedDivisionRule;
    }

protected:
    /**
     * Initial separation placement of mother/daughter cells at birth.
     * Has units of cell size at equilibrium rest length
     */
    double mMeinekeDivisionSeparation;

    /**
     * Special springs that we want to keep track of for some reason.
     * Currently used to track cells in the process of dividing
     * (which are represented as two cells joined by a shorter spring).
     */
    std::set<std::pair<CellPtr,CellPtr> > mMarkedSprings;

    /** A pointer to a division rule that is used to generate the locations of daughter cells when a cell divides.
     * This is a specialisation for centre-based models.
     */
    boost::shared_ptr<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> > mpCentreBasedDivisionRule;

    /**
     * Constructor that just takes in a mesh.
     *
     * @param rMesh the mesh for the cell population.
     */
    AbstractCentreBasedCellPopulation(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Write the current results to mpVtkMetaFile.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rDirectory  pathname of the output directory, relative to where Chaste output is stored
     */
    virtual void WriteVtkResultsToFile(const std::string& rDirectory)=0;

    /**
     * Call #AcceptCellWriter across the whole population,
     * iterating in an appropriate way for this type of
     * cell population.
     */
    virtual void AcceptCellWritersAcrossPopulation();

public:

    /**
     * Default constructor.
     *
     * @param rMesh a reference to the mesh underlying the cell population
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractCentreBasedCellPopulation(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
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
    c_vector<double, SPACE_DIM> GetLocationOfCellCentre(CellPtr pCell);

    /**
     * Get a pointer to the node corresponding to a given cell.
     *
     * @param pCell the cell
     *
     * @return address of the node
     */
    Node<SPACE_DIM>* GetNodeCorrespondingToCell(CellPtr pCell);

    /**
     * Overridden GetCellDataItemAtPdeNode() method.
     *
     * @param pdeNodeIndex index of a node in a tetrahedral mesh for use
     *         with a PDE modifier
     * @param rVariableName the name of the cell data item to get
     * @param dirichletBoundaryConditionApplies where a Dirichlet boundary condition is used
     *        (optional; defaults to false)
     * @param dirichletBoundaryValue the value of the Dirichlet boundary condition, if used
     *        (optional; defaults to 0.0)
     *
     * @return the value of a CellData item (interpolated if necessary) at a node,
     *         specified by its index in a tetrahedral mesh for use with a PDE modifier.
     * This method can be called by PDE modifier classes.
     */
    virtual double GetCellDataItemAtPdeNode(unsigned pdeNodeIndex,
                                            std::string& rVariableName,
                                            bool dirichletBoundaryConditionApplies=false,
                                            double dirichletBoundaryValue=0.0);

    /**
     * Add a new cell to the cell population.
     *
     * @param pNewCell  the cell to add
     * @param pParentCell pointer to a parent cell (if required)
     *
     * @return address of cell as it appears in the cell list
     */
    CellPtr AddCell(CellPtr pNewCell, CellPtr pParentCell=CellPtr());

    /**
     * @return a an ordered pair of pointers to two given Cells.
     * Used by the spring marking routines.
     * Elements in the returned pair are ordered by cell ID number - the
     * cell in the pair will have a smaller ID.
     *
     * @param pCell1 a Cell
     * @param pCell2 a Cell
     */
    std::pair<CellPtr,CellPtr> CreateCellPair(CellPtr pCell1, CellPtr pCell2);

    /**
     * @param rCellPair a set of pointers to Cells
     *
     * @return whether the spring between two given cells is marked.
     */
    bool IsMarkedSpring(const std::pair<CellPtr,CellPtr>& rCellPair);

    /**
     * Mark the spring between the given cells.
     *
     * @param rCellPair a set of pointers to Cells
     */
    void MarkSpring(std::pair<CellPtr,CellPtr>& rCellPair);

    /**
     * Stop marking the spring between the given cells.
     *
     * @param rCellPair a set of pointers to Cells
     */
    void UnmarkSpring(std::pair<CellPtr,CellPtr>& rCellPair);

    /**
     * Overridden IsCellAssociatedWithADeletedLocation() method.
     *
     * @param pCell the cell
     *
     * @return whether a given cell is associated with a deleted node.
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell);

    /**
     * Overridden GetNeighbouringLocationIndices() method.
     *
     * Given a cell, returns the set of location indices corresponding to neighbouring cells.
     *
     * @param pCell a cell
     * @return the set of neighbouring location indices.
     */
    virtual std::set<unsigned> GetNeighbouringLocationIndices(CellPtr pCell);

    /**
     * Checks whether a given node displacement violates the movement threshold
     * for this population. If so, a stepSizeException is generated that contains
     * a warning/error message and a suggested smaller dt that should avoid the problem.
     *
     * @param nodeIndex Index of the node in question (allows us to check whether this is a ghost or particle)
     * @param rDisplacement Movement vector of the node at this time step
     * @param dt Current time step size
     */
    virtual void CheckForStepSizeException(unsigned nodeIndex, c_vector<double,SPACE_DIM>& rDisplacement, double dt);

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
     * Method to return the connected nodes in  a centre based simulation.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @return Node pairs for force calculation.
     */
    virtual std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >& rGetNodePairs()=0;

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
     * @return The division rule that is currently being used.
     */
    boost::shared_ptr<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> > GetCentreBasedDivisionRule();

    /**
     * Set the division rule for this population.
     *
     * @param pCentreBasedDivisionRule  pointer to the new division rule
     */
    void SetCentreBasedDivisionRule(boost::shared_ptr<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> > pCentreBasedDivisionRule);

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * Overridden GetDefaultTimeStep() method.
     *
     * @return a default value for the time step to use
     * when simulating the cell population.
     *
     * A hard-coded value of 1/120 is returned. However, note that the time
     * step can be reset by calling SetDt() on the simulation object used to
     * simulate the cell population.
     */
    virtual double GetDefaultTimeStep();
};

#endif /*ABSTRACTCENTREBASEDCELLPOPULATION_HPP_*/
