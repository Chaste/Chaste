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

#ifndef ABSTRACTOFFLATTICECELLPOPULATION_HPP_
#define ABSTRACTOFFLATTICECELLPOPULATION_HPP_

#include "AbstractCellPopulation.hpp"

/**
 * An abstract facade class encapsulating an off-lattice (centre- or
 * vertex-based) cell population.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractOffLatticeCellPopulation : public AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>
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
        archive & boost::serialization::base_object<AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mDampingConstantNormal;
        archive & mDampingConstantMutant;
        archive & mAbsoluteMovementThreshold;
    }

protected:
    /**
     * Damping constant for normal cells has units of kg s^-1
     * Represented by the parameter eta in the model by Meineke et al (2001) in
     * their off-lattice model of the intestinal crypt
     * (doi:10.1046/j.0960-7722.2001.00216.x).
     */
    double mDampingConstantNormal;

    /**
     * Damping constant for mutant cells has units of kg s^-1.
     */
    double mDampingConstantMutant;

    /**
     * The absolute distance which a cell is permitted to move in one time-step.
     * Movement beyond this threshold will trigger UpdateNodeLocations to throw an exception
     *
     */
    double mAbsoluteMovementThreshold;

    /**
     * Constructor that just takes in a mesh.
     *
     * @param rMesh the mesh for the cell population.
     */
    AbstractOffLatticeCellPopulation(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

public:

    /**
     * Default constructor.
     *
     * @param rMesh a refernce to the mesh underlying the cell population
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    AbstractOffLatticeCellPopulation( AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                    std::vector<CellPtr>& rCells,
                                    const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Add a new node to the cell population.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population.
     */
    virtual unsigned AddNode(Node<SPACE_DIM>* pNewNode)=0;

    /**
     * Move the node with a given index to a new point in space.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    virtual void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM>& rNewLocation)=0;

    /**
     * This is a deprecated method and only used in NodesOnlyCellPopulationWithBuskeUpdate.
     *
     * Update the location of each node in the cell population given
     * a vector of forces on nodes and a time step over which to
     * integrate the equations of motion.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param dt time step
     */
    virtual void UpdateNodeLocations(double dt);

    /**
     * Checks whether a given node displacement violates the movement threshold
     * for this population. If so, a stepSizeException is generated that contains
     * a warning/error message and a suggested smaller dt that should avoid the problem.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex Index of the node in question (allows us to check whether this is a ghost or particle)
     * @param rDisplacement Movement vector of the node at this time step
     * @param dt Current time step size
     */
    virtual void CheckForStepSizeException(unsigned nodeIndex, c_vector<double,SPACE_DIM>& rDisplacement, double dt)=0;

    /**
     * Get the damping constant for this node - ie d in drdt = F/d.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param nodeIndex the global index of this node
     * @return the damping constant at the node.
     */
    virtual double GetDampingConstant(unsigned nodeIndex)=0;

    /**
      * Set mDampingConstantNormal.
      *
      * @param dampingConstantNormal  the new value of mDampingConstantNormal
     */
    void SetDampingConstantNormal(double dampingConstantNormal);

    /**
     * Set mDampingConstantMutant.
     *
     * @param dampingConstantMutant  the new value of mDampingConstantMutant
     */
    void SetDampingConstantMutant(double dampingConstantMutant);

    /**
     * Set mAbsoluteMovementThreshold.
     *
     * @param absoluteMovementThreshold the new value of mAbsoluteMovementThreshold
     */
    void SetAbsoluteMovementThreshold(double absoluteMovementThreshold);

    /**
     * @return mAbsoluteMovementThreshold
     */
    double GetAbsoluteMovementThreshold();

    /**
     * @return mDampingConstantNormal
     */
    double GetDampingConstantNormal();

    /**
     * @return mDampingConstantMutant
     */
    double GetDampingConstantMutant();

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#endif /*ABSTRACTOFFLATTICECELLPOPULATION_HPP_*/
