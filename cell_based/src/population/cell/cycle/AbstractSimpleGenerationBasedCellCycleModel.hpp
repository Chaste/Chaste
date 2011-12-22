/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef ABSTRACTSIMPLEGENERATIONBASEDCYCLEMODEL_HPP_
#define ABSTRACTSIMPLEGENERATIONBASEDCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractSimpleCellCycleModel.hpp"

/**
 * This class contains all the things common to simple generation-based cell cycle
 * models, i.e. models in which the length of cell cycle phases are determined
 * when the cell-cycle model is created, rather than evaluated 'on the fly'
 * by ODEs and suchlike, and in which each cell has a 'generation'.
 *
 * N.B. Whether or not the cell should actually divide may depend on
 * Wnt / Oxygen etc. in subclasses.
 */
class AbstractSimpleGenerationBasedCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        archive & mGeneration;
        archive & mMaxTransitGenerations;
    }

protected:

    /** The generation of this cell (STEM cells have a generation of 0) */
    unsigned mGeneration;

    /** How many generations a transit cell lives for before becoming fully differentiated. */
    unsigned mMaxTransitGenerations;

public:

    /**
     * Default constructor - creates an AbstractSimpleCellCycleModel.
     */
    AbstractSimpleGenerationBasedCellCycleModel();

    /**
     * Destructor.
     */
    virtual ~AbstractSimpleGenerationBasedCellCycleModel();

    /** Overridden ResetForDivision() method. */
    void ResetForDivision();

    /**
     * Set the new cell's G1 duration once it has been created after division.
     * The duration will be based on cell type.
     */
    void InitialiseDaughterCell();

    /**
     * Sets the cell's generation.
     *
     * @param generation the cell's generation
     */
    void SetGeneration(unsigned generation);

    /**
     * Returns the cell's generation.
     */
    unsigned GetGeneration() const;

    /**
     * Set mMaxTransitGenerations.
     *
     * @param maxTransitGenerations the new value of mMaxTransitGenerations
     */
    void SetMaxTransitGenerations(unsigned maxTransitGenerations);

    /**
     * @return mMaxTransitGenerations
     */
    unsigned GetMaxTransitGenerations() const;

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractSimpleGenerationBasedCellCycleModel)

#endif /*ABSTRACTSIMPLEGENERATIONBASEDCYCLEMODEL_HPP_*/
