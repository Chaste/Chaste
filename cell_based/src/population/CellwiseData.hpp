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

#ifndef CELLWISEDATA_HPP_
#define CELLWISEDATA_HPP_

#include "ChasteSerialization.hpp"
#include "SerializableSingleton.hpp"
#include <boost/serialization/vector.hpp>

#include "MeshBasedCellPopulation.hpp"

/**
 * A singleton object for storing data that certain cell-cycle models
 * need to know about, e.g. nutrient concentrations computed via some PDE
 * for use in nutrient-based cell-cycle models.
 */
template<unsigned DIM>
class CellwiseData : public SerializableSingleton<CellwiseData<DIM> >
{
    friend class TestCellwiseData;
private:

    /** The single instance of the singleton object */
    static CellwiseData* mpInstance;

    /** A pointer to a CellPopulation so a cell's node can be found */
    AbstractCellPopulation<DIM>* mpCellPopulation;

    /** Allocated memory for mData object */
    bool mAllocatedMemory;

    /** Number of variables per node to be stored */
    unsigned mNumberOfVariables;

    /** Store of the data */
    std::vector<double> mData;

    /** Helper member storing constant data. Used in tests. */
    std::vector<double> mConstantDataForTesting;

    /** Helper member storing whether mConstantDataForTesting is used. */
    bool mUseConstantDataForTesting;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        bool is_set_up = IsSetUp();
        archive & is_set_up;
        if (is_set_up)
        {
            archive & mpCellPopulation;
            archive & mAllocatedMemory;
            archive & mNumberOfVariables;
            archive & mData;
            archive & mConstantDataForTesting;
            archive & mUseConstantDataForTesting;
        }
    }

protected:

    /**
     * Protected constuctor. Not to be called, use Instance() instead.
     */
    CellwiseData();

public:

    /**
     * Get an instance of the object
     */
    static CellwiseData* Instance();

    /**
     * Destructor.
     */
    virtual ~CellwiseData();

    /**
     * Destroy the current instance. Should be called at the end of a
     * simulation.
     */
    static void Destroy();

    /**
     * Get the value of CellwiseData for a given cell and variable number.
     *
     * @param pCell the cell
     * @param variableNumber the index of CellwiseData whose value is required (defaults to zero)
     *
     * @return the value of CellwiseData.
     */
    double GetValue(CellPtr pCell, unsigned variableNumber=0);

    /**
     * Set the value for a given location index and variable number.
     *
     * @param value the value to set
     * @param locationIndex the location index
     * @param variableNumber the index of CellwiseData whose value is set (defaults to zero)
     */
    void SetValue(double value, unsigned locationIndex, unsigned variableNumber=0);

    /**
     * Set the CellPopulation. Must be called before GetValue().
     *
     * @param pCellPopulation pointer to the CellPopulation
     */
    void SetCellPopulation(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * @return reference to the CellPopulation.
     */
    AbstractCellPopulation<DIM>& rGetCellPopulation();

    /**
     * Set the number of cells and number of variables to be stored per cell. The constructor
     * assumes 1 variable so this method only really needs to be called if numVars > 1.
     *
     * @param numCells number of cells in the cell population
     * @param numVars number of variables
     * @param pCellPopulation the cell population - this is a temporary variable and will be removed on re-factor see ticket #1515
     */
    void SetNumCellsAndVars(unsigned numCells, unsigned numVars, AbstractCellPopulation<DIM>* pCellPopulation = NULL);

    /**
     * Force the data to return given values for all cells (only for testing).
     *
     * @param rValues vector of CellwiseData values
     */
    void SetConstantDataForTesting(std::vector<double>& rValues);

    /**
     * Is the instance in existence and fully set up
     */
    bool IsSetUp();

    /**
     * Reallocate size of mData according to the number of cells in the CellPopulation
     * member variable. Needed because of possible cell division and cell death over
     * the course of each time step.
     */
    void ReallocateMemory();

    /**
     * @return mNumberOfVariables
     */
    unsigned GetNumVariables();
};

#endif /*CELLWISEDATA_HPP_*/
