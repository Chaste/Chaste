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

#ifndef WNTCONCENTRATION_HPP_
#define WNTCONCENTRATION_HPP_

#include "ChasteSerialization.hpp"
#include "SerializableSingleton.hpp"
#include <boost/serialization/base_object.hpp>

#include <iostream>

#include "AbstractCellPopulation.hpp"

/**
 * Possible types of WntConcentration, currently:
 *  NONE - for testing and to remove Wnt dependence
 *  LINEAR - for cylindrical crypt model
 *  RADIAL - for crypt projection model
 */
typedef enum WntConcentrationType_
{
    NONE,
    LINEAR,
    RADIAL,
    EXPONENTIAL
} WntConcentrationType;


/**
 *  Singleton Wnt concentration object.
 */
template<unsigned DIM>
class WntConcentration : public SerializableSingleton<WntConcentration<DIM> >
{
private:

    /** Pointer to the singleton instance of WntConcentration */
    static WntConcentration* mpInstance;

    /**
     * The length of the crypt.
     */
    double mCryptLength;

    /**
     * Whether this WntConcentration object has had its crypt length set.
     */
    bool mLengthSet;

    /**
     * The type of WntConcentration current options are
     * NONE - returns zero everywhere
     * LINEAR - decreases from 1 to zero at height specified by mWntConcentrationParameter
     * RADIAL - decreases from 1 to zero at height specified by mWntConcentrationParameter
     */
    WntConcentrationType mWntType;

    /**
     * The cell population in which the WntConcentration occurs.
     */
    AbstractCellPopulation<DIM>* mpCellPopulation;

    /**
     * Whether this WntConcentration object has had its type set.
     */
    bool mTypeSet;

    /**
     * A value to return for testing purposes.
     */
    double mConstantWntValueForTesting;

    /**
     * Whether to return the testing value
     * (when false WntConcentration works with CellPopulation).
     */
    bool mUseConstantWntValueForTesting;

    /**
     * For LINEAR or RADIAL Wnt type:
     * The proportion of the crypt that has a Wnt gradient.
     * The Wnt concentration goes from one at the base to zero at this height up the crypt.
     *
     * For EXPONENTIAL Wnt type:
     * The parameter lambda in the Wnt concentration
     * Wnt = exp(-height/lambda)
     */
    double mWntConcentrationParameter;

    /**
     * Parameter a, for use in crypt projection simulations, in which the crypt
     * surface is given in cylindrical polar coordinates by z = a*r^b.
     * mCryptProjectionParameterA has no units
     */
    double mCryptProjectionParameterA;

    /**
     * Parameter b, for use in crypt projection simulations, in which the crypt
     * surface is given in cylindrical polar coordinates by z = a*r^b.
     * mCryptProjectionParameterB has no units
     */
    double mCryptProjectionParameterB;

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
        bool is_set_up = IsWntSetUp();
        archive & is_set_up;
        if (is_set_up)
        {
            archive & mCryptLength;
            archive & mLengthSet;
            archive & mWntType;
            archive & mpCellPopulation;
            archive & mTypeSet;
            archive & mConstantWntValueForTesting;
            archive & mUseConstantWntValueForTesting;
            archive & mWntConcentrationParameter;
            archive & mCryptProjectionParameterA;
            archive & mCryptProjectionParameterB;
        }
    }

protected:

    /**
     * Protected constuctor. Not to be called, use Instance() instead.
     */
    WntConcentration();

public:

    /**
     * Return a pointer to the WntConcentration object.
     * The first time this is called, the object is created.
     *
     * @return  A pointer to the singleton WntConcentration object.
     */
    static WntConcentration* Instance();

    /**
     * Destructor - frees up the singleton instance.
     */
    virtual ~WntConcentration();

    /**
     * Destroy the current WntConcentration instance.
     *  Should be called at the end of a simulation.
     */
    static void Destroy();

    /**
     * Get the Wnt level at a given height in the crypt.
     *
     * @param height the height of the cell at which we want the Wnt concentration
     * @return the Wnt concentration at this height in the crypt (dimensionless)
     */
    double GetWntLevel(double height);

    /**
     * Get the Wnt level at a given cell in the crypt. The crypt
     * must be set for this.
     *
     * @param pCell the cell at which we want the Wnt concentration
     * @return the Wnt concentration at this cell
     */
    double GetWntLevel(CellPtr pCell);

    /**
     * Get the Wnt gradient at a given location in the crypt.
     *
     * @param rLocation  the location at which we want the Wnt gradient
     */
    c_vector<double, DIM> GetWntGradient(c_vector<double, DIM>& rLocation);

    /**
     * Get the Wnt gradient at a given cell in the crypt.
     *
     * @param pCell the cell at which we want the Wnt gradient
     */
    c_vector<double, DIM> GetWntGradient(CellPtr pCell);

    /**
     * Set the crypt. Must be called before GetWntLevel().
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return reference to the CellPopulation.
     */
    AbstractCellPopulation<DIM>& rGetCellPopulation();

    /**
     * @return mCryptLength
     */
    double GetCryptLength();

    /**
     * Set mCryptLength. Must be called before GetWntLevel().
     *
     * @param cryptLength  the new value of mCryptLength
     */
    void SetCryptLength(double cryptLength);

    /**
     * Get the type of Wnt concentration.
     */
    WntConcentrationType GetType();

    /**
     * Set the type of Wnt concentration. Must be called before GetWntLevel().
     *
     * @param type the type of Wnt concentration
     */
    void SetType(WntConcentrationType type);

    /**
     * Force the Wnt concentration to return a given value for all cells.
     * Only for testing.
     *
     * @param value the constant value to set the Wnt concentration to be
     */
    void SetConstantWntValueForTesting(double value);

    /**
     * Whether a Wnt concentration has been set up.
     *
     * For archiving, and to let a CellBasedSimulation
     * find out whether whether a WntConcentration has
     * been set up or not, i.e. whether stem cells should
     * be motile.
     *
     * @return whether the Wnt concentration is set up
     */
    bool IsWntSetUp();

    /**
     * @return mWntConcentrationParameter
     */
    double GetWntConcentrationParameter();

    /**
     * Set mWntConcentrationParameter.
     *
     * @param wntConcentrationParameter the new value of mWntConcentrationParameter
     */
    void SetWntConcentrationParameter(double wntConcentrationParameter);

    /**
     * @return mCryptProjectionParameterA
     */
    double GetCryptProjectionParameterA();

    /**
     * @return mCryptProjectionParameterB
     */
    double GetCryptProjectionParameterB();

    /**
     * Set mCryptProjectionParameterA.
     *
     * @param cryptProjectionParameterA  the new value of mCryptProjectionParameterA
     */
    void SetCryptProjectionParameterA(double cryptProjectionParameterA);

    /**
     * Set mCryptProjectionParameterB.
     *
     * @param cryptProjectionParameterB  the new value of mCryptProjectionParameterB
     */
    void SetCryptProjectionParameterB(double cryptProjectionParameterB);
};

#endif /*WNTCONCENTRATION_HPP_*/
