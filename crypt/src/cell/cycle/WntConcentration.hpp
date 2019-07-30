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
     * @return the Wnt gradient at a given location in the crypt.
     *
     * @param rLocation  the location at which we want the Wnt gradient
     */
    c_vector<double, DIM> GetWntGradient(c_vector<double, DIM>& rLocation);

    /**
     * @return the Wnt gradient at a given cell in the crypt.
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
     * @return the type of Wnt concentration.
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
