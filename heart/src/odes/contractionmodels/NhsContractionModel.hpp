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
#ifndef NHSCONTRACTIONMODEL_HPP_
#define NHSCONTRACTIONMODEL_HPP_

#include "AbstractOdeBasedContractionModel.hpp"
#include "Exception.hpp"


/**
 *  NHS (Niederer, Hunter, Smith) model of active tension in cardiac cells.
 *
 *  A system of ODEs which determines the active potential, given the intracellular
 *  calcium concentration, the stretch (lambda) of the cell, and the stretch rate
 *  (dlambda_dt) of the cell.
 *
 *  The state variables are, in order: Calcium_troponin, z, Q1, Q2, Q3
 *
 *  Reference: S.A. Niederer, N.P. Smith, P.J. Hunter, "New developments in a strongly
 *  coupled cardiac electro-mechanical model" Europace 7, S118-S127
 *
 *  The active tension is returned in KPa.
 */
class NhsContractionModel  : public AbstractOdeBasedContractionModel
{
friend class TestContractionModels;

protected :
    /** The stretch. To be specified by the caller */
    double mLambda;
    /** The stretch rate. To be specified by the caller */
    double mDLambdaDt;
    /** The intracellular calcium concentration. To be specified by the caller */
    double mCalciumI;


    /** A parameter only dependent on constants and lambda, so updated whenever lambda is updated */
    double mCalciumTrop50;

    /** A constant determined from the other constrants. Set up in the constructor */
    double mK1;
    /** A constant determined from the other constrants. Set up in the constructor */
    double mK2;

    // Parameters

    /** See reference. (mMols)^-1 (ms)^-1 */
    static const double mKon;

    /** See reference. (ms)^-1 */
    static const double mKrefoff;

    /** See reference. Dimensionless */
    static const double mGamma;

    /** See reference. mMols */
    static const double mCalciumTroponinMax;

    /** See reference. (ms)^-1 */
    static const double mAlphaR1;

    /** See reference. (ms)^-1 */
    static const double mAlphaR2;

    /** See reference. Dimensionless */
    static const double mKZ;

    /** See reference. Dimensionless */
    static const unsigned mNr;

    /** See reference. Dimensionless */
    static const double mBeta1;

    /** See reference. (ms)^-1 */
    static const double mAlpha0;

    /** See reference. Dimensionless */
    static const unsigned mN;

    /** See reference. Dimensionless */
    static const double mZp;

    /** See reference. mMols */
    static const double mCalcium50ref;

    /** See reference. kPa */
    static const double mTref;

    /** See reference. Dimensionless */
    static const double mBeta0;

    /** See reference. Dimensionless */
    static const double mA;

    /** See reference. Dimensionless */
    static const double mA1;

    /** See reference. Dimensionless */
    static const double mA2;

    /** See reference. Dimensionless */
    static const double mA3;

    /** See reference. (ms)^-1 */
    static const double mAlpha1;

    /** See reference. (ms)^-1 */
    static const double mAlpha2;

    /** See reference. (ms)^-1 */
    static const double mAlpha3;

    /**
     *  Compute the calcium_trop50 concentration. This is a function of constants and
     *  lambda, so only needs to be called in the constructor or when lambda is set
     */
    void CalculateCalciumTrop50();

    /**
     *  Calculate T0. This is a function of constants, lambda and z
     *
     *  @param z
     */
    double CalculateT0(double z);

public :
    /**
     *  Constructor. Initialises all state variables to zero, lambda to 1, dlambda_dt
     *  to 0 and intracellular calcium concentration to 0
     */
    NhsContractionModel();

    /**
     *  Set the current stretch and the stretch rate of the cell/fibre
     *
     * @param lambda  current stretch
     * @param dlambdaDt  current stretch rate
     */
    void SetStretchAndStretchRate(double lambda, double dlambdaDt);

    /**
     *  Set the current intracellular calcium concentration
     *
     *  @param rInputParameters input parameters (calcium, voltage, time, of which only calcium is used)
     */
    void SetInputParameters(ContractionModelInputParameters& rInputParameters);

    /**
     *  Directly set the intracellular calcium concentration.
     *  @param calciumConcentration calcium concentration.
     */
    void SetIntracellularCalciumConcentration(double calciumConcentration);

    /**
     *  Get the current Calcium Troponin (one of the state variables) value. This
     *  may be needed if the cell model has Calcium troponin and might need overwriting
     */
    double GetCalciumTroponinValue();

    /**
     * Evaluate the derivatives of the state variables
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    /**
     *  Get the active tension, which is a function of the constants and current state variables. KILOPASCALS
     */
    double GetActiveTension();

    /**
     *  GetNextActiveTension() normally returns the active tension corresponding to the state variables
     *  that have been computed in RunDoNotUpdate. However, this only applies to when an implicit cardiac
     *  mechanics solver is used, in which case the NhsModelWithBackwardSolver should be used.
     */
    double GetNextActiveTension()
    {
        EXCEPTION("If using this in an 'explicit manner' call UpdateStateVariables() and then GetActiveTension(), otherwise use NhsModelWithBackwardSolver");
    }

   /**
    *  This model is stretch-dependent
    */
    bool IsStretchDependent()
    {
        return true;
    }

   /**
    *  This model is stretch-rate-dependent
    */
    bool IsStretchRateDependent()
    {
        return true;
    }
};
#endif /*NHSCONTRACTIONMODEL_HPP_*/
