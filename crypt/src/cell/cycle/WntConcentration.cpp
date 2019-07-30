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
#include "WntConcentration.hpp"

/** Pointer to the single instance */
template<unsigned DIM>
WntConcentration<DIM>* WntConcentration<DIM>::mpInstance = nullptr;

template<unsigned DIM>
WntConcentration<DIM>* WntConcentration<DIM>::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new WntConcentration;
    }
    return mpInstance;
}

template<unsigned DIM>
WntConcentration<DIM>::WntConcentration()
    : mCryptLength(DOUBLE_UNSET),
      mLengthSet(false),
      mWntType(NONE),
      mpCellPopulation(nullptr),
      mTypeSet(false),
      mConstantWntValueForTesting(0),
      mUseConstantWntValueForTesting(false),
      mWntConcentrationParameter(1.0),
      mCryptProjectionParameterA(0.5),
      mCryptProjectionParameterB(2.0)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == nullptr);
}

template<unsigned DIM>
WntConcentration<DIM>::~WntConcentration()
{
}

template<unsigned DIM>
void WntConcentration<DIM>::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

template<unsigned DIM>
double WntConcentration<DIM>::GetWntLevel(CellPtr pCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell-cycle models without a cell population
    {
        return mConstantWntValueForTesting;
    }

    assert(mpCellPopulation!=nullptr);
    assert(mTypeSet);
    assert(mLengthSet);

    double height;

    if (mWntType == RADIAL)
    {
        double a = GetCryptProjectionParameterA();
        double b = GetCryptProjectionParameterB();
        height = a*pow(norm_2(mpCellPopulation->GetLocationOfCellCentre(pCell)), b);
    }
    else
    {
        height = mpCellPopulation->GetLocationOfCellCentre(pCell)[DIM-1];
    }

    return GetWntLevel(height);
}

template<unsigned DIM>
c_vector<double, DIM> WntConcentration<DIM>::GetWntGradient(CellPtr pCell)
{
    if (mUseConstantWntValueForTesting)  // to test a cell and cell-cycle models without a cell population
    {
        return zero_vector<double>(DIM);
    }
    assert(mpCellPopulation!=nullptr);
    assert(mTypeSet);
    assert(mLengthSet);

    c_vector<double, DIM> location_of_cell = mpCellPopulation->GetLocationOfCellCentre(pCell);

    return GetWntGradient(location_of_cell);
}

template<unsigned DIM>
void WntConcentration<DIM>::SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation)
{
    mpCellPopulation = &rCellPopulation;
}

template<unsigned DIM>
AbstractCellPopulation<DIM>& WntConcentration<DIM>::rGetCellPopulation()
{
    return *mpCellPopulation;
}

template<unsigned DIM>
double WntConcentration<DIM>::GetCryptLength()
{
    return mCryptLength;
}

template<unsigned DIM>
void WntConcentration<DIM>::SetCryptLength(double cryptLength)
{
    assert(cryptLength > 0.0);
    if (mLengthSet==true)
    {
        EXCEPTION("Destroy has not been called");
    }

    mCryptLength = cryptLength;
    mLengthSet = true;
}

template<unsigned DIM>
WntConcentrationType WntConcentration<DIM>::GetType()
{
    return mWntType;
}

template<unsigned DIM>
void WntConcentration<DIM>::SetType(WntConcentrationType type)
{
    if (mTypeSet==true)
    {
        EXCEPTION("Destroy has not been called");
    }
    mWntType = type;
    mTypeSet = true;
}

template<unsigned DIM>
double WntConcentration<DIM>::GetWntLevel(double height)
{
    if (mWntType == NONE)
    {
        return 0.0;
    }

    // Need to call SetCryptLength first
    assert(mLengthSet);

    double wnt_level = -1.0; // Test this is changed before leaving method.

    // The first type of Wnt concentration to try
    if (mWntType==LINEAR || mWntType==RADIAL)
    {
        if ((height >= -1e-9) && (height < mWntConcentrationParameter*GetCryptLength()))
        {
            wnt_level = 1.0 - height/(mWntConcentrationParameter*GetCryptLength());
        }
        else
        {
            wnt_level = 0.0;
        }
    }

    if (mWntType==EXPONENTIAL)
    {
        if ((height >= -1e-9) && (height < GetCryptLength()))
        {
            wnt_level = exp(-height/(GetCryptLength()*mWntConcentrationParameter));
        }
        else
        {
            wnt_level = 0.0;
        }
    }

    assert(wnt_level >= 0.0);

    return wnt_level;
}

template<unsigned DIM>
c_vector<double, DIM> WntConcentration<DIM>::GetWntGradient(c_vector<double, DIM>& rLocation)
{
    c_vector<double, DIM> wnt_gradient = zero_vector<double>(DIM);

    if (mWntType!=NONE)
    {
        if (mWntType==LINEAR)
        {
            if ((rLocation[DIM-1] >= -1e-9) && (rLocation[DIM-1] < mWntConcentrationParameter*GetCryptLength()))
            {
                wnt_gradient[DIM-1] = -1.0/(mWntConcentrationParameter*GetCryptLength());
            }
        }
        else if (mWntType==RADIAL) // RADIAL Wnt concentration
        {
            double a = GetCryptProjectionParameterA();
            double b = GetCryptProjectionParameterB();
            double r = norm_2(rLocation);
            double r_critical = pow(mWntConcentrationParameter*GetCryptLength()/a, 1.0/b);

            double dwdr = 0.0;

            if (r>=-1e-9 && r<r_critical)
            {
                dwdr = -mWntConcentrationParameter*GetCryptLength()*pow(r, b-1.0)/a;
            }

            for (unsigned i=0; i<DIM; i++)
            {
                wnt_gradient[i] = rLocation[i]*dwdr/r;
            }
        }
        else
        {
            EXCEPTION("No method to calculate gradient of this Wnt type");
        }
    }
    return wnt_gradient;
}

template<unsigned DIM>
bool WntConcentration<DIM>::IsWntSetUp()
{
    bool result = false;
    if (mTypeSet && mLengthSet && mpCellPopulation!=nullptr && mWntType!=NONE)
    {
        result = true;
    }
    return result;
}

template<unsigned DIM>
void WntConcentration<DIM>::SetConstantWntValueForTesting(double value)
{
    if (value < 0)
    {
        EXCEPTION("WntConcentration<DIM>::SetConstantWntValueForTesting - Wnt value for testing should be non-negative.\n");
    }
    mConstantWntValueForTesting = value;
    mUseConstantWntValueForTesting = true;
    if (!mTypeSet)
    {
        mWntType = NONE;
    }
}

template<unsigned DIM>
double WntConcentration<DIM>::GetWntConcentrationParameter()
{
    return mWntConcentrationParameter;
}

template<unsigned DIM>
void WntConcentration<DIM>::SetWntConcentrationParameter(double wntConcentrationParameter)
{
    assert(wntConcentrationParameter > 0.0);
    mWntConcentrationParameter = wntConcentrationParameter;
}

template<unsigned DIM>
double WntConcentration<DIM>::GetCryptProjectionParameterA()
{
    return mCryptProjectionParameterA;
}

template<unsigned DIM>
double WntConcentration<DIM>::GetCryptProjectionParameterB()
{
    return mCryptProjectionParameterB;
}

template<unsigned DIM>
void WntConcentration<DIM>::SetCryptProjectionParameterA(double cryptProjectionParameterA)
{
    assert(cryptProjectionParameterA >= 0.0);
    mCryptProjectionParameterA = cryptProjectionParameterA;
}

template<unsigned DIM>
void WntConcentration<DIM>::SetCryptProjectionParameterB(double cryptProjectionParameterB)
{
    assert(cryptProjectionParameterB >= 0.0);
    mCryptProjectionParameterB = cryptProjectionParameterB;
}

// Explicit instantiation
template class WntConcentration<1>;
template class WntConcentration<2>;
template class WntConcentration<3>;
