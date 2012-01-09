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

#include "OdeSolution.hpp"

#include <sstream>

#include "ColumnDataWriter.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "VectorHelperFunctions.hpp"

OdeSolution::OdeSolution()
    : mNumberOfTimeSteps(0u),
      mpOdeSystemInformation()
{
}


unsigned OdeSolution::GetNumberOfTimeSteps() const
{
    return mNumberOfTimeSteps;
}


void OdeSolution::SetNumberOfTimeSteps(unsigned numTimeSteps)
{
    mNumberOfTimeSteps = numTimeSteps;
    mTimes.reserve(numTimeSteps+1);
    mSolutions.reserve(numTimeSteps);
}


void OdeSolution::SetOdeSystemInformation(boost::shared_ptr<const AbstractOdeSystemInformation> pOdeSystemInfo)
{
    mpOdeSystemInformation = pOdeSystemInfo;
}

std::vector<double> OdeSolution::GetVariableAtIndex(unsigned index) const
{
    std::vector<double> answer;
    answer.reserve(mTimes.size());
    double temp_number;
    for (unsigned i=0; i< mTimes.size(); ++i)
    {
        if (index < mSolutions[0].size())
        {
            temp_number = mSolutions[i][index];
        }
        else
        {
            unsigned offset = mSolutions[0].size();
            if (index - offset < mParameters.size())
            {
                temp_number = mParameters[index - offset];
            }
            else
            {
                offset += mParameters.size();
                if (index - offset < mDerivedQuantities[0].size())
                {
                    temp_number = mDerivedQuantities[i][index - offset];
                }
                else
                {
                    EXCEPTION("Invalid index passed to ""GetVariableAtIndex()"".");
                }
            }
        }
        answer.push_back(temp_number);
    }
    return answer;
}

std::vector<double> OdeSolution::GetAnyVariable(const std::string& rName) const
{
    return GetVariableAtIndex(mpOdeSystemInformation->GetAnyVariableIndex(rName));
}

std::vector<double>& OdeSolution::rGetTimes()
{
    return mTimes;
}

const std::vector<double>& OdeSolution::rGetTimes() const
{
    return mTimes;
}

std::vector<std::vector<double> >& OdeSolution::rGetSolutions()
{
    return mSolutions;
}

const std::vector<std::vector<double> >& OdeSolution::rGetSolutions() const
{
    return mSolutions;
}


template<typename VECTOR>
void OdeSolution::CalculateDerivedQuantitiesAndParameters(AbstractParameterisedSystem<VECTOR>* pOdeSystem)
{
    assert(pOdeSystem->GetSystemInformation() == mpOdeSystemInformation); // Just in case...
    rGetParameters(pOdeSystem);
    rGetDerivedQuantities(pOdeSystem);
}


template<typename VECTOR>
std::vector<double>& OdeSolution::rGetParameters(AbstractParameterisedSystem<VECTOR>* pOdeSystem)
{
    mParameters.clear();
    const unsigned num_params = pOdeSystem->GetNumberOfParameters();
    if (num_params > 0)
    {
        mParameters.reserve(num_params);
        for (unsigned i=0; i<num_params; ++i)
        {
            mParameters.push_back(pOdeSystem->GetParameter(i));
        }
    }
    return mParameters;
}


std::vector<std::vector<double> >& OdeSolution::rGetDerivedQuantities(AbstractParameterisedSystem<std::vector<double> >* pOdeSystem)
{
    assert(pOdeSystem != NULL);
    if (mDerivedQuantities.empty() && pOdeSystem->GetNumberOfDerivedQuantities() > 0)
    {
        assert(mTimes.size() == mSolutions.size()); // Paranoia
        mDerivedQuantities.reserve(mTimes.size());
        for (unsigned i=0; i<mTimes.size(); i++)
        {
            mDerivedQuantities.push_back(pOdeSystem->ComputeDerivedQuantities(mTimes[i], mSolutions[i]));
        }
    }
    return mDerivedQuantities;
}

#ifdef CHASTE_CVODE
std::vector<std::vector<double> >& OdeSolution::rGetDerivedQuantities(AbstractParameterisedSystem<N_Vector>* pOdeSystem)
{
    assert(pOdeSystem != NULL);
    if (mDerivedQuantities.empty() && pOdeSystem->GetNumberOfDerivedQuantities() > 0)
    {
        const unsigned num_solutions = mSolutions.size();
        assert(mTimes.size() == num_solutions); // Paranoia
        mDerivedQuantities.resize(mTimes.size());
        N_Vector state_vars = num_solutions > 0 ? N_VNew_Serial(mSolutions[0].size()) : NULL;
        for (unsigned i=0; i<num_solutions; i++)
        {
            CopyFromStdVector(mSolutions[i], state_vars);
            N_Vector dqs = pOdeSystem->ComputeDerivedQuantities(mTimes[i], state_vars);
            CopyToStdVector(dqs, mDerivedQuantities[i]);
            DeleteVector(dqs);
        }
        DeleteVector(state_vars);
    }
    assert(mDerivedQuantities.size()==mTimes.size());
    return mDerivedQuantities;
}
#endif // CHASTE_CVODE


void OdeSolution::WriteToFile(std::string directoryName,
                              std::string baseResultsFilename,
                              std::string timeUnits,
                              unsigned stepsPerRow,
                              bool cleanDirectory,
                              unsigned precision,
                              bool includeDerivedQuantities)
{
    assert(stepsPerRow > 0);
    assert(mTimes.size() > 0);
    assert(mTimes.size() == mSolutions.size());
    assert(mpOdeSystemInformation.get() != NULL);
    if (mpOdeSystemInformation->GetNumberOfParameters()==0 && mpOdeSystemInformation->GetNumberOfDerivedQuantities() == 0)
    {
        includeDerivedQuantities = false;
    }

    if (includeDerivedQuantities)
    {
        if ((mDerivedQuantities.empty() || mDerivedQuantities.size()!=mTimes.size()) && mParameters.empty())
        {
            EXCEPTION("You must first call ""CalculateDerivedQuantitiesAndParameters()"" in order to write derived quantities.");
        }
    }

    // Write data to a file using ColumnDataWriter
    ColumnDataWriter writer(directoryName, baseResultsFilename, cleanDirectory, precision);

    if (!PetscTools::AmMaster())
    {
        //Only the master actually writes to file
        return;
    }

    int time_var_id = writer.DefineUnlimitedDimension("Time", timeUnits);

    // Either: the ODE system should have no names&units defined, or it should
    // the same number as the number of solutions per timestep.
    assert(  mpOdeSystemInformation->rGetStateVariableNames().size()==0 ||
            (mpOdeSystemInformation->rGetStateVariableNames().size()==mSolutions[0].size()) );

    unsigned num_vars = mSolutions[0].size();
    unsigned num_params = mpOdeSystemInformation->GetNumberOfParameters();
    unsigned num_derived_quantities = mpOdeSystemInformation->GetNumberOfDerivedQuantities();

    std::vector<int> var_ids;
    var_ids.reserve(num_vars);
    if (mpOdeSystemInformation->rGetStateVariableNames().size() > 0)
    {
        for (unsigned i=0; i<num_vars; i++)
        {
            var_ids.push_back(writer.DefineVariable(mpOdeSystemInformation->rGetStateVariableNames()[i],
                                                    mpOdeSystemInformation->rGetStateVariableUnits()[i]));
        }
    }
    else
    {
        for (unsigned i=0; i<num_vars; i++)
        {
            std::stringstream string_stream;
            string_stream << "var_" << i;
            var_ids.push_back(writer.DefineVariable(string_stream.str(), ""));
        }
    }

    if (includeDerivedQuantities)
    {
        var_ids.reserve(num_vars + num_params + num_derived_quantities);
        for (unsigned i=0; i<num_params; ++i)
        {
            var_ids.push_back(writer.DefineVariable(mpOdeSystemInformation->rGetParameterNames()[i],
                                                    mpOdeSystemInformation->rGetParameterUnits()[i]));
        }
        for (unsigned i=0; i<num_derived_quantities; i++)
        {
            var_ids.push_back(writer.DefineVariable(mpOdeSystemInformation->rGetDerivedQuantityNames()[i],
                                                    mpOdeSystemInformation->rGetDerivedQuantityUnits()[i]));
        }
    }

    if (mSolverName != "")
    {
        writer.SetCommentForInfoFile("ODE SOLVER: " + mSolverName);
    }

    writer.EndDefineMode();

    for (unsigned i=0; i<mSolutions.size(); i+=stepsPerRow)
    {
        writer.PutVariable(time_var_id, mTimes[i]);
        for (unsigned j=0; j<num_vars; j++)
        {
            writer.PutVariable(var_ids[j], mSolutions[i][j]);
        }
        if (includeDerivedQuantities)
        {
            for (unsigned j=0; j<num_params; ++j)
            {
                writer.PutVariable(var_ids[j+num_vars], mParameters[j]);
            }
            for (unsigned j=0; j<num_derived_quantities; j++)
            {
                writer.PutVariable(var_ids[j+num_params+num_vars], mDerivedQuantities[i][j]);
            }
        }
        writer.AdvanceAlongUnlimitedDimension();
    }
    writer.Close();
}

//
// Explicit instantiation
//

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template std::vector<double>& OdeSolution::rGetParameters(AbstractParameterisedSystem<std::vector<double> >* pOdeSystem);
template void OdeSolution::CalculateDerivedQuantitiesAndParameters(AbstractParameterisedSystem<std::vector<double> >* pOdeSystem);

#ifdef CHASTE_CVODE
template std::vector<double>& OdeSolution::rGetParameters(AbstractParameterisedSystem<N_Vector>* pOdeSystem);
template void OdeSolution::CalculateDerivedQuantitiesAndParameters(AbstractParameterisedSystem<N_Vector>* pOdeSystem);
#endif // CHASTE_CVODE

/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
