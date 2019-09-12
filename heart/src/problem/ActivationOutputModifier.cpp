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

#include "ActivationOutputModifier.hpp"
#include "OutputFileHandler.hpp"
#include "HeartConfig.hpp"

void ActivationOutputModifier::InitialiseAtStart(DistributedVectorFactory* pVectorFactory, const std::vector<unsigned>& rNodePermutation)
{
    mLocalSize = pVectorFactory->GetLocalOwnership();
    mFirstActivitationTimes.resize(mLocalSize, -1.0);
    mFirstRecoveryTimes.resize(mLocalSize, -1.0);
    mSecondActivitationTimes.resize(mLocalSize, -1.0);
    mSecondRecoveryTimes.resize(mLocalSize, -1.0);
}

void ActivationOutputModifier::FinaliseAtEnd()
{
    //Dump out all data in a round-robin fashion
    OutputFileHandler output_handler(HeartConfig::Instance()->GetOutputDirectory(), false);

    // Belt and braces
    std::stringstream filepath_process_specific;
    filepath_process_specific << mFilename << "." << PetscTools::GetMyRank();
    out_stream file_stream_process_specific = output_handler.OpenOutputFile(filepath_process_specific.str().c_str());
    for (unsigned i=0; i<mLocalSize; i++)
    {
        (*file_stream_process_specific) << mFirstActivitationTimes[i] <<",\t"
                << mFirstRecoveryTimes[i] <<",\t"
                << mSecondActivitationTimes[i] <<",\t"
                << mSecondRecoveryTimes[i] <<"\n";
    }
    file_stream_process_specific->close();

    PetscTools::BeginRoundRobin();
    {
        out_stream file_stream = out_stream(NULL);
        // Open the file as new or append
        if (PetscTools::AmMaster())
        {
            file_stream = output_handler.OpenOutputFile(mFilename);
        }
        else
        {
            file_stream = output_handler.OpenOutputFile(mFilename, std::ios::app);
        }
        for (unsigned i=0; i<mLocalSize; i++)
        {
            (*file_stream) << mFirstActivitationTimes[i] <<",\t"
                    << mFirstRecoveryTimes[i] <<",\t"
                    << mSecondActivitationTimes[i] <<",\t"
                    << mSecondRecoveryTimes[i] <<"\n";
        }
        file_stream->close();
    }
    PetscTools::EndRoundRobin();
}

void ActivationOutputModifier::ProcessSolutionAtTimeStep(double time, Vec solution, unsigned problemDim)
{
    double* p_solution;
    VecGetArray(solution, &p_solution);
    for (unsigned local_index=0; local_index < mLocalSize; local_index++)
    {
        double v = p_solution[local_index*problemDim];
        if (v > mThreshold && mFirstActivitationTimes[local_index] < 0.0)
        {
            mFirstActivitationTimes[local_index] = time;
        }
// LCOV_EXCL_START // Continuous tests are not long enough to allow recover
        else if (mFirstActivitationTimes[local_index] >= 0.0 && mFirstRecoveryTimes[local_index] < 0.0 && v < mThreshold)
        {
            ///\todo #2570 Add to a longer running test
            mFirstRecoveryTimes[local_index] = time;
        }
        else if (mFirstRecoveryTimes[local_index] >= 0.0 && mSecondActivitationTimes[local_index] < 0.0 && v > mThreshold)
        {
            mSecondActivitationTimes[local_index] = time;
        }
        else if (mSecondActivitationTimes[local_index] >= 0.0 && mSecondRecoveryTimes[local_index] < 0.0 && v < mThreshold)
        {
            mSecondRecoveryTimes[local_index] = time;
        }
// LCOV_EXCL_STOP // Continuous tests are not long enough to allow recover
    }
    VecRestoreArray(solution, &p_solution);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ActivationOutputModifier)
