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

#include "DynamicVentilationProblem.hpp"
#include "ProgressReporter.hpp"

DynamicVentilationProblem::DynamicVentilationProblem(AbstractAcinarUnitFactory* pAcinarFactory,
                                                     const std::string& rMeshDirFilePath,
                                                     unsigned rootIndex) : mpAcinarFactory(pAcinarFactory),
                                                                           mVentilationProblem(rMeshDirFilePath, rootIndex),
                                                                           mrMesh(mVentilationProblem.rGetMesh()),
                                                                           mDt(0.01),
                                                                           mSamplingTimeStepMultiple(1u),
                                                                           mCurrentTime(0.0),
                                                                           mRootIndex(rootIndex),
                                                                           mWriteVtkOutput(false)
{
    mVentilationProblem.SetOutflowPressure(0.0);

    mpAcinarFactory->SetMesh(&mrMesh);

    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mrMesh.GetBoundaryNodeIteratorBegin();
                        iter != mrMesh.GetBoundaryNodeIteratorEnd();
                        ++iter )
    {
        if ((*iter)->GetIndex() != rootIndex)
        {
            mAcinarMap[(*iter)->GetIndex()] = mpAcinarFactory->CreateAcinarUnitForNode((*iter));
        }
    }
}

DynamicVentilationProblem::~DynamicVentilationProblem()
{
    for (std::map<unsigned, AbstractAcinarUnit*>::iterator iter = mAcinarMap.begin();
         iter != mAcinarMap.end();
         ++iter )
    {
        delete iter->second;
    }
}

MatrixVentilationProblem& DynamicVentilationProblem::rGetMatrixVentilationProblem()
{
    return mVentilationProblem;
}

std::map<unsigned, AbstractAcinarUnit*>& DynamicVentilationProblem::rGetAcinarUnitMap()
{
    return mAcinarMap;
}

void DynamicVentilationProblem::SetTimeStep(double timeStep)
{
    mDt = timeStep;
}

void DynamicVentilationProblem::SetSamplingTimeStepMultiple(unsigned timeStep)
{
    mSamplingTimeStepMultiple = timeStep;
}

void DynamicVentilationProblem::SetEndTime(double time)
{
    mEndTime = time;
}

void DynamicVentilationProblem::SetOutputDirectory(const std::string directory)
{
    mOutputDirectory = directory;
}

void DynamicVentilationProblem::SetOutputFilenamePrefix(const std::string prefix)
{
    mOutputFileNamePrefix = prefix;
}

void DynamicVentilationProblem::SetWriteVtkOutput(bool writeVtkOutput)
{
    mWriteVtkOutput = writeVtkOutput;
}

void DynamicVentilationProblem::Solve()
{
    TimeStepper time_stepper(mCurrentTime, mEndTime, mDt);

#ifdef CHASTE_VTK
    VtkMeshWriter<1, 3> vtk_writer(mOutputDirectory, mOutputFileNamePrefix, false);
#endif //CHASTE_VTK

    ProgressReporter progress_reporter(mOutputDirectory, mCurrentTime, mEndTime);

    std::vector<double> pressures(mrMesh.GetNumNodes(), -1);
    std::vector<double> fluxes(mrMesh.GetNumNodes() - 1, -1);
    std::vector<double> volumes(mrMesh.GetNumNodes(), -1);

    while (!time_stepper.IsTimeAtEnd())
    {
        //Solve coupled problem
        for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mrMesh.GetBoundaryNodeIteratorBegin();
                 iter != mrMesh.GetBoundaryNodeIteratorEnd();
                 ++iter )
        {
            if ((*iter)->GetIndex() != mRootIndex)
            {
                double pleural_pressure =  mpAcinarFactory->GetPleuralPressureForNode(time_stepper.GetNextTime(), (*iter));
                mAcinarMap[(*iter)->GetIndex()]->SetPleuralPressure(pleural_pressure);
                mAcinarMap[(*iter)->GetIndex()]->ComputeExceptFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

                mVentilationProblem.SetPressureAtBoundaryNode(*(*iter), mAcinarMap[(*iter)->GetIndex()]->GetAirwayPressure());
            }
        }

        mVentilationProblem.Solve();
        mVentilationProblem.GetSolutionAsFluxesAndPressures(fluxes, pressures);

        for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mrMesh.GetBoundaryNodeIteratorBegin();
                            iter != mrMesh.GetBoundaryNodeIteratorEnd();
                            ++iter )
        {
            if ((*iter)->GetIndex() != 0u)
            {
                unsigned boundary_element_index = (*(*iter)->rGetContainingElementIndices().begin());

                mAcinarMap[(*iter)->GetIndex()]->SetFlow(fluxes[boundary_element_index]);

                double resistance = 0.0;
                if (fluxes[boundary_element_index] != 0.0)
                {
                    resistance = std::fabs(pressures[(*iter)->GetIndex()]/fluxes[boundary_element_index]);
                }
                mAcinarMap[(*iter)->GetIndex()]->SetTerminalBronchioleResistance(resistance);
                mAcinarMap[(*iter)->GetIndex()]->UpdateFlow(time_stepper.GetTime(), time_stepper.GetNextTime());
            }
        }

        if ((time_stepper.GetTotalTimeStepsTaken() % mSamplingTimeStepMultiple) == 0u)
        {
            progress_reporter.Update(time_stepper.GetNextTime());

#ifdef CHASTE_VTK
            if (mWriteVtkOutput)
            {
                std::ostringstream suffix_name;
                suffix_name <<  "_" << std::setw(6) << std::setfill('0') << time_stepper.GetTotalTimeStepsTaken()/mSamplingTimeStepMultiple;

                vtk_writer.AddCellData("Flux"+suffix_name.str(), fluxes);
                vtk_writer.AddPointData("Pressure"+suffix_name.str(), pressures);


                for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mrMesh.GetBoundaryNodeIteratorBegin();
                                                iter != mrMesh.GetBoundaryNodeIteratorEnd();
                                                ++iter )
                {
                    if ((*iter)->GetIndex() != mRootIndex)
                    {
                        volumes[(*iter)->GetIndex()] = mAcinarMap[(*iter)->GetIndex()]->GetVolume();
                    }
                }

                vtk_writer.AddPointData("Volume"+suffix_name.str(), volumes);
            }
#endif //CHASTE_VTK
        }

        mCurrentTime = time_stepper.GetNextTime();
        time_stepper.AdvanceOneTimeStep();
    }

#ifdef CHASTE_VTK
    if (mWriteVtkOutput)
    {
        vtk_writer.WriteFilesUsingMesh(mrMesh);
    }
#endif //CHASTE_VTK
}
