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

#ifndef DYNAMICVENTILATIONPROBLEM_HPP_
#define DYNAMICVENTILATIONPROBLEM_HPP_

#include "AbstractVentilationProblem.hpp"
#include "TimeStepper.hpp"
#ifdef CHASTE_VTK
#include "VtkMeshWriter.hpp"
#endif //CHASTE_VTK
#include "AbstractAcinarUnitFactory.hpp"
#include "AbstractAcinarUnit.hpp"
#include "MatrixVentilationProblem.hpp"
#include <map>

/**
 * A class for solving dynamic one-dimensional lung ventilation problems in which each terminal of
 * the conducting airway tree is joined to an acinar "balloon" model.
 */
class DynamicVentilationProblem
{
public:

    /**
     * Create a new dynamic ventilation problem.
     *
     * @param rAcinarUnitFactory Factory class to create acinar units.
     * @param rMeshDirFilePath Path to mesh files in triangles/tetgen format
     * @param rootIndex Index of the node at the start of the trachea
     */
    DynamicVentilationProblem(AbstractAcinarUnitFactory* pAcinarFactory,
                              const std::string& rMeshDirFilePath,
                              unsigned rootIndex);

    /**
     * Destructor
     */
    ~DynamicVentilationProblem();

    /**
     * @return Reference to the ventilation problem
     */
    MatrixVentilationProblem& rGetMatrixVentilationProblem();

    /**
     * @return Reference to a map of acinar units used by this problem
     */
    std::map<unsigned, AbstractAcinarUnit*>& rGetAcinarUnitMap();

    /**
     * Set the size of the global time step
     *
     * @param timeStep The timestep used by the simulation
     */
    void SetTimeStep(double timeStep);

    /**
     * Set the ratio of the number of actual time steps to the number of time steps at which results are written to file. Default value is set to 1 by the constructor.
     *
     * @param multiple Number of discrete timesteps between output being written.
     */
    void SetSamplingTimeStepMultiple(unsigned multiple);

    /**
     * Set the finishing time of the simulation
     *
     * @param time The time to end the simulation
     */
    void SetEndTime(double time);

    /**
     * Specify the directory to write output to
     *
     * @param directory The directory to write to
     */
    void SetOutputDirectory(const std::string directory);

    /**
     * Specify the root of the output file name
     *
     * @param prefix Prefix of the output file name
     */
    void SetOutputFilenamePrefix(const std::string prefix);

    /**
     * Solves up to the next finishing time
     */
    void Solve();

    /**
     * Tell the solver to output VTK files for visualisation or not
     *
     * @param writeVtkOutput Boolean indicating whether VTK output should be written or not
     */
    void SetWriteVtkOutput(bool writeVtkOutput = true);

private:
    /**
     * Acinar factory
     */
    AbstractAcinarUnitFactory* mpAcinarFactory;

    /**
     * Ventilation solver
     */
    MatrixVentilationProblem mVentilationProblem;

    /**
     * Map between boundary node ids and acinar balloon models.
     */
    std::map<unsigned, AbstractAcinarUnit*> mAcinarMap;

    /**
     * The airway tree mesh
     */
    TetrahedralMesh<1,3>& mrMesh;

    /**
     * Time step size
     */
    double mDt;

    /**
     * Number of timesteps before output is written to file
     */
    unsigned mSamplingTimeStepMultiple;

    /**
     * The current simulation time
     */
    double mCurrentTime;

    /**
     * The time at which the simulation will end
     */
    double mEndTime;

    /**
     * The index of the node at the entrance to the trachea.
     */
    unsigned mRootIndex;

    /**
     * Directory where output will be written to
     */
    std::string mOutputDirectory;

    /**
     * Prefix of output file name.
     */
    std::string mOutputFileNamePrefix;

    /**
     * Boolean to determine whether to write out VTK output or not.
     */
    bool mWriteVtkOutput;
};

#endif /* DYNAMICVENTILATIONPROBLEM_HPP_ */
