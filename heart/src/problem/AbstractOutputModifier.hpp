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


#ifndef ABSTRACTOUTPUTMODIFIER_HPP_
#define ABSTRACTOUTPUTMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/string.hpp>
#include "ClassIsAbstract.hpp"

#include <string>
#include "AbstractTetrahedralMesh.hpp"

/**
 * A plug-in class for on-the-fly output.  This is designed so that a user can insert something in
 * order to monitor the progress of a simulation or to produce "post processed" output during the simulation.
 */
class AbstractOutputModifier
{
private:
    /** For testing */
    friend class TestMonodomainProblem;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive this modifier.  Just calls the base class version.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mFilename;
        archive & mFlushTime;
    }

protected:
    /** Constructor that does nothing, for archiving */
    AbstractOutputModifier()
        : mFilename(), mFlushTime(0.0)
    {}

    std::string mFilename; /**<The file which is eventually produced by this modifier*/

    /** Simulation time period between flushes to disk */
    double mFlushTime;

public:
    /**
     * Standard construction method contains only the name of the file which this simulation modifier should produce.
     *
     * Note that construction occurs before the Solve() loop and probably before initialisation.  This means that it
     * will not be possible to view certain data (e.g the mesh) at the time of construction
     *
     * Note the problem passes parameters in a non-templated fashion in order to keep the interface as lightweight as
     * possible.
     * @param rFilename  The file which is eventually produced by this modifier
     * @param flushTime  The (simulation) time between flushing the file to disk. Default (0) means don't flush.
     *
     */
    AbstractOutputModifier(const std::string& rFilename, double flushTime=0.0)
        : mFilename(rFilename),
          mFlushTime(flushTime)
    {}
    /**
     * Destructor should be overridden when necessary
     */
    virtual ~AbstractOutputModifier()
    {

    }

    /**
     * Initialise the modifier (open a file or make some memory) when the solve loop is starting
     *
     * Note the problem passes parameters in a non-templated fashion in order to keep the interface as lightweight as
     * possible.  That is, it might have been slicker to pass in the mesh but that would require multiple templates.
     * @param pVectorFactory  The vector factory which is associated with the calling problem's mesh
     * @param rNodePermutation The permutation associated with the calling problem's mesh (when running with parallel partitioning)
     */
    virtual void InitialiseAtStart(DistributedVectorFactory* pVectorFactory, const std::vector<unsigned>& rNodePermutation)=0;

    /**
     * Finalise the modifier (close a file or dump the calculation to disk)
     */
    virtual void FinaliseAtEnd()=0;

    /**
     * Process a solution time-step (dump a small line to file or compute the latest activation times)
     * @param time  The current simulation time
     * @param solution  A working copy of the solution at the current time-step.  This is the PETSc vector which is distributed across the processes.
     * @param problemDim  The calling problem dimension. Used here to avoid probing the size of the solution vector
     */
    virtual void ProcessSolutionAtTimeStep(double time, Vec solution, unsigned problemDim)=0;
};

CLASS_IS_ABSTRACT(AbstractOutputModifier)

#endif /* ABSTRACTOUTPUTMODIFIER_HPP_ */
