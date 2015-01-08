/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef SINGLETRACEOUTPUTMODIFIER_HPP_
#define SINGLETRACEOUTPUTMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOutputModifier.hpp"
#include "OutputFileHandler.hpp"

/**
 * Provide the trace for a transmembrane potential at a single mode of the mesh.
 * File updated during the simulation.
 */
class SingleTraceOutputModifier : public AbstractOutputModifier
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;

    unsigned mGlobalIndex; /**< The global index of the node for which the trace is to be made.*/
    unsigned mLocalIndex; /**< The local index of the node for which the trace is to be made - set to UINT_MAX if the node is not local to the process*/
    out_stream mFileStream; /**< Output file stream (remains open during solve).*/

    friend class TestOutputModifiers;

    /**
     * Archive the output modifier, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractOutputModifier>(*this);
        archive & mGlobalIndex;
        // This one need re-calculating... archive & mLocalIndex;
    }

    /** Private constructor that resets process-specific data, for archiving */
    SingleTraceOutputModifier()
        :mLocalIndex(UINT_MAX),
         mFileStream(NULL){};

public:
    /**
     * Constructor
     *
     * @param globalIndex The global index of the node which is to be output (assumes no permutation)
     * @param rFilename  The file which is eventually produced by this modifier
     *
     */
    SingleTraceOutputModifier(const std::string& rFilename, unsigned globalIndex, double flushTime=0.0)
        : AbstractOutputModifier(rFilename, flushTime),
          mGlobalIndex(globalIndex),
          mLocalIndex(UINT_MAX),
          mFileStream(NULL)
    {
    }
    /**
     * Initialise the modifier (open the file) when the solve loop is starting.
     *
     * Note the problem passes parameters in a non-templated fashion in order to keep the interface as lightweight as
     * possible.
     * @param pVectorFactory  The vector factory which is associated with the calling problem's mesh
     */
    virtual void InitialiseAtStart(DistributedVectorFactory* pVectorFactory);

    /**
     * Finalise the modifier (close the file)
     */
    virtual void FinaliseAtEnd();

    /**
     * Process a solution time-step (dump a small line to file)
     * @param time  The current simulation time
     * @param solution  A working copy of the solution at the current time-step.  This is the PETSc vector which is distributed across the processes.
     * @param problemDim  The calling problem dimension. Used here to avoid probing the size of the solution vector
     */
    virtual void ProcessSolutionAtTimeStep(double time, Vec solution, unsigned problemDim);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SingleTraceOutputModifier)

#endif // SINGLETRACEOUTPUTMODIFIER_HPP_
