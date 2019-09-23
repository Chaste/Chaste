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

#ifndef SINGLETRACEOUTPUTMODIFIER_HPP_
#define SINGLETRACEOUTPUTMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOutputModifier.hpp"
#include "OutputFileHandler.hpp"

/**
 * Provide the trace for a transmembrane potential at a single mode of the mesh.
 * File updated during the simulation.
 *
 * WARNING:  If you checkpoint this class then the file output will not be saved in the checkpoint.
 *           The state of the output file will be unchanged (up to the most recent file flush) but
 *           may be overwritten by a restarted simulation.   You will need to manually move the output file
 *           before restarting the simulation and then merge the two files together at a later point.
 */
class SingleTraceOutputModifier : public AbstractOutputModifier
{
private:
    /** For testing */
    friend class TestMonodomainProblem;
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /** The global index of the node for which the trace is to be made.
     *  This is the index *in memory at solve time*.  If you are running in parallel and you want a
     *  specific index/location in your mesh then you will need to look in the mesh permutation.
     */
    unsigned mGlobalIndex;
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
        // This one would need re-calculating, so we don't archive it... archive & mLocalIndex;
    }

    /** Private constructor that resets process-specific data, for archiving */
    SingleTraceOutputModifier()
        : mLocalIndex(UINT_MAX),
          mFileStream(NULL)
    {}

public:
    /**
     * Constructor
     *
     * @param globalIndex The global index of the node which is to be output.
     *  This is the index *in the original mesh*.  If you are running in parallel your index in the original
     *  mesh will be automatically converted to the runtime index using the node permutation in the mesh.
     *  That is, is you ask for node 5 in the original mesh then you will automatically get
     *  \code
     *     unsigned new_index_for_5 = mesh.rGetNodePermutation()[5];
     *  \endcode
     *  which is the runtime index of the same node in space.
     *
     * @param rFilename  The file which is eventually produced by this modifier
     * @param flushTime The simulation time between manual file flushes (if required)
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
     * Initialise the modifier (open a file or make some memory) when the solve loop is starting
     *
     * Note the problem passes parameters in a non-templated fashion in order to keep the interface as lightweight as
     * possible.  That is, it might have been slicker to pass in the mesh but that would require multiple templates.
     * @param pVectorFactory  The vector factory which is associated with the calling problem's mesh
     * @param rNodePermutation The permutation associated with the calling problem's mesh (when running with parallel partitioning)
     */
    virtual void InitialiseAtStart(DistributedVectorFactory* pVectorFactory, const std::vector<unsigned>& rNodePermutation);

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
