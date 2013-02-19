/*

Copyright (c) 2005-2013, University of Oxford.
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

// Serialisation headers.
#include "CheckpointArchiveTypes.hpp"

#include "ObjectCommunicator.hpp"

#include "PetscTools.hpp"
#include "ClassOfSimpleVariables.hpp"
#include <sstream>
#include <string>
#include <cstring>
#include <set>
#include <assert.h>

template<typename CLASS>
ObjectCommunicator<CLASS>::ObjectCommunicator()
{    
}

template<typename CLASS>
void ObjectCommunicator<CLASS>::SendObject(CLASS* const pObject, unsigned destinationProcess, unsigned tag)
{
    // Create an output archive
    std::stringstream ss;
    boost::archive::text_oarchive output_arch(ss);

    output_arch << pObject;

    std::string send_msg = ss.str();

    // Get + send string length
    unsigned string_length = send_msg.size();
    MPI_Send(&string_length, 1, MPI_UNSIGNED, destinationProcess, tag, PETSC_COMM_WORLD);

    // Pack it into a string of chars
    char arc_string[string_length];
    strcpy( arc_string, send_msg.c_str());

    // Send archive data
    MPI_Send(arc_string, string_length, MPI_CHAR, destinationProcess, tag, PETSC_COMM_WORLD);
}

template<typename CLASS>
CLASS* ObjectCommunicator<CLASS>::RecvObject(unsigned sourceProcess, unsigned tag, MPI_Status& status)
{
    unsigned string_length = 0;
    MPI_Recv(&string_length, 1, MPI_UNSIGNED, sourceProcess, tag, PETSC_COMM_WORLD, &status);

    char recv_string[string_length];
    MPI_Recv(&recv_string, string_length, MPI_CHAR, sourceProcess , tag, PETSC_COMM_WORLD, &status);

    // Make it into a string
    std::string rcv_msg;
    for(unsigned i=0; i<string_length; i++)
    {
        rcv_msg.push_back(recv_string[i]);
    }

    // Save
    std::stringstream ss;
    ss << rcv_msg;

    CLASS* p_recv_object;
    boost::archive::text_iarchive input_arch(ss);

    input_arch >> p_recv_object;

    return p_recv_object;
}

template<typename CLASS>
CLASS* ObjectCommunicator<CLASS>::SendRecvObject(CLASS* const pSendObject, unsigned destinationProcess, unsigned sendTag, unsigned sourceProcess, unsigned sourceTag, MPI_Status& status)
{
    // Create an output archive
    std::stringstream oss;
    boost::archive::text_oarchive output_arch(oss);

    output_arch << pSendObject;

    std::string send_msg = oss.str();

    // Get + send string length
    unsigned send_string_length = send_msg.size();
    unsigned recv_string_length;

    char arc_string[send_string_length];
    for(unsigned i=0; i<send_string_length; i++)
    {
        arc_string[i] = send_msg[i];
    }

    MPI_Sendrecv(&send_string_length, 1, MPI_UNSIGNED, destinationProcess, sendTag, &recv_string_length, 1, MPI_UNSIGNED, sourceProcess, sourceTag, PETSC_COMM_WORLD, &status);

    char recv_string[recv_string_length];

    // Send archive data
    MPI_Sendrecv(arc_string, send_string_length, MPI_CHAR, destinationProcess, sendTag, &recv_string, recv_string_length, MPI_BYTE, sourceProcess, sourceTag, PETSC_COMM_WORLD, &status);

    // Make it into a string
    std::string rcv_msg;
    for(unsigned i=0; i<recv_string_length; i++)
    {
        rcv_msg.push_back(recv_string[i]);
    }

    // Save
    std::istringstream iss(rcv_msg);

    CLASS* p_recv_object;
    boost::archive::text_iarchive input_arch(iss);

    input_arch >> p_recv_object;

    return p_recv_object;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////
template class ObjectCommunicator<ClassOfSimpleVariables>;

