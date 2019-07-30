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

#ifndef _OBJECTCOMMUNICATOR_HPP_
#define _OBJECTCOMMUNICATOR_HPP_

// Serialisation headers - must come first
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "PetscTools.hpp" // For MPI methods
#include "Exception.hpp"

const unsigned MAX_BUFFER_SIZE = 1000000;

/**
 * This is a helper class to enable classes that can be serialized to be sent using
 * PetSc MPI communication. The object is serialized in to a string of characters, and then
 * de-serialized on the receive process.
 */
template<typename CLASS>
class ObjectCommunicator
{
private:

    /** A buffer for use in asynchronous communication */
    char* mRecvBuffer;

    /** A group of buffers for use in asynchronous communication.  There's one for each process so that
     * a non-blocking send request won't accidentally overwrite a message which is actively being communicated
     * to another remote process */
    std::vector<char* > mSendBuffer;

    /** A group of strings for use in asynchronous communication.  There's one for each process so that
     * a non-blocking send request won't accidentally overwrite a message which is actively being communicated
     * to another remote process */
    std::vector<std::string> mSendString;

    /** The size of a string we are waiting for in an asynchronous receive */
    unsigned mRecvBufferLength;

    /** The size of a string we are sending */
    unsigned mSendBufferLength;

    /** An MPI_Request used in MPI_Irecv */
    MPI_Request mMpiRequest;

    /** A flag, used as a lock to ensure that we don't accidentally start overwriting the above buffer, #mRecvBuffer */
    bool mIsWriting;

public:

    /**
     * Default constructor.
     */
    ObjectCommunicator();

    /**
     * Send an object.
     *
     * @param pObject A pointer to the object to be sent
     * @param destinationProcess the index of the process to send the data to
     * @param tag a unique identifier tag for this communication
     */
    void SendObject(boost::shared_ptr<CLASS> const pObject, unsigned destinationProcess, unsigned tag);

    /**
     * Send an object.
     *
     * @param pObject A pointer to the object to be sent
     * @param destinationProcess the index of the process to send the data to
     * @param tag a unique identifier tag for this communication
     */
    void ISendObject(boost::shared_ptr<CLASS> const pObject, unsigned destinationProcess, unsigned tag);

    /**
     * Receive an object
     *
     * @param sourceProcess the process from which the data will be received
     * @param tag the unique identifier code
     * @param status pointer to the MPI status
     *
     * @return A pointer to the object returned.
     */
    boost::shared_ptr<CLASS> RecvObject(unsigned sourceProcess, unsigned tag, MPI_Status& status);

    /**
     * Post an asynchronous receive for an object
     *
     * @param sourceProcess the process from which the data will be received
     * @param tag the unique identifier code
     */
    void IRecvObject(unsigned sourceProcess, unsigned tag);

    /**
     * Obtain a proper object once a call to IRecv has completed
     *
     * @return a boost shared pointer to a receive object
     */
    boost::shared_ptr<CLASS> GetRecvObject();

     /**
     * Send and receive an object
     *
     * @param pSendObject a pointer to the object to send
     * @param destinationProcess the rank of the target process
     * @param sendTag the tag to send with.
     * @param sourceProcess the process from which the data will be received
     * @param sourceTag the tag to receive with
     * @param status a reference to an MPI_Status object.
     *
     * @return A pointer to the object returned.
     */
    boost::shared_ptr<CLASS> SendRecvObject(boost::shared_ptr<CLASS> const pSendObject, unsigned destinationProcess, unsigned sendTag, unsigned sourceProcess, unsigned sourceTag, MPI_Status& status);
};

#include <sstream>
#include <string>
#include <cstring>

// Implementation needs to be here, as CLASS could be anything
template<typename CLASS>
ObjectCommunicator<CLASS>::ObjectCommunicator()
    : mIsWriting(false)
{
    mSendBuffer.resize(PetscTools::GetNumProcs());
    mSendString.resize(PetscTools::GetNumProcs());
}

template<typename CLASS>
void ObjectCommunicator<CLASS>::SendObject(boost::shared_ptr<CLASS> const pObject, unsigned destinationProcess, unsigned tag)
{
    // Create an output archive
    std::ostringstream ss(std::ios::binary);
    boost::archive::binary_oarchive output_arch(ss);

    output_arch << pObject;

    const std::string send_msg = ss.str();

    // Get + send string length
    unsigned string_length = send_msg.size();
    MPI_Send(&string_length, 1, MPI_UNSIGNED, destinationProcess, tag, PetscTools::GetWorld());

    // Send archive data
    // The buffer is treated as const, but not specified as such by MPI_Send's signature
    char* send_buf = const_cast<char*>(send_msg.data());
    MPI_Send(send_buf, string_length, MPI_BYTE, destinationProcess, tag, PetscTools::GetWorld());
}

template<typename CLASS>
void ObjectCommunicator<CLASS>::ISendObject(boost::shared_ptr<CLASS> const pObject, unsigned destinationProcess, unsigned tag)
{
    MPI_Request request;

     // Create an output archive
    std::ostringstream ss(std::ios::binary);
    boost::archive::binary_oarchive output_arch(ss);

    output_arch << pObject;

    mSendString[destinationProcess] = ss.str();
    mSendBufferLength = mSendString[destinationProcess].size();

    // Make sure we are not going to overrun the asynchronous buffer size.
    assert(mSendBufferLength < MAX_BUFFER_SIZE);

    // Send archive data
    // The buffer is treated as const, but not specified as such by MPI_Send's signature
    mSendBuffer[destinationProcess] = const_cast<char*>(mSendString[destinationProcess].data());
    MPI_Isend(mSendBuffer[destinationProcess], mSendBufferLength, MPI_BYTE, destinationProcess, tag, PetscTools::GetWorld(), &request);
    MPI_Request_free(&request); //This is evil because it allows for another non-blocking send to overwrite the buffer
}

template<typename CLASS>
boost::shared_ptr<CLASS> ObjectCommunicator<CLASS>::RecvObject(unsigned sourceProcess, unsigned tag, MPI_Status& status)
{
    unsigned string_length = 0;
    MPI_Recv(&string_length, 1, MPI_UNSIGNED, sourceProcess, tag, PetscTools::GetWorld(), &status);

    char* recv_array = new char[string_length];
    MPI_Recv(recv_array, string_length, MPI_BYTE, sourceProcess , tag, PetscTools::GetWorld(), &status);

    // Extract a proper object from the buffer
    std::string recv_string(recv_array, string_length);
    delete[] recv_array;
    std::istringstream ss(recv_string, std::ios::binary);

    boost::shared_ptr<CLASS> p_recv_object(new CLASS);
    boost::archive::binary_iarchive input_arch(ss);

    input_arch >> p_recv_object;

    return p_recv_object;
}

template<typename CLASS>
void ObjectCommunicator<CLASS>::IRecvObject(unsigned sourceProcess, unsigned tag)
{
    assert(!mIsWriting);    // Make sure the buffers are not already being used by a previous call to IRecvObject.

    mIsWriting = true;

    mRecvBuffer = new char[MAX_BUFFER_SIZE];
    MPI_Irecv(mRecvBuffer, MAX_BUFFER_SIZE, MPI_BYTE, sourceProcess, tag, PetscTools::GetWorld(), &mMpiRequest);
}

template<typename CLASS>
boost::shared_ptr<CLASS> ObjectCommunicator<CLASS>::GetRecvObject()
{
    if (!mIsWriting)
    {
        EXCEPTION("No object to receive in ObjectCommunicator::GetRecvObject");
    }

    MPI_Status return_status;

    MPI_Wait(&mMpiRequest, &return_status);

    int recv_size;
    MPI_Get_count(&return_status, MPI_BYTE, &recv_size);

    // Extract a proper object from the buffer
    std::string recv_string(mRecvBuffer, recv_size);
    std::istringstream ss(recv_string, std::ios::binary);

    boost::shared_ptr<CLASS> p_recv_object(new CLASS);
    boost::archive::binary_iarchive input_arch(ss);

    input_arch >> p_recv_object;

    // Tidy up
    delete[] mRecvBuffer;

    mIsWriting = false;

    return p_recv_object;
}

template<typename CLASS>
boost::shared_ptr<CLASS> ObjectCommunicator<CLASS>::SendRecvObject(boost::shared_ptr<CLASS> const pSendObject, unsigned destinationProcess, unsigned sendTag, unsigned sourceProcess, unsigned sourceTag, MPI_Status& status)
{
    // Create an output archive
    std::ostringstream oss(std::ios::binary);
    boost::archive::binary_oarchive output_arch(oss);

    output_arch << pSendObject;

    std::string send_msg = oss.str();

    // Get + send string length
    unsigned send_string_length = send_msg.size();
    unsigned recv_string_length;

    MPI_Sendrecv(&send_string_length, 1, MPI_UNSIGNED, destinationProcess, sendTag, &recv_string_length, 1, MPI_UNSIGNED, sourceProcess, sourceTag, PetscTools::GetWorld(), &status);

    boost::scoped_array<char> recv_array(new char[recv_string_length]);

    // Send archive data
    char* send_buf = const_cast<char*>(send_msg.data());
    MPI_Sendrecv(send_buf, send_string_length, MPI_BYTE, destinationProcess, sendTag, recv_array.get(), recv_string_length, MPI_BYTE, sourceProcess, sourceTag, PetscTools::GetWorld(), &status);

    // Extract received object
    std::string recv_string(recv_array.get(), recv_string_length);
    std::istringstream iss(recv_string, std::ios::binary);

    boost::shared_ptr<CLASS> p_recv_object(new CLASS);
    boost::archive::binary_iarchive input_arch(iss);

    input_arch >> p_recv_object;

    return p_recv_object;
}

#endif // _OBJECTCOMMUNICATOR_HPP_
