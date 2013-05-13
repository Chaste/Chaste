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

#ifndef _OBJECTCOMMUNICATOR_HPP_
#define _OBJECTCOMMUNICATOR_HPP_

// Serialisation headers - must come first
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/shared_ptr.hpp>
#include "PetscTools.hpp" // For MPI methods
#include "Exception.hpp"

/**
 * This is a helper method to enable classes that can be serialized to be sent using
 * PetSc MPI communication. The object is serialized in to a string of characters, and then
 * de-serialized on the receive process.
 */
template<typename CLASS>
class ObjectCommunicator
{
private:

	/** A buffer for use in asyncronous communication */
	char* mRecvString;

	/** The size of a string we are waiting for in an asyncronous receive */
	unsigned mStringLength;

	/** An MPI_Request used in MPI_Irecv */
	MPI_Request mMpiRequest;

	/** A flag, used as a lock to ensure that we don't accidentally start overwriting the above buffer, mRecvString */
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
     * @param status a refernce to an MPI_Status object.
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
    MPI_Send(&string_length, 1, MPI_UNSIGNED, destinationProcess, tag, PETSC_COMM_WORLD);

    // Send archive data
    // The buffer is treated as const, but not specified as such by MPI_Send's signature
    char* send_buf = const_cast<char*>(send_msg.data());
    MPI_Send(send_buf, string_length, MPI_BYTE, destinationProcess, tag, PETSC_COMM_WORLD);
}

template<typename CLASS>
void ObjectCommunicator<CLASS>::ISendObject(boost::shared_ptr<CLASS> const pObject, unsigned destinationProcess, unsigned tag)
{
	MPI_Request request;

	 // Create an output archive
	std::ostringstream ss(std::ios::binary);
	boost::archive::binary_oarchive output_arch(ss);

	output_arch << pObject;

	const std::string send_msg = ss.str();

	// Get + send string length
	unsigned string_length = send_msg.size();
	MPI_Isend(&string_length, 1, MPI_UNSIGNED, destinationProcess, tag, PETSC_COMM_WORLD, &request);

	// Send archive data
	// The buffer is treated as const, but not specified as such by MPI_Send's signature
	char* send_buf = const_cast<char*>(send_msg.data());
	MPI_Isend(send_buf, string_length, MPI_BYTE, destinationProcess, tag, PETSC_COMM_WORLD, &request);
}

template<typename CLASS>
boost::shared_ptr<CLASS> ObjectCommunicator<CLASS>::RecvObject(unsigned sourceProcess, unsigned tag, MPI_Status& status)
{
    unsigned string_length = 0;
    MPI_Recv(&string_length, 1, MPI_UNSIGNED, sourceProcess, tag, PETSC_COMM_WORLD, &status);

    char* recv_string = new char[string_length];
    MPI_Recv(recv_string, string_length, MPI_BYTE, sourceProcess , tag, PETSC_COMM_WORLD, &status);

    // Extract a proper object from the buffer
    std::istringstream ss(std::ios::binary);
    ss.rdbuf()->pubsetbuf(recv_string, string_length);

    boost::shared_ptr<CLASS> p_recv_object(new CLASS);
    boost::archive::binary_iarchive input_arch(ss);

    input_arch >> p_recv_object;

    // Tidy up
    delete[] recv_string;

    return p_recv_object;
}

template<typename CLASS>
void ObjectCommunicator<CLASS>::IRecvObject(unsigned sourceProcess, unsigned tag)
{
	assert(!mIsWriting);	// Make sure the buffers are not already being used by a previous call to IRecvObject.

	mIsWriting = true;

	MPI_Request size_request;
    MPI_Irecv(&mStringLength, 1, MPI_UNSIGNED, sourceProcess, tag, PETSC_COMM_WORLD, &size_request);
    MPI_Wait(&size_request, MPI_STATUS_IGNORE);

    mRecvString = new char[mStringLength];
    MPI_Irecv(mRecvString, mStringLength, MPI_BYTE, sourceProcess , tag, PETSC_COMM_WORLD, &mMpiRequest);
}

template<typename CLASS>
boost::shared_ptr<CLASS> ObjectCommunicator<CLASS>::GetRecvObject()
{
	if (!mIsWriting)
	{
		EXCEPTION("No object to receive in ObjectCommunicator::GetRecvObject");
	}

    MPI_Wait(&mMpiRequest, MPI_STATUS_IGNORE);

    // Extract a proper object from the buffer
    std::istringstream ss(std::ios::binary);
    ss.rdbuf()->pubsetbuf(mRecvString, mStringLength);

    boost::shared_ptr<CLASS> p_recv_object(new CLASS);
    boost::archive::binary_iarchive input_arch(ss);

    input_arch >> p_recv_object;

    // Tidy up
    delete[] mRecvString;

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

    MPI_Sendrecv(&send_string_length, 1, MPI_UNSIGNED, destinationProcess, sendTag, &recv_string_length, 1, MPI_UNSIGNED, sourceProcess, sourceTag, PETSC_COMM_WORLD, &status);

    char recv_string[recv_string_length];

    // Send archive data
    char* send_buf = const_cast<char*>(send_msg.data());
    MPI_Sendrecv(send_buf, send_string_length, MPI_BYTE, destinationProcess, sendTag, recv_string, recv_string_length, MPI_BYTE, sourceProcess, sourceTag, PETSC_COMM_WORLD, &status);

    // Extract received object
    std::istringstream iss(std::ios::binary);
    iss.rdbuf()->pubsetbuf(recv_string, recv_string_length);

    boost::shared_ptr<CLASS> p_recv_object(new CLASS);
    boost::archive::binary_iarchive input_arch(iss);

    input_arch >> p_recv_object;

    return p_recv_object;
}

#endif // _OBJECTCOMMUNICATOR_HPP_
