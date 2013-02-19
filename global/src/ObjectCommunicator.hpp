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

#include <string>
#include <assert.h>
#include <sstream>
#include "PetscTools.hpp"

/**
 * This is a helper method to enable classes that can be serialized to be sent using
 * PetSc MPI communication. The object is serialized in to a string of characters, and then
 * de-serialized on the receive process.
 */
template<typename CLASS>
class ObjectCommunicator
{

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
    void SendObject(CLASS* const pObject, unsigned destinationProcess, unsigned tag);

    /**
     * Receive an object
     *
     * @param sourceProcess the process from which the data will be received
     * @param tag the unique identifier code
     * @param status pointer to the MPI status
     *
     * @return A pointer to the object returned.
     */
    CLASS* RecvObject(unsigned sourceProcess, unsigned tag, MPI_Status& status);
 
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
    CLASS* SendRecvObject(CLASS* const pSendObject, unsigned destinationProcess, unsigned sendTag, unsigned sourceProcess, unsigned sourceTag, MPI_Status& status);   
};

#endif // _OBJECTCOMMUNICATOR_HPP_
