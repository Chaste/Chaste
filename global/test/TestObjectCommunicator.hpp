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

#ifndef _TESTOBJECTCOMMUNICATOR_HPP_
#define _TESTOBJECTCOMMUNICATOR_HPP_

#include <cxxtest/TestSuite.h>

// Needs to be first as includes archive headers
#include "ObjectCommunicator.hpp"

#include <vector>

#include "PetscTools.hpp"
#include "ClassOfSimpleVariables.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestObjectCommunicator: public CxxTest::TestSuite
{

public:

    void TestSendingClass()
    {
        MPI_Status status;
        ObjectCommunicator<ClassOfSimpleVariables> communicator;

        if (PetscTools::AmMaster())
        {
            // Create a simple class to send
            std::vector<double> doubles(3);
            doubles[0] = 1.1;
            doubles[1] = 1.2;
            doubles[2] = 1.3;

            std::vector<bool> bools(2);
            bools[0] = true;
            bools[1] = true;

            boost::shared_ptr<ClassOfSimpleVariables> p_new_class(new ClassOfSimpleVariables(42,"hello",doubles,bools));

            // Send the class
            for (unsigned p=1; p < PetscTools::GetNumProcs(); p++)
            {
                // Arguments are object, destination, tag
                communicator.SendObject(p_new_class, p, 123);
            }

        }
        else
        {
            boost::shared_ptr<ClassOfSimpleVariables> p_recv_class;

            p_recv_class = communicator.RecvObject(0, 123, status);

            // Check that the values are correct
            TS_ASSERT_EQUALS(p_recv_class->GetNumber(),42);
            TS_ASSERT_EQUALS(p_recv_class->GetString(),"hello");
            TS_ASSERT_EQUALS(p_recv_class->GetVectorOfDoubles().size(),3u);
            TS_ASSERT_EQUALS(p_recv_class->GetVectorOfBools().size(),2u);

            TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[0],1.1,1e-12);
            TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[1],1.2,1e-12);
            TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[2],1.3,1e-12);

            TS_ASSERT(p_recv_class->GetVectorOfBools()[0]);
            TS_ASSERT(p_recv_class->GetVectorOfBools()[1]);

        }
    }

    /* This test *may* deadlock if we don't use non-blocking communication */
    void TestNonBlockingSendingClass()
    {
        ObjectCommunicator<ClassOfSimpleVariables> communicator;

        {
            // Create a simple class to send
            std::vector<double> doubles(3);
            doubles[0] = 1.1;
            doubles[1] = 1.2;
            doubles[2] = 1.3 + PetscTools::GetMyRank();

            std::vector<bool> bools(2);
            bools[0] = true;
            bools[1] = true;

            boost::shared_ptr<ClassOfSimpleVariables> p_new_class(new ClassOfSimpleVariables(42,"hello",doubles,bools));

            // Send the class
            for (unsigned p=0; p < PetscTools::GetNumProcs(); p++)
            {
                if (p != PetscTools::GetMyRank())
                {
                    p_new_class->mVectorOfDoubles[1] = 1.2 + p;
                    // Arguments are object, destination, tag
                    communicator.ISendObject(p_new_class, p, 123 + PetscTools::GetMyRank());
                }
            }
        }

        {
            boost::shared_ptr<ClassOfSimpleVariables> p_recv_class;

            for (unsigned p=0; p < PetscTools::GetNumProcs(); p++)
            {
                if (p != PetscTools::GetMyRank())
                {
                    communicator.IRecvObject(p, 123 + p);
                    p_recv_class = communicator.GetRecvObject();

                    // Check that the values are correct
                    TS_ASSERT_EQUALS(p_recv_class->GetNumber(),42);
                    TS_ASSERT_EQUALS(p_recv_class->GetString(),"hello");
                    TS_ASSERT_EQUALS(p_recv_class->GetVectorOfDoubles().size(),3u);
                    TS_ASSERT_EQUALS(p_recv_class->GetVectorOfBools().size(),2u);

                    TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[0],1.1,1e-12);
                    TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[1],1.2 + PetscTools::GetMyRank(),1e-12);
                    TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[2],1.3 + p,1e-12);

                    TS_ASSERT(p_recv_class->GetVectorOfBools()[0]);
                    TS_ASSERT(p_recv_class->GetVectorOfBools()[1]);
                }
            }
        }
        PetscTools::Barrier("Make sure that no ISendObject buffers are in use before proceeding");
    }

    /** We cannot pre-post Irecv because we need to know the (dynamic) size of the object being sent first */
    void TestNonBlockingRecvClass()
    {
        ObjectCommunicator<ClassOfSimpleVariables> communicator;

        if (PetscTools::AmMaster())
        {
            // Create a simple class to send
            std::vector<double> doubles(3);
            doubles[0] = 1.1;
            doubles[1] = 1.2;
            doubles[2] = 1.3;

            std::vector<bool> bools(2);
            bools[0] = true;
            bools[1] = true;

            boost::shared_ptr<ClassOfSimpleVariables> p_new_class(new ClassOfSimpleVariables(42,"hello",doubles,bools));

            // Send the class
            for (unsigned p=1; p < PetscTools::GetNumProcs(); p++)
            {
                // Arguments are object, destination, tag
                p_new_class->mVectorOfDoubles[1] = 1.2 + p;
                communicator.ISendObject(p_new_class, p, 123);
            }
        }
        else
        {
            boost::shared_ptr<ClassOfSimpleVariables> p_recv_class;

            TS_ASSERT_THROWS_THIS(p_recv_class = communicator.GetRecvObject(), "No object to receive in ObjectCommunicator::GetRecvObject");

            // Receive the string. This returns before the receive is complete.
            communicator.IRecvObject(0, 123);

            // Get an object from the communicator. Implicit MPI_Wait.
            p_recv_class = communicator.GetRecvObject();

            // Check that the values are correct
            TS_ASSERT_EQUALS(p_recv_class->GetNumber(),42);
            TS_ASSERT_EQUALS(p_recv_class->GetString(),"hello");
            TS_ASSERT_EQUALS(p_recv_class->GetVectorOfDoubles().size(),3u);
            TS_ASSERT_EQUALS(p_recv_class->GetVectorOfBools().size(),2u);

            TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[0],1.1,1e-12);
            TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[1],1.2+PetscTools::GetMyRank(),1e-12);
            TS_ASSERT_DELTA(p_recv_class->GetVectorOfDoubles()[2],1.3,1e-12);

            TS_ASSERT(p_recv_class->GetVectorOfBools()[0]);
            TS_ASSERT(p_recv_class->GetVectorOfBools()[1]);
        }
        PetscTools::Barrier("Make sure that no ISendObject buffers are in use before proceeding");
    }

    void TestSendRecv()
    {
        if (PetscTools::GetNumProcs() == 2)
        {
            MPI_Status status;
            ObjectCommunicator<ClassOfSimpleVariables> communicator;

            // Create a simple class to send

            std::vector<double> doubles(3);
            doubles[0] = 1.1;
            doubles[1] = 1.2;
            doubles[2] = 1.3;

            std::vector<bool> bools(2);
            bools[0] = true;
            bools[1] = true;

            boost::shared_ptr<ClassOfSimpleVariables> p_send_class(new ClassOfSimpleVariables(42,"hello",doubles,bools));

            boost::shared_ptr<ClassOfSimpleVariables> p_class(new ClassOfSimpleVariables);

            // Arguments are object, destination, tag
            p_class = communicator.SendRecvObject(p_send_class, 1-PetscTools::GetMyRank(), 123, 1-PetscTools::GetMyRank(), 123, status);

            // Check that the values are correct
            TS_ASSERT_EQUALS(p_class->GetNumber(),42);
            TS_ASSERT_EQUALS(p_class->GetString(),"hello");
            TS_ASSERT_EQUALS(p_class->GetVectorOfDoubles().size(),3u);
            TS_ASSERT_EQUALS(p_class->GetVectorOfBools().size(),2u);

            TS_ASSERT_DELTA(p_class->GetVectorOfDoubles()[0],1.1,1e-12);
            TS_ASSERT_DELTA(p_class->GetVectorOfDoubles()[1],1.2,1e-12);
            TS_ASSERT_DELTA(p_class->GetVectorOfDoubles()[2],1.3,1e-12);

            TS_ASSERT(p_class->GetVectorOfBools()[0]);
            TS_ASSERT(p_class->GetVectorOfBools()[1]);
        }
    }
};

#endif //_TESTOBJECTCOMMUNICATOR_HPP_

