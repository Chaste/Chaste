/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef _CVODECONTEXT_HPP_
#define _CVODECONTEXT_HPP_

#ifdef CHASTE_CVODE
#if CHASTE_SUNDIALS_VERSION >= 60000
#include <memory>

#include <nvector/nvector_serial.h>

/**
 * A special singleton class to manage the SUNContext object required
 * by Sundials 6.0 and above.
 *
 * This class is a singleton and an instance should be retrieved with:
 * std::shared_ptr<CvodeContextManager> ctx = CvodeContextManager::Instance();
 *
 * There is one instance per thread, since a SUNContext may not be shared between threads.
 * Ownership is shared: Sundials requires that every object created against a context is
 * freed before the context itself, but the context is created lazily on first use and so
 * would otherwise be torn down before the longer-lived objects that refer to it. Any class
 * holding Sundials resources should therefore keep a copy of the pointer returned here for
 * as long as it holds those resources, which keeps the context alive until the last user
 * has gone.
 */
class CvodeContextManager
{
private:
    /** The SUNContext **/
    SUNContext mSundialsContext;

    /**
     * Private constructor.
     * Use Instance() to access the context manager.
     */
    CvodeContextManager();

    CvodeContextManager(const CvodeContextManager&) = delete; // disable copy constructor
    CvodeContextManager& operator=(const CvodeContextManager&) = delete; // disable copy assignment

public:
    /**
     * @return a reference to the managed context object.
     *
     */
    SUNContext& GetSundialsContext();

    /**
     * @return a shared pointer to the context manager object for the calling thread.
     * The object is created the first time this method is called on that thread, and is
     * destroyed once the thread has exited and every holder of the returned pointer has
     * released it.
     */
    static std::shared_ptr<CvodeContextManager> Instance();

    /**
     * Public destructor is only called when the last holder of the shared pointer returned
     * by Instance() releases it.
     */
    ~CvodeContextManager();
};

#endif // CHASTE_SUNDIALS_VERSION >= 60000
#endif // CHASTE_CVODE

#endif // _CVODECONTEXT_HPP_
