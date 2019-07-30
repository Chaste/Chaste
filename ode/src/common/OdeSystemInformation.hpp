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


#ifndef _ODESYSTEMINFORMATION_HPP_
#define _ODESYSTEMINFORMATION_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractOdeSystemInformation.hpp"

/**
 * A concrete implementation of AbstractOdeSystemInformation, that uses templates
 * to provide an implementation for any ODE system class.
 *
 * All ODE system developers need to do is provide a specialisation of the
 * Initialise method of this class, and set mpSystemInfo in their constructor:
 *   mpSystemInfo = OdeSystemInformation<CLASS>::Instance();
 *
 * This class contains all the machinery to make it a singleton, hence providing
 * exactly one instance per value of the template parameter.
 */
template <class ODE_SYSTEM>
class OdeSystemInformation : public AbstractOdeSystemInformation
{
private:

    /**
     * The single instance of this class, for this ODE_SYSTEM.
     *
     * \todo see if using weak_ptr would work and give funkier semantics
     *   (automatically destroy the singleton when no ODE systems were using it)
     */
    static boost::shared_ptr<OdeSystemInformation<ODE_SYSTEM> > mpInstance;

protected:

    /**
     * Default constructor.
     *
     * Not user accessible - to obtain an instance of this class use the Instance
     * method.
     */
    OdeSystemInformation();

    /**
     * Copy constructor.
     */
    OdeSystemInformation(const OdeSystemInformation<ODE_SYSTEM>&);

    /**
     * @return reference to this object (as language convention)
     * Overloaded assignment operator.
     */
    OdeSystemInformation& operator= (const OdeSystemInformation<ODE_SYSTEM>&);

    /**
     * Generic implementation of Initialise, which does nothing.
     *
     * Developers should specialise this method to their ODE system.  For example,

template<>
void OdeSystemInformation<MyNewOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("Units_1");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}
     */
    void Initialise();

public:

    /**
     * @return a pointer to the singleton instance, creating it if necessary.
     */
    static boost::shared_ptr<OdeSystemInformation<ODE_SYSTEM> > Instance();
};

template<class ODE_SYSTEM>
boost::shared_ptr<OdeSystemInformation<ODE_SYSTEM> > OdeSystemInformation<ODE_SYSTEM>::Instance()
{
    if (!mpInstance)
    {
        mpInstance.reset(new OdeSystemInformation<ODE_SYSTEM>);
        mpInstance->Initialise();
    }
    return mpInstance;
}

template<class ODE_SYSTEM>
OdeSystemInformation<ODE_SYSTEM>::OdeSystemInformation()
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == nullptr);
}

template<class ODE_SYSTEM>
void OdeSystemInformation<ODE_SYSTEM>::Initialise()
{
    // does nothing; designed to be specialised
}

/**
 * Definition of the instance static member.
 */
template<class ODE_SYSTEM>
boost::shared_ptr<OdeSystemInformation<ODE_SYSTEM> > OdeSystemInformation<ODE_SYSTEM>::mpInstance;


#endif /*_ODESYSTEMINFORMATION_HPP_*/
