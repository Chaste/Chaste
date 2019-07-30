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

#ifndef COMBINEDODESYSTEMINFORMATION_HPP_
#define COMBINEDODESYSTEMINFORMATION_HPP_

#include <vector>

#include <boost/shared_ptr.hpp>

#include "AbstractOdeSystemInformation.hpp"
#include "AbstractOdeSystem.hpp"

/**
 * Provide information about a CombinedOdeSystem, by combining the
 * information provided by the subsystems.
 */
class CombinedOdeSystemInformation : public AbstractOdeSystemInformation
{
private:

    /**
     * A convenience structure to keep track of 'singleton' instances.
     */
    struct InstancePointers
    {
        /** The 'singleton' instance. */
        boost::shared_ptr<CombinedOdeSystemInformation> pInfoInstance;
        /** The subsystem information objects that contribute to this 'singleton' instance. */
        std::vector<boost::shared_ptr<const AbstractOdeSystemInformation> > subsystemInformation;
    };

    /**
     * The single instance of this class for each vector of sub-systems.
     */
    static std::vector<struct InstancePointers> msInstances;

protected:

    /**
     * Main constructor.
     *
     * Not user accessible - to obtain an instance of this class use the Instance method.
     *
     * @param rSubsystemInfo  the information objects of the ODE systems used to construct
     *     the system we are providing information about.
     */
    CombinedOdeSystemInformation(const std::vector<boost::shared_ptr<const AbstractOdeSystemInformation> >& rSubsystemInfo);

    /**
     * Copy constructor.  Not defined.
     */
    CombinedOdeSystemInformation(const CombinedOdeSystemInformation&);

    /**
     * Overloaded assignment operator.  Not defined.
     * @return reference to this object (as language convention)
     *
     */
    CombinedOdeSystemInformation& operator= (const CombinedOdeSystemInformation&);

    /**
     * We need to provide an Initialise method, but in this case all the work
     * is done by our constructor, so this method does nothing.
     */
    void Initialise();

public:

    /**
     * @return a pointer to the singleton instance, creating it if necessary.
     *
     * @param rSubsystems  the ODE systems used to construct the system we are providing information about.
     */
    static boost::shared_ptr<CombinedOdeSystemInformation> Instance(const std::vector<AbstractOdeSystem*>& rSubsystems);
};

#endif /*COMBINEDODESYSTEMINFORMATION_HPP_*/
