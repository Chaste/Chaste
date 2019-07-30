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


#ifndef MONODOMAINPROBLEM_HPP_
#define MONODOMAINPROBLEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCardiacProblem.hpp"
#include "AbstractCardiacTissue.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "MonodomainTissue.hpp"


/**
 * Class which specifies and solves a monodomain problem.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class MonodomainProblem : public AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 1>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 1> >(*this);
        archive & mpMonodomainTissue;
    }

protected:
    /** The monodomain tissue object. */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;

public:
    /** @return Created monodomain tissue. */
    AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM>* CreateCardiacTissue();

    /** @return Created suitable solver for monodomain problems. */
    AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 1>* CreateSolver();

    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     *   create cells.
     */
    MonodomainProblem(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory);

    /**
     * Constructor just used for archiving
     */
    MonodomainProblem();

    /**
     * Destructor
     */
    virtual ~MonodomainProblem();

    /** @return the monodomain PDE */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM> * GetMonodomainTissue();

    /**
     *  Print out time and max/min voltage values at current time.
     *
     * @param time  the current time
     */
    void WriteInfo(double time);

    /**
     * Define what variables are written to the primary results file.
     * @param extending  whether we are extending an existing results file
     */
    virtual void DefineWriterColumns(bool extending);

    /**
     * Write one timestep of output data to the primary results file.
     * Adds the voltage to the results.
     *
     * @param time  the current time
     * @param voltageVec  the solution vector to write
     */
    virtual void WriteOneStep(double time, Vec voltageVec);
};

#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 1, 1)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 1, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 1, 3)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 2, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainProblem, 3, 3)

#endif /*MONODOMAINPROBLEM_HPP_*/
