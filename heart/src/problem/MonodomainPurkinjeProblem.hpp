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

#ifndef MONODOMAINPURKINJEPROBLEM_HPP_
#define MONODOMAINPURKINJEPROBLEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCardiacProblem.hpp"
#include "AbstractCardiacTissue.hpp"
#include "AbstractPurkinjeCellFactory.hpp"
#include "MonodomainTissue.hpp"


/**
 * Class which specifies and solves a monodomain problem with Purkinje fibres.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class MonodomainPurkinjeProblem : public AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 2>
{

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     * @param archive
     * @param version
     * \todo Serialization of Purkinje problems is untested
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)  // LCOV_EXCL_LINE
    {
        NEVER_REACHED; // If you remove this NEVER_REACHED, then:
        // please remove the //LCOV exclude above and the ones around the default constructor in the .cpp too!


//        archive & mPurkinjeVoltageColumnId;
//        archive & boost::serialization::base_object<AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, 2> >(*this);
    }

protected:
    /** Used by the writer. */
    unsigned mPurkinjeVoltageColumnId;

    /** @return newly created our tissue object. */
    AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM>* CreateCardiacTissue();

    /** @return newly created suitable solver for monodomain problems with Purkinje. */
    AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 2>* CreateSolver();

    /** @return newly created suitable (MixedDimensionMesh) mesh for monodomain problems with Purkinje. */
    virtual void CreateMeshFromHeartConfig();

    /**
     *  Overridden method which creates initial condition using Purkinje initial voltages
     *  as well.
     *  @return newly created initial condition
     */
    Vec CreateInitialCondition();

public:
    /**
     * Constructor
     * @param pCellFactory  user defined cell factory which shows how the tissue should create cells.
     */
    MonodomainPurkinjeProblem(AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory);

    /**
     * Constructor just used for archiving
     */
    MonodomainPurkinjeProblem();

    /**
     * Destructor
     */
    virtual ~MonodomainPurkinjeProblem();

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
     * Adds the Purkinje transmembrane potential to the results.
     *
     * @param time  the current time
     * @param voltageVec  the solution vector to write
     */
    virtual void WriteOneStep(double time, Vec voltageVec);

};


#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS2(MonodomainPurkinjeProblem, 2, 2)
EXPORT_TEMPLATE_CLASS2(MonodomainPurkinjeProblem, 3, 3)



#endif // MONODOMAINPURKINJEPROBLEM_HPP_
