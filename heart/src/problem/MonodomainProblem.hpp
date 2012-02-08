/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
    /** Create our monodomain tissue. */
    AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM>* CreateCardiacTissue();

    /** Create an suitable solver for monodomain problems. */
    AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 1>* CreateSolver();

public:
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
    ~MonodomainProblem();

    /** Get the monodomain PDE */
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
     * Adds the extracellular potential to the results.
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
