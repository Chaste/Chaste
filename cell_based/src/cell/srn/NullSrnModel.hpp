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

#ifndef NULLSRNMODEL_HPP_
#define NULLSRNMODEL_HPP_

#include <vector>
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractSrnModel.hpp"
#include "SimulationTime.hpp"

class AbstractSrnModel;

/**
 * This class contains a dummy/null sub-cellular reaction network (SRN) model that can be used for
 * any cell-cycle models that do not have an SRN, or where the SRN is combined with the main cell-cycle
 * model.
 */
class NullSrnModel : public AbstractSrnModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSrnModel>(*this);
    }

    /**
     * Protected copy-constructor for use by CreateSrnModel().  The only way for external code to create a copy of a SRN model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent SRN model will have had ResetForDivision() called just before CreateSrnModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel  the SRN model to copy.
     */
    NullSrnModel(const NullSrnModel& rModel);

public:

    /**
     * Default constructor.
     *
     * Creates an NullSrnModel and calls SetBirthTime() on the
     * AbstractSrnModel to make sure that can be set 'back in time' for
     * cells which did not divide at the current time.
     */
    NullSrnModel();

    /**
     * Destructor.
     */
    virtual ~NullSrnModel();

    /*
     * Overridden SimulateToCurrentTime() method.
     *
     * Here we don't do anything as this is a null SRN model.
     */
    void SimulateToCurrentTime();

    /*
     * Overridden CreateSrnModel() method.
     *
     * Builder method to create new instances of the SRN model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     *
     * This method is called by Cell::Divide() to create a SRN
     * model for the daughter cell.  Note that the parent SRN
     * model will have had ResetForDivision() called just before
     * CreateSrnModel() is called, so performing an exact copy of the
     * parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @return new SRN model
     */
    AbstractSrnModel* CreateSrnModel();
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(NullSrnModel)

#endif /* NULLSRNMODEL_HPP_ */
