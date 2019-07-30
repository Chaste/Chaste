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

#ifndef CELLID_HPP_
#define CELLID_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"
#include "PetscTools.hpp"

class CellId;

/**
 * Cell id class.
 *
 * Each Cell owns a CellPropertyCollection, which may include a shared pointer
 * to an object of this type. When a Cell divides a new object is created with
 * the new cell id.
 *
 * The CellId object that stores the value of a the cell identifier.
 */
class CellId : public AbstractCellProperty
{
private:

    /**
     * Cell Id
     */
    unsigned mCellId;

    /** maximum cell identifier. */
    static unsigned mMaxCellId;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mCellId;
        if (!PetscTools::IsParallel()) // This is to avoid changing the static i.d. in parallel simulations
        {
            archive & mMaxCellId;
        }
    }

public:

    /**
     * Constructor.
     *
     * This doesn't do anything and AssignCellId must be called before doing anything
     */
    CellId();

    /**
     * Destructor.
     */
    virtual ~CellId();

    /**
     * This assigns the cell id to be the maximum current cell id.
     * It then increments mMaxCellId.
     */
    void AssignCellId();

    /**
     * @return the maximum value of the cell identifier
     */
    unsigned GetMaxCellId() const;

    /**
     * @return #mCellId.
     */
    unsigned GetCellId() const;

    /**
     * Reset the maximum cell id to zero.
     */
    static void ResetMaxCellId();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellId)

#endif /* CELLID_HPP_ */
