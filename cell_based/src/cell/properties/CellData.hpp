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

#ifndef CELLDATA_HPP_
#define CELLDATA_HPP_

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <vector>

#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include "Exception.hpp"

/**
 * CellData class.
 *
 * This cell property allows each cell to store one or more 'named' doubles associated with it,
 * for example corresponding to the intracellular oxygen concentration. Other classes may interrogate
 * or modify the values stored in this class.
 *
 * Within the Cell constructor, an empty CellData object is created and passed to the Cell
 * (unless there is already a CellData object present in mCellPropertyCollection).
 */
class CellData : public AbstractCellProperty
{
private:

    /**
     * The cell data.
     */
    std::map<std::string, double> mCellData;

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
        archive & mCellData;
    }

public:

    /**
     * We need the empty virtual destructor in this class to ensure Boost
     * serialization works correctly with static libraries.
     */
    virtual ~CellData();

    /**
     * This assigns the cell data.
     *
     * @param rVariableName the name of the data to be set.
     * @param data the value to set it to.
     */
    void SetItem(const std::string& rVariableName, double data);

    /**
     * @return data.
     *
     * @param rVariableName the index of the data required.
     * throws if rVariableName has not been stored
     */
    double GetItem(const std::string& rVariableName) const;

    /**
     * @return number of data items
     */
    unsigned GetNumItems() const;

    /**
     * @return all keys.
     *
     * According to STL these are sorted in lexicographical/alphabetic order (so that the ordering here is predictable).
     */
    std::vector<std::string> GetKeys() const;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellData)

#endif /* CELLDATA_HPP_ */
