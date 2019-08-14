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

#ifndef SRNCELLMODEL_HPP_
#define SRNCELLMODEL_HPP_

#include <vector>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "AbstractSrnModel.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "SimulationTime.hpp"

typedef boost::shared_ptr<AbstractSrnModel> AbstractSrnModelPtr;

/**
 * SRN model at the cell level, has representation for edges internally.
 * Also contains cell interior (cytoplasmic) SRN
 */
class SrnCellModel : public AbstractSrnModel
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
        archive & mEdgeSrnModels;
        archive & mInteriorSrnModel;
//        archive & boost::serialization::base_object<CellCycleModelOdeHandler>(*this);
//        archive & mInitialConditions;
//        archive & mStateSize;
    }

    std::vector<boost::shared_ptr<AbstractSrnModel>> mEdgeSrnModels;
    using abstractsrnmodel_t = std::vector<AbstractSrnModelPtr>;

    boost::shared_ptr<AbstractSrnModel> mInteriorSrnModel;
protected:

    /**
     * Copy constructor
     * @param rModel
     */
    SrnCellModel(const SrnCellModel &rModel);

public:

    /* Makes the class iterable which returns the individual edge SRN models */
    using iterator = abstractsrnmodel_t::iterator;
    using const_iterator = abstractsrnmodel_t::const_iterator;
    iterator begin() { return mEdgeSrnModels.begin(); }
    iterator end() { return mEdgeSrnModels.end(); }
    const_iterator begin() const { return mEdgeSrnModels.begin(); }
    const_iterator end() const { return mEdgeSrnModels.end(); }
    const_iterator cbegin() const { return mEdgeSrnModels.cbegin(); }
    const_iterator cend() const { return mEdgeSrnModels.cend(); }

    /**
     * Default constuctor.
     */
    SrnCellModel();

    /**
     * Destructor.
     */
    ~SrnCellModel();


    virtual void Initialise() override;

    virtual void SimulateToCurrentTime() override;

    virtual AbstractSrnModel* CreateSrnModel() override;

    void AddEdgeSrn(std::vector<AbstractSrnModelPtr> edgeSrn);

    void AddEdgeSrnModel(AbstractSrnModelPtr edgeSrn);

    void InsertEdgeSrn(unsigned index, AbstractSrnModelPtr edgeSrn);

    AbstractSrnModelPtr RemoveEdgeSrn(unsigned index);

    unsigned GetNumEdgeSrn();

    AbstractSrnModelPtr GetEdgeSrn(unsigned index);

    const std::vector<AbstractSrnModelPtr>& GetEdges();

    void SetInteriorSrnModel(AbstractSrnModelPtr interiorSrn);

    AbstractSrnModelPtr GetInteriorSrn();

    virtual void SetCell(CellPtr pCell) override;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SrnCellModel)

#endif /* SRNCELLMODEL_HPP_ */
