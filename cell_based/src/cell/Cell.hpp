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

#ifndef CELL_HPP_
#define CELL_HPP_

#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellMutationState.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "CellData.hpp"
#include "CellVecData.hpp"

#include "AbstractCellCycleModel.hpp"
#include "AbstractSrnModel.hpp"
#include "CellPropertyCollection.hpp"

class AbstractCellCycleModel; // Circular definition (cells need to know about cycle models and vice-versa).
class AbstractSrnModel; // Circular definition (cells need to know about subcellular reaction network models and vice-versa).
class Cell;

/** Cells shouldn't be copied - it doesn't make sense.  So all access is via this pointer type. */
typedef boost::shared_ptr<Cell> CellPtr;

/**
 * null_deleter means "doesn't delete" rather than "deletes nulls".
 *
 * Sometimes it is desirable to create a shared_ptr to an already existing object, so that the shared_ptr
 * does not attempt to destroy the object when there are no more references left. As an example, the
 * factory function:
 *
 * shared_ptr<X> createX();
 * in certain situations may need to return a pointer to a statically allocated X instance.
 *
 * The solution is to use a custom deleter that does nothing:
 */
struct null_deleter
{
    /** Does not delete */
    void operator()(void const *) const
    {
    }
};

/**
 * Cell is the basic container for all the biological information about a cell.
 * It contains the cell-cycle model and all other biological properties such as mutation
 * state, cell type, whether it is undergoing apoptosis or not.
 *
 * This class should not store any spatial information - cells are linked to space by the AbstractCellPopulation subclasses.
 */
class Cell : private boost::noncopyable, public boost::enable_shared_from_this<Cell>
{
private:

    /** Caches the result of ReadyToDivide() so Divide() can look at it. */
    bool mCanDivide;

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
        // These first four are also dealt with by {load,save}_construct_data
        archive & mCanDivide;
        archive & mpCellCycleModel;
        archive & mpSrnModel;
        archive & mUndergoingApoptosis;
        archive & mDeathTime;
        archive & mStartOfApoptosisTime;
        archive & mApoptosisTime;
        archive & mIsDead;
        archive & mIsLogged;
    }

protected:

    /** The cell's property collection. */
    CellPropertyCollection mCellPropertyCollection;

    /** The cell's cell-cycle model. */
    AbstractCellCycleModel* mpCellCycleModel;

    /** The cell's sub-cellular reaction network (SRN) model. */
    AbstractSrnModel* mpSrnModel;

    /** When the cell will/did die. */
    double mDeathTime;

    /** When the cell was commanded to start apoptosis. */
    double mStartOfApoptosisTime;

    /** The time it takes for a cell to fully undergo apoptosis. Has units of hours. */
    double mApoptosisTime;

    /** Whether the cell is currently in apoptosis - don't divide. */
    bool mUndergoingApoptosis;

    /**
     * Whether the cell is dead or not (they exist in the CellPopulation until they are
     * removed by AbstractCellPopulation::RemoveDeadCells().
     */
    bool mIsDead;

    /** Whether the cell is being tracked specially. */
    bool mIsLogged;

public:

    /**
     * Create a new cell.
     *
     * @param pMutationState the mutation state of the cell
     * @param pCellCycleModel  the cell-cycle model to use to decide when the cell divides.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     * @param pSrnModel  the sub-cellular reaction network model.
     *      This MUST be allocated using new, and will be deleted when the cell is destroyed.
     *      (Defaults to NULL and is replcaed by a new NullSrnModel in the constructr)
     * @param archiving  whether this constructor is being called by the archiver - do things slightly differently! (defaults to false)
     * @param cellPropertyCollection the cell property collection (defaults to NULL)
     */
    Cell(boost::shared_ptr<AbstractCellProperty> pMutationState,
         AbstractCellCycleModel* pCellCycleModel,
         AbstractSrnModel* pSrnModel=nullptr,
         bool archiving=false,
         CellPropertyCollection cellPropertyCollection=CellPropertyCollection());

    /**
     * Destructor, which frees the memory allocated for our cell-cycle model.
     */
    ~Cell();

    /**
     * @return the cell's proliferative type.
     */
    boost::shared_ptr<AbstractCellProliferativeType> GetCellProliferativeType() const;

    /**
     * Set the cell's proliferative type.
     *
     * @param pProliferativeType the cell's new proliferative type
     */
    void SetCellProliferativeType(boost::shared_ptr<AbstractCellProperty> pProliferativeType);

    /**
     * Set the birth time of the cell - can be negative so that your cells have an age when a simulation begins
     *
     * @param birthTime  The time the cell was born (in hours)
     */
    void SetBirthTime(double birthTime);

    /**
     * Change the cell-cycle model used. This takes effect immediately.
     *
     * @param pCellCycleModel pointer to the cell-cycle model to use
     */
    void SetCellCycleModel(AbstractCellCycleModel* pCellCycleModel);

    /**
     * @return a pointer to the Cell's cell-cycle model.
     */
    AbstractCellCycleModel* GetCellCycleModel() const;

    /**
     * Calls Initialise on the cell-cycle model associated with this cell.
     */
    void InitialiseCellCycleModel();

    /**
     * Change the SRN model used. This takes effect immediately.
     *
     * @param pSrnModel pointer to the SRN model to use
     */
    void SetSrnModel(AbstractSrnModel* pSrnModel);

    /**
     * @return a pointer to the Cell's SRN model.
     */
    AbstractSrnModel* GetSrnModel() const;

    /**
     * Calls Initialise on the SRN model associated with this cell.
     */
    void InitialiseSrnModel();

    /**
     * @return the cell's age from its cell-cycle model.
     */
    double GetAge() const;

    /**
     * @return the cell's birth time from its cell-cycle model.
     */
    double GetBirthTime() const;

    /**
     * @return the time at which apoptosis was commanded to start.
     */
    double GetStartOfApoptosisTime() const;

    /**
     * @return mApoptosisTime
     */
    double GetApoptosisTime() const;

    /**
     * Set mApoptosisTime.
     *
     * @param apoptosisTime the new value of mApoptosisTime
     */
    void SetApoptosisTime(double apoptosisTime);

    /**
     * @return the cell's current mutation state.
     */
    boost::shared_ptr<AbstractCellMutationState> GetMutationState() const;

    /**
     * Get the CellData associated with the cell.
     *
     * @return a pointer to the cell data
     */
    boost::shared_ptr<CellData> GetCellData() const;

    /**
     * Checks whether there is CellVecData associated with this cell
     *
     * @return whether the current cell stores CellVecData
     */
    bool HasCellVecData() const;

    /**
     * Get the CellVecData associated with the cell.
     *
     * @return a pointer to the cell Vec data
     */
    boost::shared_ptr<CellVecData> GetCellVecData() const;

    /**
     * Set the cell's mutation state.
     *
     * @param pMutationState the cell's new mutation state
     */
    void SetMutationState(boost::shared_ptr<AbstractCellProperty> pMutationState);

    /**
     * @return reference to #mCellPropertyCollection.
     */
    CellPropertyCollection& rGetCellPropertyCollection();

    /**
     * @return reference to #mCellPropertyCollection (used in archiving).
     */
    const CellPropertyCollection& rGetCellPropertyCollection() const;

    /**
     * Add a cell property to the cell. Use this method instead of calling
     *     rGetCellPropertyCollection().AddProperty()
     * directly, to ensure that the cell property keeps correct track of the
     * number of cells with it (if this is done).
     *
     * @param rProperty the property to add
     */
    void AddCellProperty(const boost::shared_ptr<AbstractCellProperty>& rProperty);

    /**
     * Remove a cell property of the given type. Use this method instead of
     * calling
     *     rGetCellPropertyCollection().AddProperty()
     * directly, to ensure that the cell property keeps correct track of the
     * number of cells with it (if this is done).
     */
    template<typename CLASS>
    void RemoveCellProperty()
    {
        bool cell_has_property = false;

        for (std::set<boost::shared_ptr<AbstractCellProperty> >::iterator property_iter = mCellPropertyCollection.Begin();
             property_iter != mCellPropertyCollection.End();
             ++property_iter)
        {
            if ((*property_iter)->IsType<CLASS>())
            {
                cell_has_property = true;
                (*property_iter)->DecrementCellCount();
                break;
            }
        }

        if (cell_has_property)
        {
            mCellPropertyCollection.RemoveProperty<CLASS>();
        }
    }

    /**
     * @return whether the cell property collection contains a property that has the exact type CLASS.
     * Just calls mCellPropertyCollection.HasProperty().
     */
    template<typename CLASS>
    bool HasCellProperty() const
    {
        return mCellPropertyCollection.HasProperty<CLASS>();
    }

    /**
     * @return whether this cell is ready to divide at the present simulation time.
     * MUST be called before calling Divide().
     */
    bool ReadyToDivide();

    /**
     * Divide this cell to produce a daughter cell.
     * ReadyToDivide MUST have been called at the current time, and returned true.
     *
     * @return the new daughter cell
     */
    CellPtr Divide();

    /**
     * Make the cell enter apoptosis and sets #mDeathTime using the apoptosis
     * time as defined by mApoptosisTime.
     *
     * @param setDeathTime whether we tell the cell exactly when to die (defaults to true)
     */
    void StartApoptosis(bool setDeathTime=true);

    /**
     * This labels the cell as dead, it does not delete the cell, it remains
     * in the CellPopulation until AbstractCellPopulation::RemoveDeadCells() is called.
     */
    void Kill();

    /**
     * @return whether the cell is undergoing apoptosis or not.
     */
    bool HasApoptosisBegun() const;

    /**
     * @return How long until the cell dies (if it is in apoptosis, throws an exception if not)
     */
    double GetTimeUntilDeath() const;

    /**
     * @return whether the cell is dead or undergoing apoptosis.
     */
    bool IsDead();

    /**
     * Sets a flag to perform special output on this cell only.
     */
    void SetLogged();

    /**
     * @return Whether the cell is being tracked.
     */
    bool IsLogged();

    /**
     * Give the Cell an index which it passes to its children.
     *
     * @param pCellAncestor the cell's ancestor
     */
    void SetAncestor(boost::shared_ptr<AbstractCellProperty> pCellAncestor);

    /**
     * @return The ancestor object, inherited from parents or set using the method above,
     * used for monoclonality experiments.
     */
    unsigned GetAncestor() const;

    /**
     * @return The cell identifier.
     */
    unsigned GetCellId() const;
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Cell)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Cell.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Cell * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const boost::shared_ptr<AbstractCellMutationState> p_mutation_state = t->GetMutationState();
    ar & p_mutation_state;

    const AbstractCellCycleModel* const p_cell_cycle_model = t->GetCellCycleModel();
    ar & p_cell_cycle_model;

    const AbstractSrnModel* const p_srn_model = t->GetSrnModel();
    ar & p_srn_model;

    const CellPropertyCollection& r_cell_property_collection = t->rGetCellPropertyCollection();
    ar & r_cell_property_collection;
}

/**
 * De-serialize constructor parameters and initialize a Cell.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Cell * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    ar & p_mutation_state;

    AbstractCellCycleModel* p_cell_cycle_model;
    ar & p_cell_cycle_model;

    AbstractSrnModel* p_srn_model;
    ar & p_srn_model;

    bool archiving = true;

    CellPropertyCollection cell_property_collection;
    ar & cell_property_collection;

    // Invoke inplace constructor to initialize instance
    ::new(t)Cell(p_mutation_state, p_cell_cycle_model, p_srn_model, archiving, cell_property_collection);
}
}
} // namespace ...

#endif /*CELL_HPP_*/
