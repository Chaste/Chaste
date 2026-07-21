/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef ABSTRACTCELLPROPERTY_HPP_
#define ABSTRACTCELLPROPERTY_HPP_

#include <boost/shared_ptr.hpp>
#include <type_traits>
#include "ChasteSerialization.hpp"
#include "Identifiable.hpp"

/**
 * Generates a value that is unique per template instantiation and stable for the life of
 * the process, without using RTTI. Comparing two of these is a plain pointer comparison,
 * unlike comparing std::type_info objects, which can fall back to a string comparison when
 * the types being compared live in different shared libraries.
 */
class CellPropertyTypeId
{
public:
    /**
     * @return an ID unique to CLASS, ignoring any top-level cv-qualification on CLASS (so
     * Get<CellLabel>() and Get<const CellLabel>() return the same value, matching the
     * behaviour of typeid(), which this class replaces).
     */
    template<typename CLASS>
    static const void* Get()
    {
        return GetUnqualified<typename std::remove_cv<CLASS>::type>();
    }

private:

    /**
     * @return an ID unique to CLASS. Only ever instantiated with a non-cv-qualified CLASS,
     * so that differently cv-qualified spellings of the same type share one marker.
     */
    template<typename CLASS>
    static const void* GetUnqualified()
    {
        static char s_marker;
        return &s_marker;
    }
};

/**
 * Base class for cell properties.
 *
 * Each cell has a collection of cell properties, which may express such concepts as
 * mutation states, inherited labels, or per-cell data.
 */
class AbstractCellProperty : public Identifiable
{
private:

    /**
     * The number of cells with this cell property.
     */
    unsigned mCellCount;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * We don't have any member variables yet, but let's make things easy for the future.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCellCount;
    }

public:
    /**
     * Default constructor.
     */
    AbstractCellProperty();

    /**
     * Virtual destructor, to make this class polymorphic.
     */
    virtual ~AbstractCellProperty();

    /**
     * @return whether this property is a particular class.  This tests for exact
     * run-time class identity, and doesn't match subclasses.
     *
     * Subclasses that derive from the CellProperty<CONCRETE> mixin (rather than from
     * AbstractCellProperty directly) get this comparison via a cheap pointer comparison.
     * Subclasses that don't fall back to the slower typeid-based comparison, which requires
     * a default constructor and is comparatively expensive across shared library boundaries.
     *
     * It should be called like:
     *   bool healthy = p_property->IsType<HealthyMutationState>();
     */
    template<class CLASS>
    bool IsType() const
    {
        const void* p_this_type_id = GetTypeId();
        if (p_this_type_id != nullptr)
        {
            return p_this_type_id == CellPropertyTypeId::Get<CLASS>();
        }

        CLASS ref_obj;
        return IsSame(&ref_obj);
    }

    /**
     * @return whether this property is an instance of a particular class, or any subclass
     * of that class.
     *
     * It should be called like:
     *   bool is_mutation = p_property->IsSubType<AbstractCellMutationState>();
     */
    template<class BASECLASS>
    bool IsSubType() const
    {
        // We put a const_cast in here so users don't have to worry about whether the
        // property object they are testing is const or not.
        BASECLASS* p_subclass = dynamic_cast<BASECLASS*>(const_cast<AbstractCellProperty*>(this));
        return (p_subclass != nullptr);
    }

    /**
     * @return whether this property is the same class as another.
     * @param pOther  the property to compare against.
     */
    bool IsSame(const AbstractCellProperty* pOther) const;

    /**
     * @return whether this property is the same as another.
     * @param pOther  the property to compare against.
     */
    bool IsSame(boost::shared_ptr<const AbstractCellProperty> pOther) const;

    /**
     * Increment #mCellCount.
     */
    void IncrementCellCount();

    /**
     * Decrement #mCellCount.
     */
    void DecrementCellCount();

    /**
     * @return #mCellCount
     */
    unsigned GetCellCount() const;

protected:

    /**
     * @return an ID unique to this object's concrete class, for use by IsType(), or nullptr
     * if the concrete class hasn't opted into the fast comparison path (see CellProperty).
     */
    virtual const void* GetTypeId() const
    {
        return nullptr;
    }
};

/**
 * Mixin that gives a concrete cell property class a fast, RTTI-free identity for IsType() to
 * compare against, by overriding GetTypeId(). Concrete property classes should derive from
 * this instead of from AbstractCellProperty (or another abstract property base) directly, e.g.
 *   class CellLabel : public CellProperty<CellLabel> { ... };
 *   class WildTypeCellMutationState : public CellProperty<WildTypeCellMutationState, AbstractCellMutationState> { ... };
 *
 * @tparam CONCRETE the concrete (most-derived) class inheriting this mixin.
 * @tparam BASE the class CONCRETE would otherwise have derived from directly.
 */
template<typename CONCRETE, typename BASE = AbstractCellProperty>
class CellProperty : public BASE
{
public:

    // Inherit BASE's constructors
    using BASE::BASE;

protected:

    const void* GetTypeId() const override
    {
        return CellPropertyTypeId::Get<CONCRETE>();
    }
};

#endif /* ABSTRACTCELLPROPERTY_HPP_ */
