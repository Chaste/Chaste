/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef SERIALIZABLESINGLETON_HPP_
#define SERIALIZABLESINGLETON_HPP_

#include <boost/utility.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>

/**
 * This is a "wrapper" allowing more straightforward serialization of singleton classes.
 * Any singleton class which needs to be serialized should inherit from this base.  It
 * provides both part of the "singleton-ness" (by inheriting from boost::noncopyable),
 * and also a method GetSerializationWrapper().  Users of the singleton which wish to
 * serialize it should not do so directly.  Instead, they should call GetSerializationWrapper
 * and serialize the returned pointer.  Doing so will ensure that only a single global
 * instance of the singleton is maintained when loading from an archive.
 *
 * Note that if this is not done, and the singleton is serialized directly via the
 * instance pointer, then objects loaded from the archive will refer to a different instance
 * of the singleton from other code!
 *
 * Usage examples:
 *
 * For saving:
 * \code
            SerializableSingleton<RandomNumberGenerator>* const p_wrapper = p_gen->GetSerializationWrapper();
            output_arch << p_wrapper;
 * \endcode
 *
 * For loading:
 * \code
            SerializableSingleton<RandomNumberGenerator>* p_wrapper;
            input_arch >> p_wrapper;
 * \endcode
 *
 * Within a serialize method:
 * \code
            SerializableSingleton<RandomNumberGenerator>* p_wrapper = p_gen->GetSerializationWrapper();
            archive & p_wrapper;
 * \endcode
 *
 * Note that immediately after a load the wrapper pointer loaded into becomes invalid; call
 * GetSerializationWrapper again if you need a new wrapper for a subsequent save.
 */
template<class SINGLETON_CLASS>
class SerializableSingleton : boost::noncopyable
{
public:
    /**
     * Get the wrapper object to use to serialize the related singleton instance.
     */
    SerializableSingleton<SINGLETON_CLASS>* GetSerializationWrapper() const
    {
        return const_cast<SerializableSingleton<SINGLETON_CLASS>*>(this);
    }

private:
    friend class boost::serialization::access;

    /**
     * Save the wrapped singleton.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        SINGLETON_CLASS* p_instance = SINGLETON_CLASS::Instance();
        archive & *p_instance;
        archive & p_instance;
    }

    /**
     * Load the wrapped singleton.
     *
     * @param archive the archive
     * @param version the saved version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        SINGLETON_CLASS* p_instance = SINGLETON_CLASS::Instance();
        archive & *p_instance;
        archive & p_instance;

        // To avoid a memory leak, we now need to delete this (temporary) wrapper, and
        // ideally tell the archive that it should be using the instance's instead.
        // Unfortunately reset_object_address doesn't work from within load(), so we
        // just create a new wrapper each time; not too bad.
        //archive.reset_object_address(p_instance->GetSerializationWrapper(), this);
        delete this;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

#endif // SERIALIZABLESINGLETON_HPP_
