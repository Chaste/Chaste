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
class SerializableSingleton : private boost::noncopyable
{
public:
    /**
     * @return the wrapper object to use to serialize the related singleton instance.
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
