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

#ifndef FORTESTARCHIVING_HPP_
#define FORTESTARCHIVING_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "ClassIsAbstract.hpp"
#include "Identifiable.hpp"

class BaseClass
{
public:
    unsigned mTagInBaseClass;
    BaseClass();
    virtual ~BaseClass();
    virtual void Hello()=0;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mTagInBaseClass;
    }
};

CLASS_IS_ABSTRACT(BaseClass)

// See http://www.boost.org/libs/serialization/doc/index.html
class ParentClass;

class ChildClass : public BaseClass
{
public:
    unsigned mTag;
    ParentClass* mpParent;
    ChildClass();
    void SetParent(ParentClass* pParent);

    void Hello();

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<BaseClass>(*this);
        archive & mTag;
    }
};

class SubChildClass : public ChildClass
{
public:
    unsigned mSubTag;

    SubChildClass();

    void Hello();

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<ChildClass>(*this);
        archive & mSubTag;
    }
};

/**
 * This is an identifiable child class which has not been registered for serialization.
 * It is to cover the warnings given when a class cannot be identified.
 */
class BadIdentifiable : public Identifiable
{
public:

    BadIdentifiable() : Identifiable()
    {

    }

    ~BadIdentifiable()
    {

    }

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ChildClass)
CHASTE_CLASS_EXPORT(SubChildClass)

class ParentClass
{
public:
    unsigned mTag;
    ChildClass* mpChild;
    ParentClass(ChildClass* pChild);

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mTag;
    }
};

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance of ParentClass.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const ParentClass * t, const unsigned int file_version)
{
    ar << t->mpChild;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a ParentClass instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, ParentClass * t, const unsigned int file_version)
{
    /*
     * It doesn't actually matter what values we pass to our standard
     * constructor, provided they are valid parameter values, since the
     * state loaded later from the archive will overwrite their effect
     * in this case.
     */

    // Invoke inplace constructor to initialize instance of ParentClass.
    ChildClass* p_child;
    ar >> p_child;
    ::new(t)ParentClass(p_child);
}
}
} // namespace ...


#endif /* FORTESTARCHIVING_HPP_ */
