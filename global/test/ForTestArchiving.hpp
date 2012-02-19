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

#ifndef FORTESTARCHIVING_HPP_
#define FORTESTARCHIVING_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "ClassIsAbstract.hpp"

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

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ChildClass)

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
