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

#ifndef TESTARCHIVING_HPP_
#define TESTARCHIVING_HPP_

#include "CheckpointArchiveTypes.hpp"

#include <fstream>
#include <cxxtest/TestSuite.h>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "OutputFileHandler.hpp"
#include "Warnings.hpp"
#include "ClassOfSimpleVariables.hpp"
#include "ForTestArchiving.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * Contains some good examples of how to checkpoint things.
 * Most of these would need changing to work in parallel, we have some easy ways of doing this:
 * see https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/BoostSerialization for details.
 */
class TestArchiving : public CxxTest::TestSuite
{
public:
    void TestArchiveSimpleVars()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "simple_vars.arch";

        // Create an output archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            std::vector<double> doubles(3);
            doubles[0] = 1.1;
            doubles[1] = 1.2;
            doubles[2] = 1.3;

            std::vector<bool> bools(2);
            bools[0] = true;
            bools[1] = true;

            ClassOfSimpleVariables i(42,"hello",doubles,bools);

            // Cast to const
            output_arch << static_cast<const ClassOfSimpleVariables&>(i);
        }

        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            std::vector<double> bad_doubles(1);
            bad_doubles[0] = 10.3;

            std::vector<bool> bad_bools(1);
            bad_bools[0] = false;

            ClassOfSimpleVariables j(0,"bye",bad_doubles,bad_bools);

            // Read the archive
            input_arch >> j;

            // Check that the values are correct
            TS_ASSERT_EQUALS(j.GetNumber(),42);
            TS_ASSERT_EQUALS(j.GetString(),"hello");
            TS_ASSERT_EQUALS(j.GetVectorOfDoubles().size(),3u);
            TS_ASSERT_EQUALS(j.GetVectorOfBools().size(),2u);

            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[0],1.1,1e-12);
            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[1],1.2,1e-12);
            TS_ASSERT_DELTA(j.GetVectorOfDoubles()[2],1.3,1e-12);

            TS_ASSERT(j.GetVectorOfBools()[0]);
            TS_ASSERT(j.GetVectorOfBools()[1]);
        }
    }

    void TestArchivingLinkedChildAndParent()
    {
        /*
         * This test is an abstraction of archiving a cyclically linked parent-child pair.
         * The parent represents a Cell and the child represents an AbstractCellCycleModel.
         */

        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "linked_classes.arch";

        // Save
        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ChildClass* p_child = new ChildClass;
            ParentClass* p_parent = new ParentClass(p_child);

            p_child->mTag = 11;
            p_parent->mTag = 10;

            ParentClass* const p_parent_for_archiving = p_parent;
            //ChildClass* const p_child_for_archiving = p_child;

            //output_arch << p_child_for_archiving;
            output_arch << p_parent_for_archiving;

            delete p_child;
            delete p_parent;
        }

        // Load
        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            //ChildClass* p_child;
            ParentClass* p_parent;

            input_arch >> p_parent;

            TS_ASSERT_EQUALS(p_parent->mTag, 10u);
            TS_ASSERT_EQUALS(p_parent->mpChild->mTag, 11u);
            TS_ASSERT_EQUALS(p_parent->mpChild->mpParent, p_parent);

            // Tidy up
            delete p_parent->mpChild;
            delete p_parent;
        }
    }

    void TestArchivingSetOfSetOfPointers()
    {
        /*
         * This test is an abstraction of archiving a set of sets of pointers and a list of objects.
         * Note that the list.push_back method uses the copy constructor. This is why we iterate
         * through the list to generate the pointers to populate the set.
         */

        std::vector<double> doubles;
        std::vector<bool> bools;

        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "pointer_set.arch";

        // Save
        {
            // Create aClassOfSimpleVariablesn output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ClassOfSimpleVariables one(42, "hello", doubles,bools);
            ClassOfSimpleVariables two(256, "goodbye", doubles,bools);
            ClassOfSimpleVariables three(1, "not used in set", doubles,bools);

            std::list<ClassOfSimpleVariables> a_list;
            std::set<ClassOfSimpleVariables*> a_set;
            a_list.push_back(one);
            a_set.insert( &(a_list.back()) );
            a_list.push_back(two);
            a_set.insert( &(a_list.back()) );
            a_list.push_back(three);

            std::set<std::set<ClassOfSimpleVariables*> > wrapper_set;
            wrapper_set.insert(a_set);

            output_arch << static_cast<const std::list<ClassOfSimpleVariables>&>(a_list);
            output_arch << static_cast<const std::set<std::set<ClassOfSimpleVariables*> >&>(wrapper_set);
        }

        // Load
        {
            std::set<std::set<ClassOfSimpleVariables*> > wrapper_set;
            std::list<ClassOfSimpleVariables> a_list;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            TS_ASSERT_EQUALS(wrapper_set.size(), 0u);
            input_arch >> a_list;
            input_arch >> wrapper_set;
            TS_ASSERT_EQUALS(wrapper_set.size(), 1u);
            const std::set<ClassOfSimpleVariables*>& a_set = *(wrapper_set.begin());
            TS_ASSERT_EQUALS(a_set.size(), 2u);

            ClassOfSimpleVariables* p_one_in_set = NULL;
            ClassOfSimpleVariables* p_two_in_set = NULL;
            for (std::set<ClassOfSimpleVariables*>::iterator it = a_set.begin();
                 it != a_set.end();
                 ++it)
            {
                   ClassOfSimpleVariables* p_class = *(it);
                   if (p_class->GetNumber() == 42)
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 42);
                        TS_ASSERT_EQUALS(p_class->GetString(), "hello");
                        p_one_in_set = p_class;
                   }
                   else
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 256);
                        TS_ASSERT_EQUALS(p_class->GetString(), "goodbye");
                        p_two_in_set = p_class;
                   }
            }

            ClassOfSimpleVariables* p_one_in_list = NULL;
            ClassOfSimpleVariables* p_two_in_list = NULL;
            for (std::list<ClassOfSimpleVariables>::iterator it = a_list.begin();
                 it != a_list.end();
                 ++it)
            {
                   ClassOfSimpleVariables* p_class = &(*it);
                   if (p_class->GetNumber() == 42)
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 42);
                        TS_ASSERT_EQUALS(p_class->GetString(), "hello");
                        p_one_in_list = p_class;
                   }
                   else if (p_class->GetNumber() == 256)
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 256);
                        TS_ASSERT_EQUALS(p_class->GetString(), "goodbye");
                        p_two_in_list = p_class;
                   }
                   else
                   {
                        TS_ASSERT_EQUALS(p_class->GetNumber(), 1);
                        TS_ASSERT_EQUALS(p_class->GetString(), "not used in set");
                   }
            }
            TS_ASSERT_DIFFERS(p_one_in_list, (void*)NULL);
            TS_ASSERT_DIFFERS(p_two_in_list, (void*)NULL);
            TS_ASSERT_EQUALS(p_one_in_list, p_one_in_set);
            TS_ASSERT_EQUALS(p_two_in_list, p_two_in_set);
        }
    }

    void TestArchivingBoostSharedPtrToChild()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "shared_ptr.arch";

        // Save
        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            boost::shared_ptr<ChildClass> p_child(new ChildClass);

            p_child->mTag = 11;
            p_child->mTagInBaseClass = 3;

            boost::shared_ptr<ChildClass> const p_child_for_archiving = p_child;

            output_arch << p_child_for_archiving;
        }

        // Load
        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            boost::shared_ptr<ChildClass> p_child;

            input_arch >> p_child;

            TS_ASSERT_EQUALS(p_child->mTag, 11u);
            TS_ASSERT_EQUALS(p_child->mTagInBaseClass, 3u);
        }
    }

    void TestArchivingBoostSharedPtrToChildUsingBaseClass()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "shared_ptr_abs.arch";

        // Save
        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            boost::shared_ptr<BaseClass> p_base(new ChildClass());

            p_base->mTagInBaseClass = 6;

            output_arch << p_base;
        }

        // Load
        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            boost::shared_ptr<BaseClass> p_base;

            input_arch >> p_base;

            TS_ASSERT_EQUALS(p_base->mTagInBaseClass, 6u);
        }
    }

    void TestArchivingSubChild()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "subchild.arch";

        // Save
        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            boost::shared_ptr<BaseClass> p_base(new SubChildClass());

            p_base->mTagInBaseClass = 6;

            output_arch << p_base;
        }

        // Load
        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            boost::shared_ptr<BaseClass> p_base;

            input_arch >> p_base;

            TS_ASSERT_EQUALS(p_base->mTagInBaseClass, 6u);
            TS_ASSERT(dynamic_cast<ChildClass*>(p_base.get()));
            TS_ASSERT(dynamic_cast<SubChildClass*>(p_base.get()));
        }
    }

    /**
     * HOW_TO_TAG General/Archiving
     * Use a binary rather than ascii boost archive format, for speed and smaller file sizes.
     *
     * (NB this file is actually larger, but when there's lots of data it should be smaller!)
     *
     * Note that a binary archive could become machine architecture-specific.
     * But it can be helpful to use them for speed. In this case we have had
     * success in storing an ascii one, loading and re-saving as binary.
     * Loading code could then look for the binary one, and revert to ascii if
     * it is not there.
     *
     * The test below is identical to the one above, apart from the two lines indicated.
     */
    void TestUsingABinaryArchive()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "subchild_binary.arch";

        // Save
        {
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str(), std::ios::binary);
            boost::archive::binary_oarchive output_arch(ofs); // LINE CHANGED!

            boost::shared_ptr<BaseClass> p_base(new SubChildClass());

            p_base->mTagInBaseClass = 6;

            output_arch << p_base;
        }

        // Load
        {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::binary_iarchive input_arch(ifs); // LINE CHANGED!

            boost::shared_ptr<BaseClass> p_base;

            input_arch >> p_base;

            TS_ASSERT_EQUALS(p_base->mTagInBaseClass, 6u);
            TS_ASSERT(dynamic_cast<ChildClass*>(p_base.get()));
            TS_ASSERT(dynamic_cast<SubChildClass*>(p_base.get()));
        }
    }

    void TestUndentifiableClass()
    {
        // This Identifiable child class is not registered for serialization, it is expected to give a warning when GetIdentifiable is called.
        BadIdentifiable bad_indentifiable;

        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        std::string class_name = bad_indentifiable.GetIdentifier();
        TS_ASSERT_EQUALS(class_name, "UnknownClass-ReflectionFailed");
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);

        // Clean up
        Warnings::QuietDestroy();
    }
};

#endif /*TESTARCHIVING_HPP_*/
