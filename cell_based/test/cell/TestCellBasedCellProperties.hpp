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

#ifndef TESTCELLBASEDCELLPROPERTIES_HPP_
#define TESTCELLBASEDCELLPROPERTIES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "CellId.hpp"
#include "CellData.hpp"

#include "CellPropertyRegistry.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCellBasedCellProperties : public AbstractCellBasedTestSuite
{
public:

    void TestCellIdMethods()
    {
        // Resetting the Maximum cell Id to zero (to account for previous tests) is done in the Setup method...
        // CellId::ResetMaxCellId();

        MAKE_PTR(CellId, p_cell_id);

        TS_ASSERT_THROWS_THIS(p_cell_id->GetCellId(), "AssignCellId must be called before using the CellID");
        TS_ASSERT_THROWS_THIS(p_cell_id->GetMaxCellId(), "AssignCellId must be called before using the CellID");

        // The ID is not assigned until this method is called.
        p_cell_id->AssignCellId();

        TS_ASSERT_EQUALS(p_cell_id->GetCellId(), 0u);
        TS_ASSERT_EQUALS(p_cell_id->GetMaxCellId(), 1u);
    }

    void TestArchiveCellId()
    {
        MAKE_PTR(CellId, p_extra_cell_id);
        p_extra_cell_id->AssignCellId();
        // In this test the Max Cell Id starts at 2 because we just made a cell and assigned to it

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell_id.arch";

        // Archive Cell ID
        {
            CellId* p_cell_id = new CellId();
            p_cell_id->AssignCellId();

            TS_ASSERT_EQUALS(p_cell_id->GetCellId(), 1u);
            TS_ASSERT_EQUALS(p_cell_id->GetMaxCellId(), 2u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_cell_id = p_cell_id;
            output_arch << p_const_cell_id;

            delete p_cell_id;
        }

        // Restore cell ID
        {
            AbstractCellProperty* p_cell_id;

            // Restore the Cell ID
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_id;

            CellId* p_real_cell_id = dynamic_cast<CellId*>(p_cell_id);
            TS_ASSERT(p_real_cell_id != NULL);

            TS_ASSERT_EQUALS(p_real_cell_id->GetCellId(), 1u);
            TS_ASSERT_EQUALS(p_real_cell_id->GetMaxCellId(), 2u);

            // Tidy up
            delete p_cell_id;
        }
    }

    void TestCellDataMethods()
    {
        MAKE_PTR(CellData, p_cell_data);

        TS_ASSERT_THROWS_THIS(p_cell_data->GetItem("thing1"), "The item thing1 is not stored");

        p_cell_data->SetItem("thing1", 1.0);
        p_cell_data->SetItem("thing2", 2.0);
        p_cell_data->SetItem("thing3", 3.0);

        TS_ASSERT_DELTA(p_cell_data->GetItem("thing1"), 1.0, 1e-8);
        TS_ASSERT_DELTA(p_cell_data->GetItem("thing2"), 2.0, 1e-8);
        TS_ASSERT_DELTA(p_cell_data->GetItem("thing3"), 3.0, 1e-8);
        TS_ASSERT_EQUALS(p_cell_data->GetNumItems(), 3u);
    }

    void TestArchiveCellData()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "cell_data.arch";

        // Archive Cell data
        {
            CellData* p_cell_data = new CellData;

            p_cell_data->SetItem("thing1",1.0);
            p_cell_data->SetItem("thing2",2.0);

            TS_ASSERT_DELTA(p_cell_data->GetItem("thing1"), 1.0, 1e-8);
            TS_ASSERT_DELTA(p_cell_data->GetItem("thing2"), 2.0, 1e-8);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_cell_data = p_cell_data;
            output_arch << p_const_cell_data;

            delete p_cell_data;
        }

        // Restore Cell data
        {
            AbstractCellProperty* p_cell_data;

            // Restore the Cell data
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_data;

            CellData* p_real_cell_data = dynamic_cast<CellData*>(p_cell_data);
            TS_ASSERT(p_real_cell_data != NULL);

            TS_ASSERT_DELTA(p_real_cell_data->GetItem("thing1"), 1.0, 1e-8);
            TS_ASSERT_DELTA(p_real_cell_data->GetItem("thing2"), 2.0, 1e-8);

            // Tidy up
            delete p_cell_data;
        }
    }
};

#endif /* TESTCELLBASEDCELLPROPERTIES_HPP_ */
