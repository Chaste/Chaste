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

#ifndef TESTDYNAMICALLYLOADEDCELLMODELS_HPP_
#define TESTDYNAMICALLYLOADEDCELLMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/assign.hpp>
#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveLocationInfo.hpp"
#endif // CHASTE_CAN_CHECKPOINT_DLLS

#include "RunAndCheckIonicModels.hpp"
#include "DynamicLoadingHelperFunctions.hpp"

#include "DynamicCellModelLoader.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "CellMLLoader.hpp"
#include "CellMLToSharedLibraryConverter.hpp"

#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ChasteBuildRoot.hpp"
#include "HeartConfig.hpp"
#include "FileFinder.hpp"
#include "HeartFileFinder.hpp"
#include "ChasteSyscalls.hpp"

#include "AbstractDynamicallyLoadableEntity.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "AbstractCvodeCell.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDynamicallyLoadedCellModels : public CxxTest::TestSuite
{
private:

    void RunLr91Test(DynamicCellModelLoader& rLoader,
                     unsigned vIndex=4u,
                     bool testTables=false,
                     double tolerance=1e-3,
                     double tableTestV=-100000)
    {
        AbstractCardiacCellInterface* p_cell = CreateLr91CellFromLoader(rLoader, vIndex);
        SimulateLr91AndCompare(p_cell, tolerance);

        if (testTables)
        {
            double v = p_cell->GetVoltage();
            p_cell->SetVoltage(tableTestV);
            TS_ASSERT_THROWS_CONTAINS(p_cell->GetIIonic(), "membrane_voltage outside lookup table range");
            p_cell->SetVoltage(v);
        }

        delete p_cell;
    }

    void SimulateLr91AndCompare(AbstractCardiacCellInterface* pCell,
                                double tolerance=1e-3)
    {
        double end_time = 1000.0; //One second in milliseconds
        // Solve and write to file
        clock_t ck_start = clock();

        // Don't use RunOdeSolverWithIonicModel() as this is hardcoded to AbstractCardiacCells
        OdeSolution solution1 = pCell->Compute(0.0, end_time);
        solution1.WriteToFile("TestIonicModels", "DynamicallyLoadableLr91", "ms", 100, false, 4);

        clock_t ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tSolve time: " << forward << std::endl;

        // Compare with 'normal' LR91 model results
        CheckCellModelResults("DynamicallyLoadableLr91", "Lr91DelayedStim", tolerance);

        // Test GetIIonic against hardcoded result from TestIonicModels.hpp
        // Don't use RunOdeSolverWithIonicModel() as this is hardcoded to AbstractCardiacCells
        OdeSolution solution2 = pCell->Compute(0.0, 60.0);
        solution2.WriteToFile("TestIonicModels", "DynamicallyLoadableLr91GetIIonic", "ms", 100, false, 4);

        // For coverage
        TS_ASSERT_DELTA(pCell->GetIntracellularCalciumConcentration(), 0.0012, 1e-4);
        TS_ASSERT_DELTA(pCell->GetIIonic(), 1.9411, tolerance);
    }

    AbstractCardiacCellInterface* CreateLr91CellFromLoader(DynamicCellModelLoader& rLoader,
                                                           unsigned vIndex=4u)
    {
        AbstractCardiacCellInterface* p_cell = CreateCellWithStandardStimulus(rLoader);
        TS_ASSERT_EQUALS(p_cell->GetVoltageIndex(), vIndex);
        return p_cell;
    }

    mode_t ResetMode(FileFinder& rFile, mode_t mode)
    {
        struct stat our_stats;
        int retcode = stat(rFile.GetAbsolutePath().c_str(), &our_stats);
        EXCEPT_IF_NOT(retcode == 0);
        chmod(rFile.GetAbsolutePath().c_str(), mode);
        return our_stats.st_mode;
    }

    std::string mArchivingDirName;
    FileFinder mArchivingModel;
    void SaveSoForArchivingTest(FileFinder& rSoFile)
    {
        mArchivingDirName = "TestDynamicallyLoadedCellModelsArchiving";
        OutputFileHandler handler(mArchivingDirName);
        mArchivingModel = handler.CopyFileTo(rSoFile);
    }

public:
    /**
     * This test demonstrates the easiest way to load a single cell model from CellML,
     * using the CellMLLoader class.
     */
    void TestCellmlLoaderClass()
    {
        FileFinder cellml_file("heart/src/odes/cellml/LuoRudy1991.cellml", RelativeTo::ChasteSourceRoot);
        // Stimulus to use for simulation, so it matches other tests in this suite
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new SimpleStimulus(-25.5, 2.0, 50.0));
        {
            OutputFileHandler handler("TestCardiacCellMLLoader");
            // Note that the --cvode flag will be ignored since we call the LoadCardiacCell method.
            std::vector<std::string> options = boost::assign::list_of("--cvode")("--expose-annotated-variables");
            CellMLLoader loader(cellml_file, handler, options);
            boost::shared_ptr<AbstractCardiacCell> p_cell = loader.LoadCardiacCell();
            TS_ASSERT_EQUALS(p_cell->GetSystemName(), "luo_rudy_1991");
            p_cell->SetStimulusFunction(p_stimulus);
            SimulateLr91AndCompare(p_cell.get());

            // Are sources from the conversion preserved?
            FileFinder cpp_file(handler.GetRelativePath() + "/LuoRudy1991.cpp", RelativeTo::ChasteTestOutput);
            TS_ASSERT(cpp_file.Exists());
            TS_ASSERT(cpp_file.IsNewerThan(cellml_file));
            FileFinder hpp_file(handler.GetRelativePath() + "/LuoRudy1991.hpp", RelativeTo::ChasteTestOutput);
            TS_ASSERT(hpp_file.Exists());
            TS_ASSERT(hpp_file.IsNewerThan(cellml_file));

           // We can't now call LoadCvodeCell on this loader
#ifdef CHASTE_CVODE
            TS_ASSERT_THROWS_THIS(loader.LoadCvodeCell(),
                                  "You cannot call both LoadCvodeCell and LoadCardiacCell on the same CellMLLoader.");
#endif
        }
#ifdef CHASTE_CVODE
        {
            OutputFileHandler handler("TestCvodeCellMLLoader");
            std::vector<std::string> options = boost::assign::list_of("--expose-annotated-variables");
            CellMLLoader loader(cellml_file, handler, options);
            boost::shared_ptr<AbstractCvodeCell> p_cell = loader.LoadCvodeCell();
            TS_ASSERT_EQUALS(p_cell->GetSystemName(), "luo_rudy_1991");
            p_cell->SetStimulusFunction(p_stimulus);
            SimulateLr91AndCompare(p_cell.get(), 1.0); // Large tolerance due to different ODE solver

            // We can't now call LoadCardiacCell on this loader
            TS_ASSERT_THROWS_THIS(loader.LoadCardiacCell(),
                                  "You cannot call both LoadCvodeCell and LoadCardiacCell on the same CellMLLoader.");
        }
#endif
    }

    /**
     * This is based on TestOdeSolverForLR91WithDelayedSimpleStimulus from
     * TestIonicModels.hpp.
     */
    void TestDynamicallyLoadedLr91()
    {
        // Load the cell model dynamically
        std::string model_name = "libDynamicallyLoadableLr91.";
        model_name  += CellMLToSharedLibraryConverter::msSoSuffix;
        // All tests use the registry, as not doing so can lead to segfaults...
        DynamicCellModelLoaderPtr p_loader = DynamicModelLoaderRegistry::Instance()->GetLoader(
            ChasteComponentBuildDir("heart") + "dynamic/" + model_name);
        RunLr91Test(*p_loader);

        // The .so also gets copied into the source folder
        DynamicCellModelLoaderPtr p_loader2 = DynamicModelLoaderRegistry::Instance()->GetLoader(
            std::string(ChasteBuildRootDir()) + "heart/dynamic/" + model_name);
        RunLr91Test(*p_loader2);
    }

    void TestExceptions()
    {
        // Try loading a .so that doesn't exist
        std::string file_name = "non-existent-file-we-hope";
        TS_ASSERT_THROWS_CONTAINS(DynamicCellModelLoader::Create(file_name),
                                  "Unable to load .so file '" + file_name + "':");

        // Try loading a .so that doesn't define a cell model
        file_name = "libNotACellModel.";
        file_name += CellMLToSharedLibraryConverter::msSoSuffix;
        TS_ASSERT_THROWS_CONTAINS(DynamicCellModelLoader::Create(ChasteComponentBuildDir("heart") + "dynamic/" + file_name),
                                  "Failed to load cell creation function from .so file");
    }

    void TestCellmlConverter()
    {
        // Copy CellML file into output dir
        std::string dirname = "TestCellmlConverter";
        OutputFileHandler handler(dirname);
        FileFinder cellml_file_src("heart/dynamic/luo_rudy_1991_dyn.cellml", RelativeTo::ChasteSourceRoot);

        CellMLToSharedLibraryConverter converter(true);

        // Convert a real CellML file
        FileFinder cellml_file = handler.CopyFileTo(cellml_file_src);
        TS_ASSERT(cellml_file.Exists());
        FileFinder so_file(dirname + "/libluo_rudy_1991_dyn."+CellMLToSharedLibraryConverter::msSoSuffix, RelativeTo::ChasteTestOutput);
        TS_ASSERT(!so_file.Exists());
        DynamicCellModelLoaderPtr p_loader = converter.Convert(cellml_file);
        SaveSoForArchivingTest(so_file);
        TS_ASSERT(so_file.Exists());
        TS_ASSERT(so_file.IsNewerThan(cellml_file));
        // Converting a .so should be a "no-op"
        DynamicCellModelLoaderPtr p_loader2 = converter.Convert(so_file);
        TS_ASSERT(so_file.Exists());
        TS_ASSERT(p_loader2 == p_loader);
        RunLr91Test(*p_loader, 0u);

        // Cover exceptions
        std::string file_name = "test.no.file";
        FileFinder no_exist(dirname + "/" + file_name, RelativeTo::ChasteTestOutput);
        TS_ASSERT_THROWS_THIS(converter.Convert(no_exist), "Dynamically loadable cell model '"
                              + no_exist.GetAbsolutePath() + "' does not exist.");

        file_name = "test";
        FileFinder no_ext(dirname + "/" + file_name, RelativeTo::ChasteTestOutput);
        PetscTools::Barrier("TestCellmlConverter_pre_touch");
        if (PetscTools::AmMaster())
        {
            out_stream fp = handler.OpenOutputFile(file_name);
            fp->close();
        }
        PetscTools::Barrier("TestCellmlConverter_post_touch");

        TS_ASSERT_THROWS_THIS(converter.Convert(no_ext), "File does not have an extension: " + no_ext.GetAbsolutePath());
        FileFinder unsupp_ext("global/src/FileFinder.hpp", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_THIS(converter.Convert(unsupp_ext), "Unsupported extension '.hpp' of file '"
                              + unsupp_ext.GetAbsolutePath() + "'; must be .so, .dylib or .cellml");


        TRY_IF_MASTER( so_file.Remove() );

        TS_ASSERT_THROWS_THIS(converter.Convert(cellml_file, false),
                              "Unable to convert .cellml to .so unless called collectively, due to possible race conditions.");

        // This one is tricky!
        mode_t old_mode = ResetMode(cellml_file, 0);
        TS_ASSERT_THROWS_CONTAINS(converter.Convert(cellml_file),
                                  "Conversion of CellML to Chaste shared object failed.");
        ResetMode(cellml_file, old_mode);

        // What if the Chaste build tree is missing?
        FileFinder::FakePath(RelativeTo::ChasteBuildRoot, "/tmp/not-a-chaste-source-tree");
        TS_ASSERT_THROWS_THIS(converter.Convert(cellml_file),
                              "No Chaste build tree found at '/tmp/not-a-chaste-source-tree' - you need the source to use CellML models directly in Chaste.");
        FileFinder::StopFaking();

        // Or a required project is missing?
        {
            FileFinder build_root("", RelativeTo::ChasteBuildRoot);
            CellMLToSharedLibraryConverter failed_converter(true, "not_a_project");
            TS_ASSERT_THROWS_CONTAINS(failed_converter.Convert(cellml_file),
                "Unable to convert CellML model: required Chaste component 'not_a_project' does not exist in '"
                                      + build_root.GetAbsolutePath() + "'.");
        }

        // Or we can't create the temp folder for some other reason?
        {
            OutputFileHandler fake_chaste_tree("FakeChaste/heart");
            FileFinder::FakePath(RelativeTo::ChasteBuildRoot, fake_chaste_tree.GetChasteTestOutputDirectory() + "FakeChaste");
            TS_ASSERT_THROWS_CONTAINS(converter.Convert(cellml_file),
                                      "Failed to create temporary folder '");
            FileFinder::StopFaking();
        }
    }

    void TestCellmlConverterWithOptions()
    {
        // Copy CellML file into output dir
        std::string dirname = "TestCellmlConverterWithOptions";
        std::string model = "LuoRudy1991";
        OutputFileHandler handler(dirname + "/plain");
        FileFinder cellml_file("heart/src/odes/cellml/" + model + ".cellml", RelativeTo::ChasteSourceRoot);
        FileFinder copied_file = handler.CopyFileTo(cellml_file);

        // Do the conversions preserving generated sources
        CellMLToSharedLibraryConverter converter(true);

        // Create options file & convert
        std::vector<std::string> args;
        args.push_back("--opt");
        converter.CreateOptionsFile(handler, model, args);
        // Ensure that conversion works if CWD != ChasteSourceRoot
        EXPECT0(chdir, "heart");
        DynamicCellModelLoaderPtr p_loader = converter.Convert(copied_file);
        EXPECT0(chdir, "..");
        RunLr91Test(*p_loader, 0u, true, 0.01); // Implementation of lookup tables has improved...
        // Check the sources exist
        TS_ASSERT(handler.FindFile(model + ".cpp").Exists());
        TS_ASSERT(handler.FindFile(model + ".hpp").Exists());

        {
            // Backward Euler
            args[0] = "--backward-euler";
            OutputFileHandler handler2(dirname + "/BE");
            FileFinder copied_file2 = handler2.CopyFileTo(cellml_file);
            FileFinder maple_output_file("heart/src/odes/cellml/LuoRudy1991.out", RelativeTo::ChasteSourceRoot);
            handler2.CopyFileTo(maple_output_file);
            converter.CreateOptionsFile(handler2, model, args);
            p_loader = converter.Convert(copied_file2);
            RunLr91Test(*p_loader, 0u, true, 0.3);
        }
#ifdef CHASTE_CVODE
        {
            // With a for_model section and Cvode
            args[0] = "--opt";
            args.push_back("--cvode");
            OutputFileHandler handler3(dirname + "/CO");
            FileFinder copied_file3 = handler3.CopyFileTo(cellml_file);
            std::string for_model = std::string("<for_model id='luo_rudy_1991'><lookup_tables><lookup_table>")
                    + "<var type='config-name'>transmembrane_potential</var>"
                    + "<max>69.9999</max>"
                    + "</lookup_table></lookup_tables></for_model>\n";
            converter.CreateOptionsFile(handler3, model, args, for_model);
            p_loader = converter.Convert(copied_file3);
            RunLr91Test(*p_loader, 0u, true, 1, 70); // Large tolerance due to different ODE solver
        }
#endif
    }

    void TestArchiving()
    {
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
        // Check the previous test has left us a .so file
        TS_ASSERT(mArchivingModel.Exists());

        // Get a loader for the .so and load a cell model
        CellMLToSharedLibraryConverter converter;
        DynamicCellModelLoaderPtr p_loader = converter.Convert(mArchivingModel);
        AbstractCardiacCellInterface* p_cell = CreateLr91CellFromLoader(*p_loader, 0u);

        // Archive it
        OutputFileHandler handler(mArchivingDirName, false);
        handler.SetArchiveDirectory();
        std::string archive_filename1 = ArchiveLocationInfo::GetProcessUniqueFilePath("first-save.arch");
        {
            AbstractCardiacCellInterface* const p_const_cell = p_cell;
            std::ofstream ofs(archive_filename1.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            ///\todo #2417 this archiving throws exception on Mac OSX
            try
            {
                output_arch << p_const_cell;
            }
            catch(boost::archive::archive_exception& boost_exception)
            {
                TS_ASSERT_EQUALS(boost_exception.code, boost::archive::archive_exception::unregistered_class);
                TS_FAIL("Archiving cell models in unavailable.  Please refer to  #2417");
                //Bail out
                return;
            }
        }

        // Load from archive
        AbstractCardiacCellInterface* p_loaded_cell1;
        {
            std::ifstream ifs(archive_filename1.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_loaded_cell1;
        }

        // Archive the un-archived model
        std::string archive_filename2 = ArchiveLocationInfo::GetProcessUniqueFilePath("second-save.arch");
        {
            AbstractCardiacCellInterface* const p_const_cell = p_loaded_cell1;
            std::ofstream ofs(archive_filename2.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_const_cell;
        }

        // Load from the new archive
        AbstractCardiacCellInterface* p_loaded_cell2;
        {
            std::ifstream ifs(archive_filename2.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_loaded_cell2;
        }

        // Check simulations of both loaded cells
        SimulateLr91AndCompare(p_loaded_cell1);
        delete p_loaded_cell1;
        SimulateLr91AndCompare(p_loaded_cell2);
        delete p_loaded_cell2;
        delete p_cell;
#else
        std::cout << "Note: this test can only actually test anything on Boost>=1.37 (on non-Mac systems #2417)." << std::endl;
#endif // CHASTE_CAN_CHECKPOINT_DLLS
    }
};

#endif /* TESTDYNAMICALLYLOADEDCELLMODELS_HPP_ */
