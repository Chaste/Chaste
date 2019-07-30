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

#include "CellMLLoader.hpp"

#include <algorithm>
#include "EulerIvpOdeSolver.hpp"

CellMLLoader::CellMLLoader(const FileFinder& rCellMLFile, const OutputFileHandler& rOutputFileHandler, const std::vector<std::string>& rOptions)
    : mCellMLFile(rCellMLFile),
      mOutputFileHandler(rOutputFileHandler),
      mOptions(rOptions),
      mpConverter(new CellMLToSharedLibraryConverter(true)), // Builds using just the 'heart' component
      mUseCvode(boost::logic::indeterminate)
{
}

AbstractCardiacCellInterface* CellMLLoader::LoadCellMLFile(bool makeCvodeCell)
{
    std::string model_name = mCellMLFile.GetLeafNameNoExtension();
    FileFinder copied_model = mOutputFileHandler.FindFile(mCellMLFile.GetLeafName());

    // If this is the first call, set up to do the conversion
    if (boost::logic::indeterminate(mUseCvode))
    {
        mUseCvode = makeCvodeCell;

        // Remove or add the "--cvode" option to get the desired cell type
        std::vector<std::string>::iterator it = std::find(mOptions.begin(), mOptions.end(), "--cvode");
        if (!makeCvodeCell && it != mOptions.end())
        {
            mOptions.erase(it);
        }
        if (makeCvodeCell && it == mOptions.end())
        {
            mOptions.push_back("--cvode");
        }

        // Create an options file and put it in the output directory with the CellML file
        mpConverter->CreateOptionsFile(mOutputFileHandler, model_name, mOptions);
        mOutputFileHandler.CopyFileTo(mCellMLFile);
        // Copy a .out file also if present
        FileFinder out_file(model_name + ".out", mCellMLFile);
        if (out_file.Exists())
        {
            mOutputFileHandler.CopyFileTo(out_file);
        }
    }
    // If however we've made a cell before, check that we're making the same type this time
    else if (makeCvodeCell != mUseCvode)
    {
        EXCEPTION("You cannot call both LoadCvodeCell and LoadCardiacCell on the same CellMLLoader.");
    }

    // Convert the CellML to a shared library (no-op if shared library exists)
    DynamicCellModelLoaderPtr p_loader = mpConverter->Convert(copied_model);

    // Use the shared library to load a concrete cell
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    boost::shared_ptr<EulerIvpOdeSolver> p_solver;
    if (!makeCvodeCell)
    {
        // Put in a Forward Euler solver as a default for normal cells
        p_solver.reset(new EulerIvpOdeSolver);
    }
    AbstractCardiacCellInterface* p_loaded_cell = p_loader->CreateCell(p_solver, p_stimulus);

    // If the CellML file has a default stimulus we may as well use it.
    if (p_loaded_cell->HasCellMLDefaultStimulus())
    {
        p_loaded_cell->UseCellMLDefaultStimulus();
    }
    return p_loaded_cell;
}

boost::shared_ptr<AbstractCardiacCell> CellMLLoader::LoadCardiacCell(void)
{
    AbstractCardiacCellInterface* p_loaded_cell = LoadCellMLFile(false);
    boost::shared_ptr<AbstractCardiacCell> p_model(dynamic_cast<AbstractCardiacCell*>(p_loaded_cell));
    return p_model;
}

#ifdef CHASTE_CVODE
boost::shared_ptr<AbstractCvodeCell> CellMLLoader::LoadCvodeCell(void)
{
    AbstractCardiacCellInterface* p_loaded_cell = LoadCellMLFile(true);
    boost::shared_ptr<AbstractCvodeCell> p_model(dynamic_cast<AbstractCvodeCell*>(p_loaded_cell));
    return p_model;
}
#endif


