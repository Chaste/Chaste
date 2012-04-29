/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef CELLMLLOADER_HPP_
#define CELLMLLOADER_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>

#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "AbstractCardiacCell.hpp"
#include "CellMLToSharedLibraryConverter.hpp"

#ifdef CHASTE_CVODE
#include "AbstractCvodeCell.hpp"
#endif

/**
 * A helper class which will dynamically load a cellML file and
 * provide a method to get a pointer to an AbstractCardiacCell
 * (or AbstractCvodeCell if CVODE is enabled).
 */
class CellMLLoader
{
private:
    /** The location of the cellml file to convert */
    FileFinder mCellMLFile;

    /** The location of an output folder to put the converted file and shared library in */
    OutputFileHandler mOutputFileHandler;

    /** A vector of options to be passed to the cellML converter (PyCML) e.g. "--expose-annotated-variables" */
    std::vector<std::string> mOptions;

    /** The converter we will use */
    boost::shared_ptr<CellMLToSharedLibraryConverter> mpConverter;

    /**
     * A method to make a new cell model which is then cast in the public methods to the relevant type.
     * It copies the cellML file to the #mOutputFileHandler location,
     * writes an options file, converts to .hpp and .cpp, makes a shared library, and loads the cell.
     *
     * @param makeCvodeCell  Whether this cell should be an AbstractCvodeCell (false = AbstractCardiacCell)
     * @return a pointer to an AbstractCardiacCellInterface which can be cast as required.
     */
    AbstractCardiacCellInterface* LoadCellMLFile(bool makeCvodeCell);

public:
    /**
     * Constructor
     *
     * @param rCellMLFile  The location of a CellML file to load on the fly
     * @param rOutputFileHandler  An output directory to do the conversion in
     * @param options  Any options to be passed to PyCML e.g. "--expose-annotated-variables"
     */
    CellMLLoader(const FileFinder& rCellMLFile, const OutputFileHandler& rOutputFileHandler, const std::vector<std::string>& options);

    /**
     * Make an AbstractCardiacCell, makes a (Forward)Euler solver, and uses a default CellML stimulus if present.
     * @return a pointer to the cell
     */
    boost::shared_ptr<AbstractCardiacCell> LoadCardiacCellFromCellML(void);

#ifdef CHASTE_CVODE
    /**
     * Make an AbstractCvodeCell, uses a default CellML stimulus if present.
     * @return a pointer to the cell
     */
    boost::shared_ptr<AbstractCvodeCell> LoadCvodeCellFromCellML(void);
#endif
};

#endif // CELLMLLOADER_HPP_
