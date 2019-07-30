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

#ifndef ABSTRACTCELLBASEDWRITER_HPP_
#define ABSTRACTCELLBASEDWRITER_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "Identifiable.hpp"
#include "OutputFileHandler.hpp"

/**
 * Abstract class for a writer that takes data from an AbstractCellPopulation and writes it to file.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractCellBasedWriter : public Identifiable
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mFileName;
    }

protected:

    /** The name of the output file. */
    std::string mFileName;

    /** An output stream for writing data. */
    out_stream mpOutStream;

public:

    /**
     * Constructor.
     *
     * @param rFileName the name of the file to write to.
     */
    AbstractCellBasedWriter(const std::string& rFileName);

    /**
     * Virtual destructor.
     */
    virtual ~AbstractCellBasedWriter();

    /**
     * Close mpOutStream.
     */
    void CloseFile();

    /**
     * Open mpOutStream for writing.
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     */
    virtual void OpenOutputFile(OutputFileHandler& rOutputFileHandler);

    /**
     * Open mpOutStream for appending.
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     *
     */
    void OpenOutputFileForAppend(OutputFileHandler& rOutputFileHandler);

    /**
     * Write the current time stamp to mpOutStream.
     */
    virtual void WriteTimeStamp();

    /**
     * Add a newline character to mpOutStream.
     */
    virtual void WriteNewline();

    /**
     * Set the output file name.
     * This method allows the user to change mFileName from
     * its default value, which is set in each subclass's
     * constructor.
     *
     * @param fileName the output file name
     */
    void SetFileName(std::string fileName);

    /**
     * @return the output file name.
     */
    std::string GetFileName();
};

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractCellBasedWriter)

#endif /*ABSTRACTCELLBASEDWRITER_HPP_*/
