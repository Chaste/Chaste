/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef ABSTRACTCELLPOPULATIONWRITER_HPP_
#define ABSTRACTCELLPOPULATIONWRITER_HPP_

#include "OutputFileHandler.hpp"

// All cell populations
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <string>

template<unsigned SPACE_DIM>
class NodeBasedCellPopualtion;

/**
 * A class written using the visitor pattern for writing node location from a cell population to file.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractCellPopulationWriter
{
private:

    /** The directory in which to write the file. */
    std::string mDirectory;

protected:

    /** The name of node location files */
    std::string mFileName;

    /** An output stream for writing data */
    out_stream mpOutStream;

public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
    AbstractCellPopulationWriter(std::string directory);

    /** A virtual destructor */
    virtual ~AbstractCellPopulationWriter();

    /** Close the stream file */
    void CloseFile();

    /**
     * Open the out stream for writing
     *
     * @param cleanOutputDirectory whether to clean the output directory before writing.
     */
    virtual void OpenOutputFile(bool cleanOutputDirectory = true);

    /**
     * Open the out stream for appending.
     */
    void OpenOutputFileForAppend();

    /**
     * Write the current time stamp to the file.
     */
    void WriteTimeStamp();

    /**
     * Add a newline character to the stream.
     */
    void WriteNewline();

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the mesh based cell population population to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the mesh based cell population population to visit.
     */
    virtual void Visit(MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the mutliple-ca based cell population population to visit.
     */
    virtual void Visit(MultipleCaBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the node based cell population population to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the node based cell population population to visit.
     */
    virtual void Visit(NodeBasedCellPopulationWithBuskeUpdate<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the node based cell population population to visit.
     */
    virtual void Visit(NodeBasedCellPopulationWithParticles<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the potts based cell population population to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the vertex based cell population population to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;
};

#endif /*ABSTRACTCELLPOPULATIONWRITER_HPP_*/
