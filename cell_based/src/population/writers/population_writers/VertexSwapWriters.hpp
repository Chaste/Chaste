/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef VERTEXSWAPWRITERS_HPP_
#define VERTEXSWAPWRITERS_HPP_

#include "OutputFileHandler.hpp"
#include "AbstractCellPopulation.hpp"

#include "AbstractCellPopulationWriter.hpp"

#include <string>

/** A class written using the visitor pattern for writing the location of T3 swaps to file. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexT3SwapLocationsWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{

public:

    /**
     * Default constructor.
     *
     * @param directory the path to the directory in to which this class should write.
     */
    VertexT3SwapLocationsWriter(std::string directory)
        : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>(directory)
    {
        this->mFileName = "T3SwapLocations.dat";
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
        std::vector< c_vector<double, SPACE_DIM> > t3_swap_locations = pCellPopulation->rGetMesh().GetLocationsOfT3Swaps();

        *this->mpOutStream << t3_swap_locations.size() << "\t";
        for (unsigned index = 0;  index < t3_swap_locations.size(); index++)
        {
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                *this->mpOutStream << t3_swap_locations[index][i] << "\t";
            }
        }

        pCellPopulation->rGetMesh().ClearLocationsOfT3Swaps();
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the MultipleCaBasedCellPopulation to visit.
     */
    virtual void Visit(MultipleCaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
    }
};

/** Writer for T1 swap locations.*/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexT1SwapLocationsWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{

public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
    VertexT1SwapLocationsWriter(std::string directory)
        : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>(directory)
    {
            this->mFileName = "T1SwapLocations.dat";
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
        std::vector< c_vector<double, SPACE_DIM> > t1_swap_locations = pCellPopulation->rGetMesh().GetLocationsOfT1Swaps();

        *this->mpOutStream << t1_swap_locations.size() << "\t";
        for (unsigned index = 0;  index < t1_swap_locations.size(); index++)
        {
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                *this->mpOutStream << t1_swap_locations[index][i] << "\t";
            }
        }

        pCellPopulation->rGetMesh().ClearLocationsOfT1Swaps();
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the MultipleCaBasedCellPopulation to visit.
     */
    virtual void Visit(MultipleCaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
    }

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
    {
    }
};

#endif /*VERTEXSWAPWRITERS_HPP_*/
