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

#ifndef ABSTRACTCELLWRITER_HPP_
#define ABSTRACTCELLWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedWriter.hpp"
#include "Cell.hpp"
#include "UblasVectorInclude.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCellPopulation;

/**
 * An abstract class for a writer that visits individual cells of a population and writes their data.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractCellWriter : public AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mOutputScalarData;
        archive & mOutputVectorData;
        archive & mVtkCellDataName;
        archive & mVtkVectorCellDataName;
    }

protected:

    /** Whether to output scalar data for VTK using GetCellDataForVtkOutput(). Default true. */
    bool mOutputScalarData;

    /** Whether to output scalar data for VTK using GetVectorCellDataForVtkOutput(). Default false. */
    bool mOutputVectorData;

    /** The name of the cell data used in VTK output. */
    std::string mVtkCellDataName;

    /** The name of the vector cell data used in VTK output. */
    std::string mVtkVectorCellDataName;

public:

    /**
     * Default constructor.
     * @param rFileName the name of the file to write to.
     */
    AbstractCellWriter(const std::string& rFileName);

    /**
     * Get a double associated with a cell. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * By default this method returns a DOUBLE_UNSET, but it may be overridden in subclasses
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell.
     *
     * @return data associated with the cell
     */
    virtual double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Get a c_vector associated with a cell. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * By default this method returns a c_vector of DOUBLE_UNSET, but it may be overridden in subclasses
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell.
     *
     * @return data associated with the cell
     */
    virtual c_vector<double, SPACE_DIM> GetVectorCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit a cell and write its data.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)=0;

    /**
     * Get whether to invoke GetCellDataForVtkOutput()
     *
     * @return mOutputScalarData
     */
    bool GetOutputScalarData();

    /**
     * Get whether to invoke GetVectorCellDataForVtkOutput()
     *
     * @return mOutputScalarData
     */
    bool GetOutputVectorData();

    /**
     * Set the name of the scalar cell data used in VTK output.
     * This method allows the user to change mVtkCellDataName from
     * its default value, which is set in each subclass's
     * constructor.
     *
     * @param vtkCellDataName the name of the VTK field
     */
    void SetVtkCellDataName(std::string vtkCellDataName);

    /**
     * Set the name of the vector cell data used in VTK output.
     * This method allows the user to change mVtkCellDataName from
     * its default value, which is set in each subclass's
     * constructor.
     *
     * @param vtkCellDataName the name of the VTK field
     */
    void SetVtkVectorCellDataName(std::string vtkCellDataName);

    /**
     * @return the name of the scalar cell data used in VTK output.
     */
    std::string GetVtkCellDataName();

    /**
     * @return the name of the vector cell data used in VTK output.
     */
    std::string GetVtkVectorCellDataName();
};

#endif /*ABSTRACTCELLWRITER_HPP_*/
