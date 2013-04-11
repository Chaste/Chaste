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

#ifndef CELLWRITERS_HPP_
#define CELLWRITERS_HPP_

#include "AbstractCellWriter.hpp"
#include "AbstractOdeBasedCellCycleModel.hpp"

/** A class written using the visitor pattern for writing cell proliferative types to file */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellProliferativeTypesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellProliferativeTypesWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
		this->mFileName = "results.vizcelltypes";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
    	*this->mpOutStream << pCell->GetCellProliferativeType()->GetColour() << " ";
    }

};

/** A class written using the visitor pattern for writing cell ages to file. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellAgesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellAgesWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
		this->mFileName = "cellages.dat";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        // Write location index corresponding to cell
        *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

        // Write cell location
        c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);

        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << cell_location[i] << " ";
        }

        // Write cell age
        *this->mpOutStream << pCell->GetAge() << " ";
    }
};

/** A class written using the visitor pattern for writing cell ancestors */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellAncestorWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellAncestorWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
			this->mFileName = "results.vizancestors";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        unsigned ancestor = pCell->GetAncestor();
        if (ancestor == UNSIGNED_UNSET)
        {
            // Set the file to -1 to mark this case.
    		ancestor = 1;
            *this->mpOutStream << "-";
        }
        *this->mpOutStream << ancestor << " ";
    }
};

/** A class written using the visitor pattern for writing cell ids from a cell population to file. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellIdWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellIdWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
			this->mFileName = "loggedcell.dat";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        unsigned cell_id = pCell->GetCellId();
        unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
        *this->mpOutStream << " " << cell_id << " " << location_index;

        c_vector<double, SPACE_DIM> coords = pCellPopulation->GetLocationOfCellCentre(pCell);
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << " " << coords[i];
        }
    }
};

/** A class written using the visitor pattern for writing cell locations to file. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellLocationWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellLocationWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
		this->mFileName = "results.vizlocations";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
    	unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);

    	*this->mpOutStream << location_index << " ";
    }
};

/** A class written using the visitor pattern for writing the cell proliferative phase to file. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellProliferativePhasesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellProliferativePhasesWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
		this->mFileName = "results.vizcellphases";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
    	*this->mpOutStream << pCell->GetCellCycleModel()->GetCurrentCellCyclePhase() << " ";
    }
};

/** A class written using the visitor pattern for writing cell variables to file. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellVariablesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellVariablesWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
		this->mFileName = "cellvariables.dat";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
    	if (dynamic_cast<AbstractOdeBasedCellCycleModel*>(pCell->GetCellCycleModel()))
    	{
    		// Write location index corresponding to cell
    		*this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

    		// Write cell variables
    		std::vector<double> proteins = (static_cast<AbstractOdeBasedCellCycleModel*>(pCell->GetCellCycleModel()))->GetProteinConcentrations();
    		for (unsigned i=0; i<proteins.size(); i++)
    		{
    				*this->mpOutStream << proteins[i] << " ";
    		}
    	}
    }
};

/** A class written using the visitor pattern for writing cell volumes to file. */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellVolumesWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Default constructor
     * @param directory the path to the directory in to which this class should write.
     */
	CellVolumesWriter(std::string directory)
		: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(directory)
	{
		this->mFileName = "cellareas.dat";
	}

    /**
     * Visit a cell and write its data.
     *
     * @param pCell the cell to write
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
    	unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    	unsigned cell_id = pCell->GetCellId();
    	c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    	double volume = pCellPopulation->GetVolumeOfCell(pCell);

    	*this->mpOutStream << location_index << " " << cell_id << " ";
    	for (unsigned i=0; i<SPACE_DIM; i++)
    	{
    			*this->mpOutStream << centre_location[i] << " ";
    	}
    	*this->mpOutStream << volume << " ";
    }
};

#endif /*CELLWRITERS_HPP_*/
