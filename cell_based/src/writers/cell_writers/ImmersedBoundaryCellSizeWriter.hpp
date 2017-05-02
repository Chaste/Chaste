//
// Created by bartmanski on 02/05/17.
//

#ifndef CHASTE_IMMERSEDBOUNDARYCELLSIZEWRITER_HPP
#define CHASTE_IMMERSEDBOUNDARYCELLSIZEWRITER_HPP

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A class written using the visitor pattern for writing the size of each cell in the cell population.
 *
 * The output file is called ib_cell_size.dat by default. If VTK is switched on, then the writer also specifies the
 * VTK output for each cell, which is stored in the VTK cell data "Cell Size" by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryCellSizeWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
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
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    ImmersedBoundaryCellSizeWriter();

    /**
     * Overridden GetCellDataForVtkOutput() method.
     *
     * Get an unsigned integer associated with a cell. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return data associated with the cell
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Overridden VisitCell() method.
     *
     * Visit a cell and write its size.
     *
     * Outputs a line of space-separated values of the form:
     * ...[location index] [cell id] [x-pos] [y-pos] [z-pos] [cell size] ...
     * with [y-pos] and [z-pos] included for 2 and 3 dimensional simulations, respectively.
     *
     * This is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ImmersedBoundaryCellSizeWriter)

#endif //CHASTE_IMMERSEDBOUNDARYCELLSIZEWRITER_HPP
