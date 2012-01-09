/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef CRYPTSIMULATION2D_HPP_
#define CRYPTSIMULATION2D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "WntConcentration.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CryptSimulationBoundaryCondition.hpp"

/**
 * A 2D crypt simulation object. For more details on the crypt geometry, see the
 * papers by van Leeuwen et al (2009) [doi:10.1111/j.1365-2184.2009.00627.x] and
 * Osborne et al (2010) [doi:10.1098/rsta.2010.0173].
 */
class CryptSimulation2d : public OffLatticeSimulation<2>
{
    // Allow tests to access private members, in order to test computation of private functions e.g. DoCellBirth()
    friend class TestCryptSimulation2dWithMeshBasedCellPopulation;
    friend class TestCryptSimulation2dWithVertexBasedCellPopulation;

protected:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variable.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<2> >(*this);

        SerializableSingleton<WntConcentration<2> >* p_wnt_wrapper = WntConcentration<2>::Instance()->GetSerializationWrapper();
        archive & p_wnt_wrapper;

        archive & mWriteBetaCatenin;
    }

    /**
     * Whether the simulation includes the cell-cycle models
     * VanLeeuwen2009WntSwatCellCycleModelHypothesisOne or
     * VanLeeuwen2009WntSwatCellCycleModelHypothesisOne, and
     * hence whether beta catenin results are written to file.
     */
    bool mWriteBetaCatenin;

    /** The file that the values of beta catenin is written out to. */
    out_stream mVizBetaCateninResultsFile;

    /**
     * Helper member that stores whether we are using a MeshBasedCellPopulationWithGhostNodes
     * (if not, then we are assumed to be using a VertexBasedCellPopulation).
     */
    bool mUsingMeshBasedCellPopulation;

    /**
     * In the case of a MeshBasedCellPopulationWithGhostNodes, this method
     * calculates the new locations of a dividing cell's cell centres. The
     * node correspond to the dividing cell is moved a bit and the co-ordinates
     * of the new node are returned. This is done by drawing a random
     * direction (0->2PI) and placing the parent and daughter nodes in
     * opposing directions along this axis.
     *
     * In the case of a VertexBasedCellPopulation, by default this method
     * returns the zero vector. If the parent cell is a stem cell, then
     * this method returns the vector (0,1). This is then used by the
     * VertexBasedCellPopulation method AddCell() as the axis along which
     * the cell divides.
     *
     * @param pParentCell the parent cell
     *
     * @return daughter_coords the coordinates for the daughter cell.
     */
    c_vector<double, 2> CalculateCellDivisionVector(CellPtr pParentCell);

    /**
     * Use an output file handler to create a beta catenin results file.
     */
    void SetupWriteBetaCatenin();

    /**
     * Write beta catenin results to file.
     *
     * @param time the current time
     */
    virtual void WriteBetaCatenin(double time);

    /**
     * Overridden SetupSolve() method.
     *
     * Write initial beta catenin results to file if required.
     */
    void SetupSolve();

    /**
     * Overridden PostSolve() method.
     *
     * Write current beta catenin results to file if required.
     */
    void PostSolve();

    /**
     * Overridden AfterSolve() method.
     *
     * Closes beta catenin results file if required, then calls
     * the base class method.
     */
    void AfterSolve();

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to false when loading from an archive)
     */
    CryptSimulation2d(AbstractCellPopulation<2>& rCellPopulation,
                      bool deleteCellPopulationInDestructor=false,
                      bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the CryptSimulationBoundaryCondition.
     */
    virtual ~CryptSimulation2d();

    /**
     * Set method for mUseJiggledBottomCells.
     */
    void UseJiggledBottomCells();

    /**
     * Sets the Ancestor index of all the cells at the bottom in order,
     * can be used to trace clonal populations.
     */
    void SetBottomCellAncestors();

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation2d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation2d * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation2d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation2d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar & p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSimulation2d(*p_cell_population, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION2D_HPP_*/
