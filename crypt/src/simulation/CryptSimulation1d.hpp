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

#ifndef CRYPTSIMULATION1D_HPP_
#define CRYPTSIMULATION1D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "WntConcentration.hpp"
#include "CryptSimulationBoundaryCondition.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"

/**
 * A 1D crypt simulation object. The model is a simplified version of a 2D crypt model
 * developed by Meineke et al (doi:10.1046/j.0960-7722.2001.00216.x).
 */
class CryptSimulation1d : public OffLatticeSimulation<1>
{
    // Allow tests to access private members to test private functions
    friend class TestCryptSimulation1d;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<1> >(*this);

        SerializableSingleton<WntConcentration<1> >* p_wnt_wrapper = WntConcentration<1>::Instance()->GetSerializationWrapper();
        archive & p_wnt_wrapper;
    }

    /** Helper member that is a static cast of the cell population. */
    MeshBasedCellPopulation<1>* mpStaticCastCellPopulation;

    /**
     * Calculates the new locations of a dividing cell's cell centres.
     * Moves the dividing node a bit and returns co-ordinates for the new node.
     * It does this by picking a random direction (0->2PI) and placing the parent
     * and daughter in opposing directions on this axis.
     *
     * @param pParentCell the parent cell
     *
     * @return daughter_coords the coordinates for the daughter cell.
     */
    c_vector<double, 1> CalculateCellDivisionVector(CellPtr pParentCell);

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to false when loading from an archive)
     */
    CryptSimulation1d(AbstractCellPopulation<1>& rCellPopulation,
                      bool deleteCellPopulationInDestructor=false,
                      bool initialiseCells=true);

    /**
     * Destructor.
     *
     * This frees the CryptSimulationBoundaryCondition.
     */
    virtual ~CryptSimulation1d();

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
CHASTE_CLASS_EXPORT(CryptSimulation1d)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CryptSimulation1d.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CryptSimulation1d * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<1>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CryptSimulation1d.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CryptSimulation1d * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<1>* p_cell_population;
    ar & p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)CryptSimulation1d(*p_cell_population, true, false);
}
}
} // namespace

#endif /*CRYPTSIMULATION1D_HPP_*/
