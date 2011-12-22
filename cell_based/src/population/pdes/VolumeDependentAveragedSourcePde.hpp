/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_
#define VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NodeBasedCellPopulation.hpp"
#include "AveragedSourcePde.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractLinearEllipticPde.hpp"

/**
 *  A PDE which calculates the source term by adding the number of cells
 *  in the element containing that point and scaling by the element area.
 */
template<unsigned DIM>
class VolumeDependentAveragedSourcePde : public AveragedSourcePde<DIM>
{
    friend class TestCellBasedPdes;

private:

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the PDE and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AveragedSourcePde<DIM> >(*this);
    }

    /** Static cast of the NodeBasedCellPopulation. **/
    NodeBasedCellPopulation<DIM>* mpStaticCastCellPopulation;

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param coefficient the coefficient of consumption of nutrient by cells (defaults to 0.0)
     */
    VolumeDependentAveragedSourcePde(AbstractCellPopulation<DIM>& rCellPopulation, double coefficient=0.0);

    /**
     * Set up the source terms.
     *
     * @param rCoarseMesh reference to the coarse mesh
     * @param pCellPdeElementMap optional pointer to the map from cells to coarse elements
     */
    void SetupSourceTerms(TetrahedralMesh<DIM,DIM>& rCoarseMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=NULL);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeDependentAveragedSourcePde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VolumeDependentAveragedSourcePde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VolumeDependentAveragedSourcePde<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a VolumeDependentAveragedSourcePde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VolumeDependentAveragedSourcePde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)VolumeDependentAveragedSourcePde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*VOLUMEDEPENDENTAVERAGEDSOURCEPDE_HPP_*/
