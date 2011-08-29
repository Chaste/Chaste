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

#ifndef NODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_
#define NODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_

#include "NodeBasedCellPopulation.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


/**
 * A NodeBasedCellPopulationWithBuskeUpdate is a CellPopulation consisting of only nodes in space with associated cells.
 * There are no elements and no mesh. The cell positions are updated according to Buske et al (2011) (doi:10.1371/journal.pcbi.1001045).
 */
template<unsigned DIM>
class NodeBasedCellPopulationWithBuskeUpdate : public NodeBasedCellPopulation<DIM>
{
    friend class TestNodeBasedCellPopulationWithBuskeUpdate;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the nodes is handled by load/save_construct_data,
     * so we don't actually have to do anything here except delegate to the base class.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<NodeBasedCellPopulation<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * Note that the cell population will take responsibility for freeing the memory used by the nodes.
     *
     * @param rMesh a mutable nodes-only mesh
     * @param rCells a vector of cells
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param deleteMesh whether to delete nodes-only mesh in destructor
     */
    NodeBasedCellPopulationWithBuskeUpdate(NodesOnlyMesh<DIM>& rMesh,
                                           std::vector<CellPtr>& rCells,
                                           const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                           bool deleteMesh=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable nodes-only mesh
     */
    NodeBasedCellPopulationWithBuskeUpdate(NodesOnlyMesh<DIM>& rMesh);

    /**
     * Overridden UpdateNodeLocations() method. To use the update from Buske et al (2011) (doi:10.1371/journal.pcbi.1001045).
     *
     * @param rNodeForces a vector containing the force on each node in the cell population
     * @param dt the time step
     */
    virtual void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt);

    /**
     * Outputs CellPopulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithBuskeUpdate)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a NodeBasedCellPopulationWithBuskeUpdate.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const NodeBasedCellPopulationWithBuskeUpdate<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const NodesOnlyMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a NodeBasedCellPopulationWithBuskeUpdate.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, NodeBasedCellPopulationWithBuskeUpdate<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    NodesOnlyMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)NodeBasedCellPopulationWithBuskeUpdate<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*NODEBASEDCELLPOPULATIONWITHBUSKEUPDATE_HPP_*/
