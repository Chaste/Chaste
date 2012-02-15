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

#ifndef NODEBASEDCELLPOPULATIONWITHPARTICLES_HPP_
#define NODEBASEDCELLPOPULATIONWITHPARTICLES_HPP_

#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A NodeBasedCellPopulationWithParticles is a NodeBasedCellPopulation with cells interacting with
 * particles, such as ECM or fluid particles.
 */
template<unsigned DIM>
class NodeBasedCellPopulationWithParticles : public NodeBasedCellPopulation<DIM>
{
    friend class TestNodeBasedCellPopulationWithParticles;

protected:
    /**
     * Check consistency of our internal data structures.
     */
    void Validate();

private:

    /** Records whether a node is a particle or not */
    std::vector<bool> mIsParticle;

    /**
     * Set the particles by taking in a set of which nodes indices are particles.
     *
     * @param rParticleIndices set of node indices corresponding to particles
     */
    void SetParticles(const std::set<unsigned>& rParticleIndices);

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
        archive & mIsParticle;
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
    NodeBasedCellPopulationWithParticles(NodesOnlyMesh<DIM>& rMesh,
                            std::vector<CellPtr>& rCells,
                            const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                            bool deleteMesh=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a mutable nodes-only mesh
     */
    NodeBasedCellPopulationWithParticles(NodesOnlyMesh<DIM>& rMesh);

	/**
	 * Update particle positions.
	 *
	 * @param rNodeForces
	 * @param dt
	 */
	void UpdateParticlePositions(const std::vector< c_vector<double, DIM> >& rNodeForces,
					  			 double dt);

    /**
     * Update mIsParticle if required by a remesh.
     *
     * @param rMap A map between node indices before and after remesh
     */
    void UpdateParticlesAfterReMesh(NodeMap& rMap);

    /**
     * Overridden UpdateNodeLocation() method.
     *
     * Update the location of each node in the cell population given
     * a two vectors of forces on cells and particles and a time step over which
     * to integrate the equations of motion.
     *
     * @param rNodeForcesOnCells  forces on cells
     * @param rNodeForcesOnParticles  forces on particles
     * @param dt  time step
     */
    void UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt);

    /**
	 * @return mIsParticle, a vector that stores which nodes are particles.
	 */
	std::vector<bool>& rGetParticles();

	/**
	 * IsParticle() method.
	 *
	 * Find if a given node is a particle
	 *
	 * @param index the global index of a specified node
	 *
	 * @return whether the node is a particle
	 */
	bool IsParticle(unsigned index);

	/**
	 * @return the indices of those nodes that are particles.
	 */
	std::set<unsigned> GetParticleIndices();


	/**
	 * Overridden AddCell() method.
	 *
	 * Add a new cell to the cell population and update mIsParticle.
	 *
	 * @param pNewCell  the cell to add
	 * @param rCellDivisionVector  the position in space at which to put it
	 * @param pParentCell pointer to a parent cell  - this is required for
	 *  mesh-based cell populations
	 * @return address of cell as it appears in the cell list (internal of this method uses a copy constructor along the way)
	 */
	CellPtr AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell);

	/**
	 * Overridden WriteVtkResultsToFile method.
	 */
	void WriteVtkResultsToFile();

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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithParticles)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a NodeBasedCellPopulationWithParticles.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const NodeBasedCellPopulationWithParticles<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const NodesOnlyMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a NodeBasedCellPopulationWithParticles.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, NodeBasedCellPopulationWithParticles<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    NodesOnlyMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)NodeBasedCellPopulationWithParticles<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*NODEBASEDCELLPOPULATIONWITHPARTICLES_HPP_*/
