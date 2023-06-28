/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef SEMBASEDCELLPOPULATION_HPP_
#define SEMBASEDCELLPOPULATION_HPP_

#include "AbstractOffLatticeCellPopulation.hpp"
#include "SemMesh.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A facade class encapsulating a subcellular element method (SEM) based 
 * cell population.
 *
 * Contains a group of cells and maintains the associations between CellPtrs and 
 * SemElements in the SemMesh.
 */
template<unsigned DIM>
class SemBasedCellPopulation : public AbstractOffLatticeCellPopulation<DIM>
{
private:

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this cell population has been de-serialized.
     */
    bool mDeleteMesh;

    /**
     * A static cast of the AbstractMesh from AbstractCellPopulation for use in 
     * this class.
     */
    SemMesh<DIM>* mpSemMesh;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by 
     * load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOffLatticeCellPopulation<DIM> >(*this);
    }

    /**
     * Check the consistency of internal data structures.
     * Each SemElement must have a CellPtr associated with it.
     */
    void Validate() override;

public:

    /**
     * Create a new cell population facade from a SemMesh and collection of 
     * CellPtrs.
     *
     * There must be precisely one CellPtr for each SemElement in the SemMesh.
     *
     * @param rMesh reference to a SemMesh
     * @param rCells reference to a vector of CellPtrs
     * @param deleteMesh set to true if you want the cell population to free the 
     *                   mesh memory on destruction (defaults to false)
     * @param validate whether to validate the cell population when it is 
     *                 created (defaults to true)
     * @param locationIndices an optional vector of location indices that 
     *                        correspond to real cells
     */
    SemBasedCellPopulation(SemMesh<DIM>& rMesh,
                           std::vector<CellPtr>& rCells,
                           bool deleteMesh=false,
                           bool validate=true,
                           const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Constructor for use by boost serialization only.
     *
     * @param rMesh a SemMesh
     */
    SemBasedCellPopulation(SemMesh<DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~SemBasedCellPopulation();

    /**
     * @return reference to mpSemMesh.
     */
    SemMesh<DIM>& rGetMesh();

    /**
     * @return const reference to mpSemMesh (used in archiving).
     */
    const SemMesh<DIM>& rGetMesh() const;

    /**
     * Get a particular SemElement.
     *
     * @param elementIndex the global index of the SemElement
     *
     * @return a pointer to the SemElement.
     */
    SemElement<DIM>* GetElement(unsigned elementIndex);

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes() override;

    /**
     * Overridden GetLocationOfCellCentre() method.
     *
     * Find the centre of mass of a given cell (assuming uniform density).
     * Note that, as there is no guarantee of convexity, this may lie
     * outside the SemElement corresponding to the cell.
     *
     * @param pCell pointer to a cell in the population
     *
     * @return the location of the centre of mass of the SemElement 
     *         corresponding to this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell) override;

    /**
     * Overridden GetNode() method.
     *
     * @param index global index of the specified node
     *
     * @return a pointer to the node.
     */
    Node<DIM>* GetNode(unsigned index) override;

    /**
     * Overridden GetNeighbouringLocationIndices() method.
     *
     * Given a cell, returns the set of location indices corresponding to 
     * neighbouring cells.
     *
     * @param pCell pointer to a cell in the population
     * @return the set of neighbouring location indices.
     */
    std::set<unsigned> GetNeighbouringLocationIndices(CellPtr pCell) override;

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the cell population.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population
     */
    unsigned AddNode(Node<DIM>* pNewNode) override;

    /**
     * Overridden SetNode() method.
     *
     * Move the node with a given index to a new point in space.
     *
     * @param index the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned index, ChastePoint<DIM>& rNewLocation) override;

    /**
     * Get a pointer to the SemElement corresponding to a given CellPtr.
     *
     * @param pCell pointer to a cell in the population
     *
     * @return pointer to the SemElement.
     */
    SemElement<DIM>* GetElementCorrespondingToCell(CellPtr pCell);

    /**
     * Overridden GetVolumeOfCell() method.
     *
     * @param pCell pointer to a cell in the population
     * 
     * @return volume via associated SemElement.
     */
    double GetVolumeOfCell(CellPtr pCell) override;

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SemBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SemBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const SemMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a VertexBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SemBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    SemMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)SemBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...

#endif /*SEMBASEDCELLPOPULATION_HPP_*/