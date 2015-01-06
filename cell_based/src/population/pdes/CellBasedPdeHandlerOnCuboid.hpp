/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef CELLBASEDPDEHANDLERONCUBOID_HPP_
#define CELLBASEDPDEHANDLERONCUBOID_HPP_

#include <map>
#include <memory>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>

#include "CellBasedPdeHandler.hpp"

/**
 * A helper class, containing code for handling the numerical solution of one or more PDEs
 * (using the finite element method) associated with a cell-based simulation object.
 *
 * By letting AbstractCellBasedSimulation have a pointer to an object of this type as a
 * member variable, we separate out all PDE-related functionality into this class, and thus
 * obviate the need for specialized cell-based simulation subclasses.
 */
template<unsigned DIM>
class CellBasedPdeHandlerOnCuboid : public CellBasedPdeHandler<DIM>
{
    // Allow tests to access private members, in order to test computation of private functions
    friend class TestCellBasedPdeHandler;
    friend class TestOffLatticeSimulationWithPdes;
    friend class TestOnLatticeSimulationWithPdes;

private:

    /** Container for pointers to boundary conditions that are passed into the boundary condition containers. */
    std::vector<ConstBoundaryCondition<DIM>* > mConstBoundaryConditions;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<CellBasedPdeHandler<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to a cell population
     * @param deleteMemberPointersInDestructor whether to delete member pointers in the destructor (defaults to false)
     */
    CellBasedPdeHandlerOnCuboid(AbstractCellPopulation<DIM>* pCellPopulation, bool deleteMemberPointersInDestructor=false);

    /**
     * Destructor.
     */
    virtual ~CellBasedPdeHandlerOnCuboid();

    /**
     * Overridden ConstructBoundaryConditionsContainer method to implement
     * different boundary conditions on each face of the cuboid.
     *
     * @param pPdeAndBc a pointer to the PDE and BCs
     * @param pMesh the mesh on which to solve the PDE
     * @return The full boundary conditions container
     */
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ConstructBoundaryConditionsContainer(
            PdeAndBoundaryConditions<DIM>* pPdeAndBc,
            TetrahedralMesh<DIM,DIM>* pMesh);

    /**
     * Output parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellBasedPdeHandlerOnCuboid)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CellBasedPdeHandler.
 *
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CellBasedPdeHandlerOnCuboid<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = t->GetCellPopulation();
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CellBasedPdeHandler.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CellBasedPdeHandlerOnCuboid<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CellBasedPdeHandlerOnCuboid<DIM>(p_cell_population, true);
}
}
}

#endif /*CELLBASEDPDEHANDLERONCUBOID_HPP_*/
