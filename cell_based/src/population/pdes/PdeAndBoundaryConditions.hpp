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

#ifndef PDEANDBOUNDARYCONDITIONS_HPP_
#define PDEANDBOUNDARYCONDITIONS_HPP_

#include "ChasteSerialization.hpp"
#include "AbstractBoundaryCondition.hpp"
#include "ArchiveLocationInfo.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "AveragedSourcePde.hpp"
#include "PetscTools.hpp"
#include "FileFinder.hpp"

/**
 * A helper class for use in cell-based simulations with PDEs. The class
 * contains a pointer to a linear elliptic PDE, which is to be solved
 * on the domain defined by the cell population. The class also contains
 * information describing the boundary condition that is to be imposed
 * when solving the PDE. Currently we allow Neumann (imposed flux) or
 * Dirichlet (imposed value) boundary conditions. The boundary condition
 * may be constant on the boundary or vary spatially and/or temporally.
 * In cell-based simulations with PDEs, one or more of these objects are
 * accessed via the CellBasedHandler class.
 */
template<unsigned DIM>
class PdeAndBoundaryConditions
{
    friend class TestPdeAndBoundaryConditions;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the PDE object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Note that archiving of mSolution is handled by the methods save/load_construct_data
        archive & mpPde;
        archive & mpBoundaryCondition;
        archive & mIsNeumannBoundaryCondition;
    }

    /** Pointer to a linear elliptic PDE object. */
    AbstractLinearEllipticPde<DIM,DIM>* mpPde;

    /** Pointer to a boundary condition object. */
    AbstractBoundaryCondition<DIM>* mpBoundaryCondition;

    /** Whether the boundary condition is Neumann (false corresponds to a Dirichlet boundary condition). */
    bool mIsNeumannBoundaryCondition;

    /** The solution to the PDE problem, for use as an initial guess when solving at the next time step. */
    Vec mSolution;

    /** Whether to delete member pointers in the destructor (used in archiving). */
    bool mDeleteMemberPointersInDestructor;

public:

    /**
     * Constructor.
     *
     * @param pPde A pointer to a linear elliptic PDE object (defaults to NULL)
     * @param pBoundaryCondition A pointer to an abstract boundary condition
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param solution A solution vector (defaults to NULL)
     * @param deleteMemberPointersInDestructor whether to delete member pointers in the destructor
     *     (defaults to false)
     */
    PdeAndBoundaryConditions(AbstractLinearEllipticPde<DIM,DIM>* pPde=NULL,
                             AbstractBoundaryCondition<DIM>* pBoundaryCondition=NULL,
                             bool isNeumannBoundaryCondition=true,
                             Vec solution=NULL,
                             bool deleteMemberPointersInDestructor=false);

    /**
     * Destructor.
     */
    ~PdeAndBoundaryConditions();

    /**
     * @return mpPde
     */
    AbstractLinearEllipticPde<DIM,DIM>* GetPde();

    /**
     * @return mpBoundaryCondition
     */
    AbstractBoundaryCondition<DIM>* GetBoundaryCondition() const;

    /**
     * @return mSolution
     */
    Vec GetSolution();

    /**
     * @return mSolution (used in archiving)
     */
    Vec GetSolution() const;

    /**
     * Set mSolution.
     *
     * @param solution the current solution
     */
    void SetSolution(Vec solution);

    /**
     * @return mIsNeumannBoundaryCondition
     */
    bool IsNeumannBoundaryCondition();

    /**
     * @return whether the PDE is of type AveragedSourcePde
     */
    bool HasAveragedSourcePde();

    /**
     * Call VecDestroy on mSolution.
     */
    void DestroySolution();

    /**
     * In the case where mpPde is of type AveragedSourcePde, set the source terms
     * using the information in the given mesh.
     *
     * @param pMesh Pointer to a tetrahedral mesh
     * @param pCellPdeElementMap map between cells and elements, from CellBasedPdeHandler
     */
    void SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM,DIM>* pMesh, std::map< CellPtr, unsigned >* pCellPdeElementMap=NULL);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PdeAndBoundaryConditions)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const PdeAndBoundaryConditions<DIM> * t, const unsigned int file_version)
{
    if (t->GetSolution())
    {
        std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
        PetscTools::DumpPetscObject(t->GetSolution(), archive_filename);
    }
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, PdeAndBoundaryConditions<DIM> * t, const unsigned int file_version)
{
    Vec solution = NULL;

    std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
    FileFinder file_finder(archive_filename, RelativeTo::Absolute);

    if (file_finder.Exists())
    {
        PetscTools::ReadPetscObject(solution, archive_filename);
    }

    ::new(t)PdeAndBoundaryConditions<DIM>(NULL, NULL, false, solution, true);
}
}
} // namespace ...

#endif /* PDEANDBOUNDARYCONDITIONS_HPP_ */
