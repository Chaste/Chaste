/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef ELLIPTICPDEANDBOUNDARYCONDITIONS_HPP_
#define ELLIPTICPDEANDBOUNDARYCONDITIONS_HPP_

#include "AbstractPdeAndBoundaryConditions.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearEllipticPde.hpp"

/**
 * A helper class for use in cell-based simulations with PDEs. The class
 * contains a pointer to a linear elliptic PDE, which is to be solved
 * on the domain defined by the cell population. The class also contains
 * information describing the boundary condition that is to be imposed
 * when solving the PDE. Currently we allow Neumann (imposed flux) or
 * Dirichlet (imposed value) boundary conditions. The boundary condition
 * may be constant on the boundary or vary spatially and/or temporally.
 * In cell-based simulations with PDEs, one or more of these objects are
 * accessed via the CellBasedPdeHandler class.
 *
 * \todo Rename as EllipticPdeAndBoundaryConditions and move to cell_based/src/cell_based/pde (#2687)
 *       and merge, or create shared parent class, with ParabolicPdeAndBoundaryConditions (#2767)
 */
template<unsigned DIM>
class EllipticPdeAndBoundaryConditions : public AbstractPdeAndBoundaryConditions
{
    friend class TestEllipticPdeAndBoundaryConditions;

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
        archive & boost::serialization::base_object<AbstractPdeAndBoundaryConditions<DIM> >(*this);
        archive & mpPde;
    }

    /** Pointer to a linear elliptic PDE object. */
    AbstractLinearEllipticPde<DIM,DIM>* mpPde;

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
     * @return whether the PDE is of type AveragedSourceEllipticPde
     */
    bool HasAveragedSourcePde();

    /**
     * In the case where mpPde is of type AveragedSourceEllipticPde, set the source terms
     * using the information in the given mesh.
     *
     * @param pMesh Pointer to a tetrahedral mesh
     * @param pCellPdeElementMap map between cells and elements, from CellBasedPdeHandler
     */
    void SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM,DIM>* pMesh, std::map< CellPtr, unsigned >* pCellPdeElementMap=NULL);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipticPdeAndBoundaryConditions)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const EllipticPdeAndBoundaryConditions<DIM> * t, const unsigned int file_version)
{
    if (t->GetSolution())
    {
        std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
        PetscTools::DumpPetscObject(t->GetSolution(), archive_filename);
    }
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, EllipticPdeAndBoundaryConditions<DIM> * t, const unsigned int file_version)
{
    Vec solution = NULL;

    std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
    FileFinder file_finder(archive_filename, RelativeTo::Absolute);

    if (file_finder.Exists())
    {
        PetscTools::ReadPetscObject(solution, archive_filename);
    }

    ::new(t)EllipticPdeAndBoundaryConditions<DIM>(NULL, NULL, false, solution, true);
}
}
} // namespace ...

#endif /* ELLIPTICPDEANDBOUNDARYCONDITIONS_HPP_ */
