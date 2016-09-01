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

#ifndef PDEANDBOUNDARYCONDITIONS_HPP_
#define PDEANDBOUNDARYCONDITIONS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearPde.hpp"
#include "AbstractBoundaryCondition.hpp"
#include "ArchiveLocationInfo.hpp"
#include "TetrahedralMesh.hpp"
#include "Cell.hpp"
#include "PetscTools.hpp"
#include "FileFinder.hpp"

/**
 * Abstract class storing shared functionality of PdeAndBoundaryConditions 
 * and PdeAndBoundaryConditions.
 */
template<unsigned DIM>
class PdeAndBoundaryConditions
{
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
        // Note that archiving of mpPde and mSolution is handled by the methods save/load_construct_data
        archive & mpPde;
        archive & mpBoundaryCondition;
        archive & mIsNeumannBoundaryCondition;
        archive & mDependentVariableName;
    }

protected:

    /** Pointer to a linear elliptic PDE object. */
    AbstractLinearPde<DIM,DIM>* mpPde;

    /** Pointer to a boundary condition object. */
    AbstractBoundaryCondition<DIM>* mpBoundaryCondition;

    /** Whether the boundary condition is Neumann (false corresponds to a Dirichlet boundary condition). */
    bool mIsNeumannBoundaryCondition;

    /**
     * The solution to the PDE problem, for use as an initial guess when solving at the next time step.
     *
     * \todo Consider removing this member and associated methods (#2687)
     */
    Vec mSolution;

    /** Whether to delete member pointers in the destructor (used in archiving). */
    bool mDeleteMemberPointersInDestructor;

    /**
     * For use in PDEs where we know what the quantity for which we are solving is called,
     * e.g. oxygen concentration.
     */
    std::string mDependentVariableName;

public:

    /**
     * Constructor.
     *
     * @param pPde A pointer to a linear PDE object (defaults to NULL)
     * @param pBoundaryCondition A pointer to an abstract boundary condition
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param solution A solution vector (defaults to NULL)
     * @param deleteMemberPointersInDestructor whether to delete member pointers in the destructor
     *     (defaults to false)
     */
    PdeAndBoundaryConditions(AbstractLinearPde<DIM,DIM>* pPde=NULL,
                                     AbstractBoundaryCondition<DIM>* pBoundaryCondition=NULL,
                                     bool isNeumannBoundaryCondition=true,
                                     Vec solution=NULL,
                                     bool deleteMemberPointersInDestructor=false);

    /**
     * Destructor.
     */
    ~PdeAndBoundaryConditions();

    /**
     * @return mpPde (used in archiving)
     */
    const AbstractLinearPde<DIM,DIM>* GetConstPde() const;

    /**
     * @return mpPde (used in archiving)
     */
    AbstractLinearPde<DIM,DIM>* GetPde() const;

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
     * Call PetscTools::Destroy on mSolution.
     */
    void DestroySolution();

    /**
     * Set the name of the dependent variable.
     *
     * @param rName the name.
     */
    void SetDependentVariableName(const std::string& rName);

    /**
     * Get the name of the dependent variable.
     *
     * @return the name
     */
    std::string& rGetDependentVariableName();

    /**
     * @return whether the PDE has an averaged source
     */
    bool HasAveragedSourcePde();

    /**
     * In the case where the PDE has an averaged source, set the source terms
     * using the information in the given mesh.
     *
     * @param pMesh Pointer to a tetrahedral mesh
     * @param pCellPdeElementMap map between cells and elements
     */
    void SetUpSourceTermsForAveragedSourcePde(TetrahedralMesh<DIM,DIM>* pMesh, std::map<CellPtr, unsigned>* pCellPdeElementMap=NULL);
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
