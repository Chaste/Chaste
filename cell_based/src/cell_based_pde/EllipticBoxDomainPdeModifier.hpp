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

#ifndef ELLIPTICBOXDOMAINPDEMODIFIER_HPP_
#define ELLIPTICBOXDOMAINPDEMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractBoxDomainPdeModifier.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "PetscTools.hpp"
#include "FileFinder.hpp"

/**
 * A modifier class in which an elliptic PDE is solved on a box domain and the results are stored in CellData.
 *
 * \todo Improve documentation (#2687)
 */
template<unsigned DIM>
class EllipticBoxDomainPdeModifier : public AbstractBoxDomainPdeModifier<DIM>
{
    friend class TestBoxDomainPdeModifiers;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractBoxDomainPdeModifier<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pPdeAndBcs shared pointer to a PDE object with associated boundary conditions (default to NULL)
     * @param pMeshCuboid pointer to a ChasteCuboid specifying the outer boundary for the FE mesh (defaults to NULL)
     * @param stepSize step size to be used in the FE mesh (defaults to 1.0, i.e. the default cell size)
     * @param solution solution vector (defaults to NULL)
     */
    EllipticBoxDomainPdeModifier(boost::shared_ptr<PdeAndBoundaryConditions<DIM> > pPdeAndBcs=boost::shared_ptr<PdeAndBoundaryConditions<DIM> >(),
                                 ChasteCuboid<DIM>* pMeshCuboid=NULL,
                                 double stepSize=1.0,
                                 Vec solution=NULL);

    /**
     * Destructor.
     */
    virtual ~EllipticBoxDomainPdeModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @param rCellPopulation reference to the cell population
     *
     * @return the full boundary conditions container
     */
    virtual std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipticBoxDomainPdeModifier)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const EllipticBoxDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    if (t->GetSolution())
    {
        std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
        PetscTools::DumpPetscObject(t->GetSolution(), archive_filename);
    }
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, EllipticBoxDomainPdeModifier<DIM> * t, const unsigned int file_version)
{
    Vec solution = NULL;

    std::string archive_filename = ArchiveLocationInfo::GetArchiveDirectory() + "solution.vec";
    FileFinder file_finder(archive_filename, RelativeTo::Absolute);

    if (file_finder.Exists())
    {
        PetscTools::ReadPetscObject(solution, archive_filename);
    }

    ::new(t)EllipticBoxDomainPdeModifier<DIM>(boost::shared_ptr<PdeAndBoundaryConditions<DIM> >(), NULL, 1.0, solution);
}
}
} // namespace ...

#endif /*ELLIPTICBOXDOMAINPDEMODIFIER_HPP_*/
