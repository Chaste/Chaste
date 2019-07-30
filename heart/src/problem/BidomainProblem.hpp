/*

Copyright (c) 2005-2019, University of Oxford.
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


#ifndef BIDOMAINPROBLEM_HPP_
#define BIDOMAINPROBLEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCardiacProblem.hpp"
#include "AbstractCardiacTissue.hpp"
#include "AbstractBidomainSolver.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "Electrodes.hpp"
#include "BidomainTissue.hpp"
#include "HeartRegionCodes.hpp"
#include "DistributedTetrahedralMesh.hpp"

/**
 * Class which specifies and solves a bidomain problem.
 *
 * The solution vector is of the form:
 * (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N),
 * where V_j is the voltage at node j and phi_j is the
 * extracellular potential at node j.
 */
template<unsigned DIM>
class BidomainProblem : public AbstractCardiacProblem<DIM,DIM, 2>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;

    // #1082
    friend class TestPCTwoLevelsBlockDiagonal;

    /**
     * Save the member variables to an archive.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractCardiacProblem<DIM, DIM, 2> >(*this);
        archive & mpBidomainTissue;
        //archive & mExtracelluarColumnId; // Created by InitialiseWriter, called from Solve
        archive & mRowForAverageOfPhiZeroed;
        archive & mHasBath;
        archive & mpElectrodes;
    }

    /**
     * Load the member variables from an archive.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacProblem<DIM, DIM, 2> >(*this);
        archive & mpBidomainTissue;
        //archive & mExtracelluarColumnId; // Created by InitialiseWriter, called from Solve
        archive & mRowForAverageOfPhiZeroed;
        archive & mHasBath;
        archive & mpElectrodes;

        if (mHasBath)
        {
            // We only save element annotations, so annotate bath nodes from these
            AnalyseMeshForBath();
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    friend class TestBidomainWithBathProblem;
    friend class TestCardiacSimulationArchiver;

protected:
    /** The bidomain PDE */
    BidomainTissue<DIM>* mpBidomainTissue;

    /** Nodes at which the extracellular voltage is fixed to zero (replicated) */
    std::vector<unsigned> mFixedExtracellularPotentialNodes;
   /** Used by the writer */
    unsigned mExtracelluarColumnId;
    /**
     * Another method of resolving the singularity in the bidomain equations.
     * Specifies a row of the matrix at which to impose an extra condition.
     */
    unsigned mRowForAverageOfPhiZeroed;

    /** Whether the mesh has a bath, ie whether this is a bath simulation */
    bool mHasBath;

    /** Electrodes used to provide a shock */
    boost::shared_ptr<Electrodes<DIM> > mpElectrodes;

    /**
     *  Create normal initial condition but overwrite V to zero for bath nodes, if
     *  there are any.
     *  @return the newly created intial conditions vector
     */
    Vec CreateInitialCondition();

    /**
     * Annotate bath nodes with the correct region code, if a bath is present.
     * Will throw if #mHasBath is set but no bath is present in the mesh.
     */
    void AnalyseMeshForBath();

    /**
     * We need to save the solver that is being used to be able to switch
     * off the electrodes (by adding default boundary conditions to the
     * solver)
     */
    AbstractBidomainSolver<DIM,DIM>* mpSolver;

    /** @return a newly created bidomain PDE object */
    virtual AbstractCardiacTissue<DIM> *CreateCardiacTissue();

    /** @return a newly created suitable bidomain solver */
    virtual AbstractDynamicLinearPdeSolver<DIM,DIM,2>* CreateSolver();

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     *   create cells.
     * @param hasBath Whether the simulation has a bath (if this is true, all elements with
     *   attribute = 1 will be set to be bath elements (the rest should have
     *   attribute = 0)).
     *
     */
    BidomainProblem(AbstractCardiacCellFactory<DIM>* pCellFactory, bool hasBath=false);

    /**
     * Constructor just used for archiving
     */
    BidomainProblem();

    /**
     *  Set the nodes at which phi_e (the extracellular potential) is fixed to
     *  zero. This does not necessarily have to be called. If it is not, phi_e
     *  is only defined up to a constant.
     *
     *  @param nodes  the nodes to be fixed.
     *
     *  @note currently, the value of phi_e at the fixed nodes cannot be set to be
     *  anything other than zero.
     */
    void SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes);

    /**
     * Set which row of the linear system should be used to enforce the
     * condition that the average of phi_e is zero.  If not called, this
     * condition will not be used.
     *
     * @param node  the mesh node index giving the row at which to impose the constraint
     */
    void SetNodeForAverageOfPhiZeroed(unsigned node);

    /**
     *  @return the pde. Can only be called after Initialise()
     */
    BidomainTissue<DIM>* GetBidomainTissue();

    /**
     *  Print out time and max/min voltage/phi_e values at current time.
     * @param time  current time.
     */
    void WriteInfo(double time);

    /**
     * Define what variables are written to the primary results file.
     * Adds the extracellular potential.
     * @param extending  whether we are extending an existing results file
     */
    virtual void DefineWriterColumns(bool extending);

    /**
     * Write one timestep of output data to the primary results file.
     * Adds the extracellular potential to the results.
     *
     * @param time  the current time
     * @param voltageVec  the solution vector to write
     */
    virtual void WriteOneStep(double time, Vec voltageVec);

    /**
     * Performs a series of checks before solving.
     * Checks that a suitable method of resolving the singularity is being used.
     */
    void PreSolveChecks();

    /**
     *  Called at beginning of each time step in the main time-loop in
     *  AbstractCardiacProblem::Solve(). Overloaded here to switch on
     *  the electrodes (if there are any).
     *
     * @param time  the current time
     */
    void AtBeginningOfTimestep(double time);

    /**
     *  Called at end of each time step in the main time-loop in
     *  AbstractCardiacProblem::Solve(). Overloaded here to switch off
     *  the electrodes (if there are any).
     *
     * @param time  the current time
     */
    void OnEndOfTimestep(double time);

    /**
     * Method to fill in a vector of additional stopping times consisting of the times when
     * the electrodes are to be turned on and off
     *
     * @param rAdditionalStoppingTimes reference to vector that will contain the on and off times
     *  for the electrodes
     */
    void SetUpAdditionalStoppingTimes(std::vector<double>& rAdditionalStoppingTimes);

    /**
     * Used when loading a set of archives written by a parallel simulation onto a single process.
     * Loads data from the given process-specific archive (written by a non-master process) and
     * merges it into our data.
     *
     * @param archive  the archive to load
     * @param version  the archive file version
     *
     * \note The process-specific archives currently contain the following data.  If the layout changes,
     * then this method will need to be altered, since it hard-codes knowledge of the order in
     * which things are archived.
     *
     *  -# Stuff known by AbstractCardiacProblem
     *  -# #mpElectrodes->mpBoundaryConditionsContainer
     *
     * This gets called by AbstractCardiacProblem::LoadExtraArchive when it's done the generic stuff.
     */
    template<class Archive>
    void LoadExtraArchiveForBidomain(Archive & archive, unsigned version);

    /**
     * @return whether this is a bidomain problem with bath or not
     */
    bool GetHasBath();

    /**
     *  Set an electrode object (which provides boundary conditions). Only
     *  valid if there is a bath.
     */
    void SetElectrodes();
};

#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainProblem)


template<unsigned DIM>
template<class Archive>
void BidomainProblem<DIM>::LoadExtraArchiveForBidomain(Archive & archive, unsigned version)
{
    // Not all bidomain problems have electrodes...
    if (mpElectrodes)
    {
        // Electrodes will always have a BCC object by this point
        assert(mpElectrodes->GetBoundaryConditionsContainer());
        // We might need to get some of the boundary conditions from this archive, but it might be
        // the case that the problem's BCC is the same object as the Electrodes' (if they are turned
        // on) in which case we need to do 'nothing'.
        boost::shared_ptr<BoundaryConditionsContainer<DIM, DIM, 2> > p_bcc;
        archive >> p_bcc;
        if (mpElectrodes->GetBoundaryConditionsContainer() != this->mpBoundaryConditionsContainer)
        {
            // The BCs will only actually be different if using a distributed tetrahedral mesh
            DistributedTetrahedralMesh<DIM,DIM>* p_dist_mesh = dynamic_cast<DistributedTetrahedralMesh<DIM,DIM>*>(this->mpMesh);
            if (p_dist_mesh)
            {
                mpElectrodes->GetBoundaryConditionsContainer()->MergeFromArchive(archive, this->mpMesh);
            }
            else
            {
                // Load into the temporary container, which will get thrown away shortly
                p_bcc->LoadFromArchive(archive, this->mpMesh);
                /// \todo #1159 sanity check that the contents of p_bcc and mpElectrodes->GetBoundaryConditionsContainer() match.
            }
        }
    }
}

#endif /*BIDOMAINPROBLEM_HPP_*/
