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


#ifndef ABSTRACTCARDIACPROBLEM_HPP_
#define ABSTRACTCARDIACPROBLEM_HPP_

#include <string>
#include <vector>
#include <cassert>
#include <climits>
#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "ClassIsAbstract.hpp"
#include "ChasteSerializationVersion.hpp"

#include "AbstractTetrahedralMesh.hpp"
#include "AbstractCardiacCell.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCardiacTissue.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "DistributedVectorFactory.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5DataWriter.hpp"
#include "Warnings.hpp"
#include "AbstractOutputModifier.hpp"
/*
 * Archiving extravaganza:
 *
 * We archive mesh and tissue through a pointer to an abstract class.  All the potential concrete
 * classes need to be included here, so they are registered with boost.  If not, boost won't be
 * able to find the archiving methods of the concrete class and will throw the following
 * exception:
 *
 *       terminate called after throwing an instance of 'boost::archive::archive_exception'
 *       what():  unregistered class
 *
 * No member variable is defined to be of any of these clases, so removing them won't
 * produce any compiler error.  The exception above will occur at runtime.
 *
 * This might not be even necessary in certain cases, if the file is included implicitly by another
 * header file or by the test itself.  It's safer though.
 */
#include "DistributedTetrahedralMesh.hpp"
#include "TetrahedralMesh.hpp" //May be needed for unarchiving a mesh
#include "MonodomainTissue.hpp"
#include "BidomainTissue.hpp"

/**
 * Empty untemplated base class so that CardiacSimulation::GetSavedProblem can work.
 */
class AbstractUntemplatedCardiacProblem : private boost::noncopyable
{
public:
    /** Virtual destructor to make this class polymorphic */
    virtual ~AbstractUntemplatedCardiacProblem()
    {}
};

/**
 * Base class for cardiac problems;
 * contains code generic to mono-/bi-domain and bidomain-with-bath.
 *
 * This class contains the tissue (PDEs and 'cells' ODEs),
 * boundary conditions, and postprocessing/results writers.
 *
 * It is called by CardiacSimulation, which is the outer wrapper class
 * for running a cardiac simulation, used by the executable.
 *
 * Many non-standard simulations will use this class directly,
 * and this is the preferred method for non-executable users.
 * See tutorials for usage.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractCardiacProblem : public AbstractUntemplatedCardiacProblem
{
    friend class TestBidomainWithBath;
    friend class TestMonodomainProblem;
    friend class TestCardiacSimulationArchiver;

    /** To save typing */
    typedef typename boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> >
        BccType;

private:
    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Save the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        if (version >= 1)
        {
            const unsigned element_dim=ELEMENT_DIM;
            archive & element_dim;
            const unsigned space_dim=SPACE_DIM;
            archive & space_dim;
            const unsigned problem_dim=PROBLEM_DIM;
            archive & problem_dim;
        }
        archive & mMeshFilename;
        archive & mpMesh;
        //archive & mAllocatedMemoryForMesh; // Mesh is deleted by AbstractCardiacTissue

        // We shouldn't ever have to save the old version
        assert(version >= 2);
//        {
//            bool use_matrix_based_assembly = true;
//            archive & use_matrix_based_assembly;
//        }

        archive & mWriteInfo;
        archive & mPrintOutput;
        archive & mNodesToOutput;
        //archive & mVoltageColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mExtraVariablesId; // Created by InitialiseWriter, called from Solve
        //archive & mTimeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mNodeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mpWriter; // Created by InitialiseWriter, called from Solve
        archive & mpCardiacTissue;
        //archive & mpSolver; // Only exists during calls to the Solve method
        bool has_solution = (mSolution != NULL);
        archive & has_solution;
        if (has_solution)
        {
            /// \todo #1317 code for saving/loading mSolution is PROBLEM_DIM specific, move it into the save/load methods for Mono and BidomainProblem.
            /// Note that extended_bidomain has its own version of this code.
            Hdf5DataWriter writer(*mpMesh->GetDistributedVectorFactory(), ArchiveLocationInfo::GetArchiveRelativePath(), "AbstractCardiacProblem_mSolution", false);
            writer.DefineFixedDimension(mpMesh->GetDistributedVectorFactory()->GetProblemSize());
            writer.DefineUnlimitedDimension("Time", "msec", 1);

            int vm_col = writer.DefineVariable("Vm","mV");

            if (PROBLEM_DIM==1)
            {
                writer.EndDefineMode();
                writer.PutUnlimitedVariable(0.0);
                writer.PutVector(vm_col, mSolution);
            }

            if (PROBLEM_DIM==2)
            {
                int phie_col = writer.DefineVariable("Phie","mV");
                std::vector<int> variable_ids;
                variable_ids.push_back(vm_col);
                variable_ids.push_back(phie_col);
                writer.EndDefineMode();
                writer.PutUnlimitedVariable(0.0);
                writer.PutStripedVector(variable_ids, mSolution);
            }

            writer.Close();

        }
        archive & mCurrentTime;

        // Save boundary conditions
        SaveBoundaryConditions(archive, mpMesh, mpBoundaryConditionsContainer);
        SaveBoundaryConditions(archive, mpMesh, mpDefaultBoundaryConditionsContainer);

        if (version >= 3)
        {
            archive & mOutputModifiers;
        }

        if (version >= 4)
        {
            archive & mUseHdf5DataWriterCache;
            archive & mHdf5DataWriterChunkSizeAndAlignment;
        }
    }

    /**
     * Load the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        if (version >= 1)
        {
            unsigned element_dim;
            unsigned space_dim;
            unsigned problem_dim;
            archive & element_dim;
            archive & space_dim;
            archive & problem_dim;
            if ((element_dim != ELEMENT_DIM) ||(space_dim != SPACE_DIM) ||(problem_dim != PROBLEM_DIM))
            {
                /*If we carry on from this point then the mesh produced by unarchiving from the
                 * archive is templated as AbstractTetrahedralMesh<element_dim, space_dim>
                 * which doesn't match AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*  mpMesh.
                 * Boost will through away the unarchived one, without deleting it properly and
                 * then set mpMesh=NULL.  We need to avoid this happening by bailing out.
                 */
                EXCEPTION("Failed to load from checkpoint because the dimensions of the archive do not match the object it's being read into.");
            }
        }
        archive & mMeshFilename;
        archive & mpMesh;
        assert(mpMesh != NULL); //If NULL then loading mesh has failed without an exception so Boost has given up on the mesh.  This would happen if a 2-dimensional mesh was successfully unarchived but mpMesh was expecting a 3-d mesh etc.
        //archive & mAllocatedMemoryForMesh; // Will always be true after a load

        if (version < 2)
        {
            bool use_matrix_based_assembly;
            archive & use_matrix_based_assembly;
        }

        archive & mWriteInfo;
        archive & mPrintOutput;
        archive & mNodesToOutput;
        //archive & mVoltageColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mExtraVariablesId; // Created by InitialiseWriter, called from Solve
        //archive & mTimeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mNodeColumnId; // Created by InitialiseWriter, called from Solve
        //archive & mpWriter; // Created by InitialiseWriter, called from Solve
        archive & mpCardiacTissue;
        //archive & mpSolver; // Only exists during calls to the Solve method
        bool has_solution;
        archive & has_solution;
        if ((has_solution) && PROBLEM_DIM < 3)
        {
            /// \todo #1317 code for saving/loading mSolution is PROBLEM_DIM specific, move it into the save/load methods for Mono and BidomainProblem.  (ExtendedBidomain has its own already.)
            /// \todo #1317 is there a reason we can't use PETSc's load/save vector functionality?
            mSolution = mpMesh->GetDistributedVectorFactory()->CreateVec(PROBLEM_DIM);
            DistributedVector mSolution_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(mSolution);

            Vec vm = mpMesh->GetDistributedVectorFactory()->CreateVec();
            Vec phie = mpMesh->GetDistributedVectorFactory()->CreateVec();

            std::string archive_dir = ArchiveLocationInfo::GetArchiveRelativePath();
            Hdf5DataReader reader(archive_dir, "AbstractCardiacProblem_mSolution", !FileFinder::IsAbsolutePath(archive_dir));
            reader.GetVariableOverNodes(vm, "Vm", 0);

            if (PROBLEM_DIM==1)
            {
                //reader.Close(); // no need to call close explicitly, done in the destructor

                DistributedVector vm_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(vm);

                DistributedVector::Stripe mSolution_vm(mSolution_distri,0);

                for (DistributedVector::Iterator index = mSolution_distri.Begin();
                     index != mSolution_distri.End();
                     ++index)
                {
                    mSolution_vm[index] = vm_distri[index];
                }
            }

            if (PROBLEM_DIM==2)
            {
                reader.GetVariableOverNodes(phie, "Phie", 0);
                //reader.Close(); // no need to call close explicitly, done in the destructor

                DistributedVector vm_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(vm);
                DistributedVector phie_distri = mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(phie);

                DistributedVector::Stripe mSolution_vm(mSolution_distri,0);
                DistributedVector::Stripe mSolution_phie(mSolution_distri,1);

                for (DistributedVector::Iterator index = mSolution_distri.Begin();
                     index != mSolution_distri.End();
                     ++index)
                {
                    mSolution_vm[index] = vm_distri[index];
                    mSolution_phie[index] = phie_distri[index];
                }
            }

            mSolution_distri.Restore();

            PetscTools::Destroy(vm);
            PetscTools::Destroy(phie);

        }
        archive & mCurrentTime;

        // Load boundary conditions
        mpBoundaryConditionsContainer = LoadBoundaryConditions(archive, mpMesh);
        mpDefaultBoundaryConditionsContainer = LoadBoundaryConditions(archive, mpMesh);

        if (version >= 3)
        {
            archive & mOutputModifiers;
        }

        if (version >= 4)
        {
            archive & mUseHdf5DataWriterCache;
            archive & mHdf5DataWriterChunkSizeAndAlignment;
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * Serialization helper method to save a boundary conditions container.
     *
     * @param archive  the archive to save to
     * @param pMesh  the mesh boundary conditions are defined on
     * @param pBcc  the container to save
     */
    template<class Archive>
    void SaveBoundaryConditions(Archive & archive,
                                AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                BccType pBcc) const
    {
        (*ProcessSpecificArchive<Archive>::Get()) & pBcc;
    }

    /**
     * Serialization helper method to load a boundary conditions container.
     *
     * @param archive  the archive to load from
     * @param pMesh  the mesh boundary conditions are to be defined on
     * @return  the loaded container
     */
    template<class Archive>
    BccType LoadBoundaryConditions(
            Archive & archive,
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    {
        // Load pointer from archive
        BccType p_bcc;
        (*ProcessSpecificArchive<Archive>::Get()) & p_bcc;

        // Fill in the conditions, if we have a container and it's not already full
        if (p_bcc)
        {
            p_bcc->LoadFromArchive(*ProcessSpecificArchive<Archive>::Get(), pMesh);
        }

        return p_bcc;
    }

protected:
    /** Meshes can be read from file or instantiated and passed directly to this
     *  class, this is for the former */
    std::string mMeshFilename;

    /** Whether this problem class has created the mesh itself, as opposed to being given it */
    bool mAllocatedMemoryForMesh;
    /** Whether to print some statistics (max/min voltage) to screen during the simulation */
    bool mWriteInfo;
    /** Whether to write any output at all */
    bool mPrintOutput;

    /** If only outputing voltage for selected nodes, which nodes to output at */
    std::vector<unsigned> mNodesToOutput;

    /** Used by the writer */
    unsigned mVoltageColumnId;
    /** List of extra variables to be written to HDF5 file */
    std::vector<unsigned> mExtraVariablesId;
    /** Used by the writer */
    unsigned mTimeColumnId;
    /** Used by the writer */
    unsigned mNodeColumnId;

    /** The monodomain or bidomain pde */
    AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* mpCardiacTissue;

    /** Boundary conditions container used in the simulation */
    BccType mpBoundaryConditionsContainer;
    /** It is convenient to also have a separate variable for default (zero-Neumann) boundary conditions */
    BccType mpDefaultBoundaryConditionsContainer;
    /** The PDE solver */
    AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* mpSolver;
    /** The cell factory creates the cells for each node */
    AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* mpCellFactory;
    /** The mesh. Can either by passed in, or the mesh filename can be set */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /** The current solution vector, of the form [V_0 .. V_N ] for monodomain and
     *  [V_0 phi_0 .. V_N phi_N] for bidomain */
    Vec mSolution;

    /**
     * The current simulation time.
     *
     * This is used to be able to restart simulations at a point other than time zero,
     * either because of repeated calls to Solve (with increased simulation duration)
     * or because of restarting from a checkpoint.
     */
    double mCurrentTime;

    /** Adaptivity controller (defaults to NULL). */
    AbstractTimeAdaptivityController* mpTimeAdaptivityController;

    /**
     * Subclasses must override this method to create a PDE object of the appropriate type.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     * @return a pointer to the newly created tissue
     */
    virtual AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* CreateCardiacTissue() =0;

    /**
     * Subclasses must override this method to create a suitable solver object.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     * @return pointer to newly created PDE solver
     */
    virtual AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>* CreateSolver() =0;

    /**
     * Subclasses must override this method to create a suitable mesh object.
     *
     * Only needed if the subclass needs something other than a DistributedTetrahedralMesh.
     *
     * This class will take responsibility for freeing the object when it is finished with.
     */
    virtual void CreateMeshFromHeartConfig();

    /**
     * CardiacElectroMechanicsProblem needs access to #mpWriter.
     */
    template<unsigned DIM, unsigned ELEC_PROB_DIM>
    friend class CardiacElectroMechanicsProblem;

    /**
     * The object to use to write results to disk.
     */
    Hdf5DataWriter* mpWriter;

    /**
     * Whether to instruct the writer to cache writes.
     */
    bool mUseHdf5DataWriterCache;

    /**
     * Size to pass to Hdf5DataWriter for chunk size and alignment.
     */
    hsize_t mHdf5DataWriterChunkSizeAndAlignment;

    /**
     * A vector of user-defined output modifiers which may be used to produce lightweight on the fly output
     */
    std::vector<boost::shared_ptr<AbstractOutputModifier> > mOutputModifiers;

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     * create cells.
     */
    AbstractCardiacProblem(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory);

    /**
     * Constructor used by archiving.
     */
    AbstractCardiacProblem();

    /**
     *  Destructor
     */
    virtual ~AbstractCardiacProblem();

    /**
     * Initialise the system, once parameters have been set up.
     *
     * Must be called before first calling Solve().  If loading from a checkpoint,
     * do NOT call this method, as it can also be used to reset the problem to
     * perform another simulation from time 0.
     */
    void Initialise();

    /**
     *  Set a file from which the nodes for each processor are read
     *
     * @param rFilename
     */
    void SetNodesPerProcessorFilename(const std::string& rFilename);

    /**
     *  Set the boundary conditions container.
     *  @param pBcc is a pointer to a boundary conditions container
     */
    void SetBoundaryConditionsContainer(BccType pBcc);

    /**
     *  Performs a series of checks before solving.
     *  It checks whether the cardiac pde has been defined,
     *  whether the simulation time is greater than zero and
     *  whether the output directory is specified (or the output is set not to be produced).
     *  It throws exceptions if any of the above checks fails.
     */
    virtual void PreSolveChecks();

    /**
     *
     * This method sets the initial condition for the PDE by getting the
     * voltages (V) from the cell models at the nodes.
     *
     * If the problem dimension is two (Bidomain) the second variable (phi_e) is set to zero.
     *
     * This is virtual so BidomainProblem can overwrite V to zero for bath nodes, if
     * there are any.
     *
     * @return the newly created intial conditions vector
     * \todo Perhaps this should be a method of AbstractCardiacTissue??
     */
    virtual Vec CreateInitialCondition();

    /**
     * This only needs to be called if a mesh filename has not been set.
     *
     * @param pMesh  the mesh object to use
     */
    void SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     *  Set whether the simulation will generate results files.
     *
     * @param rPrintOutput
     */
    void PrintOutput(bool rPrintOutput);

    /**
     *  Set whether extra info will be written to stdout during computation.
     *
     * @param writeInfo
     */
    void SetWriteInfo(bool writeInfo = true);

    /**
     *  @return the final solution vector. This vector is distributed over all processes.
     *
     *  In case of Bidomain, this is of length 2*numNodes, and of the form
     *  (V_1, phi_1, V_2, phi_2, ......, V_N, phi_N).
     *  where V_j is the voltage at node j and phi_j is the
     *  extracellular potential at node j.
     *
     *  Use with caution since we don't want to alter the state of the PETSc vector.
     */
    Vec GetSolution();

    /**
     * @return the solution vector, wrapped in a DistributedVector.
     *
     * See also GetSolution.
     */
    DistributedVector GetSolutionDistributedVector();

    /**
     * @return the current time of the simulation
     */
    double GetCurrentTime();

    /**
     * @return the mesh used
     */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM> & rGetMesh();

    /**
     * @return the cardiac tissue object used
     */
    AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* GetTissue();

    /**
     *  First performs some checks by calling  the PreSolveChecks method.
     *  It creates an solver to which it passes the boundary conditions specified by the user
     *  (otherwise it passes the defauls bcc).
     *   It then calls the Solve method on the solver class.
     *   It also handles the output, if necessary.
     *
     * @note This method is collective, and hence must be called by all processes.
     */
    void Solve();

    /**
     * Closes the files where the solution is stored and,
     * if specified so (as it is by default), converts the output to Meshalyzer format
     * by calling the WriteFilesUsingMesh method in the MeshalyzerWriter class.
     *
     * @note This method is collective, and hence must be called by all processes.
     */
    void CloseFilesAndPostProcess();

    /**
     * Write informative details about the progress of the simulation to standard output.
     *
     * Implemented only in subclasses.
     *
     * @param time  the current time
     */
    virtual void WriteInfo(double time)=0;

    /**
     * Define what variables are written to the primary results file.
     * @param extending  whether we are extending an existing results file
     */
    virtual void DefineWriterColumns(bool extending);

    /**
     * Define the user specified variables to be written to the primary results file
     * @param extending  whether we are extending an existing results file
     */
    void DefineExtraVariablesWriterColumns(bool extending);

    /**
     * Write one timestep of output data to the primary results file.
     *
     * @param time  the current time
     * @param voltageVec  the solution vector to write
     */
    virtual void WriteOneStep(double time, Vec voltageVec) = 0;

    /**
     * Write one timestep of output data for the extra variables to the primary results file.
     */
    void WriteExtraVariablesOneStep();

    /**
     * It creates and initialises the hdf writer from the Hdf5DataWriter class.
     * It passes the output directory and file name to it.
     * It is called by Solve(), if the output needs to be generated.
     *
     * This method will try to open an existing .h5 file for extension if
     * we are loading from an archive and one is present.
     *
     * @return whether the writer is outputting to an existing file.
     */
    bool InitialiseWriter();

    /**
     * Set whether to use caching in the Hdf5DataWriter. This tells the
     * Hdf5DataWriter to write only whole chunks to disk, rather than every
     * print timestep, which is much faster when running with many processes.
     * @param useCache Whether to use the cache
     */
    void SetUseHdf5DataWriterCache(bool useCache=true);

    /**
     * Set Hdf5DataWriter target chunk size and alignment parameters.
     *
     * The most likely use case for this is setting it to the stripe size on a
     * striped filesystem. This results in the writer choosing chunk dimensions
     * that should efficiently fit within a stripe, AND padding the chunks so
     * each chunk fits in one stripe.
     *
     * For example, if your filesystem uses 1 M stripes, call this method with
     * argument 1048576 (or 0x100000 if you like round numbers).
     *
     * NOTE: The alignment parameter is only used for NEW HDF5 files. The chunk
     *       size is only used for NEW datasets (e.g. results or postprocessing)
     *       and NOT when EXTENDING a dataset.
     *       In other words, when adding a dataset to a file the alignment
     *       parameter will be ignored. When adding to a dataset, both will be
     *       ignored. etc.
     * @param size size in bytes to use for target chunk size and alignment
     */
    void SetHdf5DataWriterTargetChunkSizeAndAlignment(hsize_t size);

    /**
     * Specifies which nodes in the mesh to output. This method must be called before InitialiseWriter,
     * otherwise all nodes will still be output. If this method is called when extending an existing
     * HDF5 file, it will be ignored.
     *
     * @param rNodesToOutput is a reference to a vector with the indexes of the nodes
     * where the output is desired.
     * If empty, the output will be for all the nodes in the mesh.
     */
    void SetOutputNodes(std::vector<unsigned> & rNodesToOutput);

    /**
     * @return a newly created data reader configured to read the results we've been outputting.
     */
    Hdf5DataReader GetDataReader();

    /**
     *  Called at beginning of each time step in the main time-loop in
     *  Solve(). Empty implementation but can be overloaded by child
     *  classes.
     *
     * @param time  the current time
     */
    virtual void AtBeginningOfTimestep(double time)
    {}

    /**
     *  Called at end of each time step in the main time-loop in
     *  Solve(). Empty implementation but can be overloaded by child
     *  classes.
     *
     * @param time  the current time
     */
    virtual void OnEndOfTimestep(double time)
    {}

    /**
     * Allow subclasses to define additional 'stopping times' for the printing
     * time step loop.  This allows bidomain simulations to specify exactly
     * when the Electrodes should be turned on or off.
     *
     * @param rAdditionalStoppingTimes  to be filled in with the additional stopping times
     */
    virtual void SetUpAdditionalStoppingTimes(std::vector<double>& rAdditionalStoppingTimes)
    {}

    ///\todo #1704 add default adaptivity controller and allow the user just to call with true
    // and no controller, in which case the default is used.

    /**
     *  Set whether (or not) to use a time adaptivity controller
     *  @param useAdaptivity whether to use adaptivity
     *  @param pController The controller (only relevant if useAdaptivity==true)
     */
    void SetUseTimeAdaptivityController(bool useAdaptivity,
                                        AbstractTimeAdaptivityController* pController = NULL);

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
     *  -# (via #mpMesh) DistributedVectorFactory*
     *  -# (via #mpCardiacTissue LoadCardiacCells) DistributedVectorFactory*
     *  -# (via #mpCardiacTissue LoadCardiacCells) number_of_cells and sequence of AbstractCardiacCell*, possibly with Purkinje interleaved
     *  -# (via #mpCardiacTissue) DistributedVectorFactory*
     *  -# #mpBoundaryConditionsContainer
     *  -# #mpDefaultBoundaryConditionsContainer
     *  -# (if we're a BidomainProblem) stuff in BidomainProblem::LoadExtraArchiveForBidomain
     */
    template<class Archive>
    void LoadExtraArchive(Archive & archive, unsigned version);

    /**
     * @return whether there's bath defined in this problem
     */
    virtual bool GetHasBath();

    /**
     *  Set an electrode object (which provides boundary conditions). Only
     *  valid if there is a bath.
     */
    virtual void SetElectrodes();

    /**
     * Add an output modifier onto a list of such objects.  These will be processed in the order in which they have been given.
     * The modifier should not be destroyed before the solve loop has completed
     * @param pOutputModifier  Pointer to the modifier to be added
     */
    void AddOutputModifier( boost::shared_ptr<AbstractOutputModifier> pOutputModifier)
    {
        mOutputModifiers.push_back(pOutputModifier);
    }
};

TEMPLATED_CLASS_IS_ABSTRACT_3_UNSIGNED(AbstractCardiacProblem)


template<unsigned DIM>
class BidomainProblem;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
template<class Archive>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::LoadExtraArchive(Archive & archive, unsigned version)
{
    // The vector factory must be loaded, but isn't needed for anything.
    DistributedVectorFactory* p_mesh_factory;
    archive >> p_mesh_factory;

    // How many processes were used by the saving simulation?
    DistributedVectorFactory* p_original_factory = p_mesh_factory->GetOriginalFactory();
    unsigned orig_num_procs = 1;
    if (p_original_factory)
    {
        orig_num_procs = p_original_factory->GetNumProcs();
    }

    // The cardiac cells - load only the cells we actually own
    mpCardiacTissue->LoadCardiacCells(archive, version);

    {
        DistributedVectorFactory* p_pde_factory;
        archive >> p_pde_factory;
        assert(p_pde_factory == p_mesh_factory); // Paranoia...
    }
    // We no longer need this vector factory, since we already have our own.
    delete p_mesh_factory;

    // The boundary conditions
    BccType p_bcc;
    archive >> p_bcc;
    if (p_bcc)
    {
        if (!mpBoundaryConditionsContainer)
        {
            mpBoundaryConditionsContainer = p_bcc;
            mpBoundaryConditionsContainer->LoadFromArchive(archive, mpMesh);
        }
        else
        {
            // If the mesh which was archived was a TetrahedralMesh then we have all the boundary conditions
            // in every process-specific archive.  We no longer test for this.
// LCOV_EXCL_START
            if (!dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>*>(mpMesh) && orig_num_procs > 1)
            {
                // The correct way to do this should be:
                // p_bcc->LoadFromArchive(archive, mpMesh);
                WARNING("Loading from a parallel archive which used a non-distributed mesh.  This scenario should work but is not fully tested.");
            }
// LCOV_EXCL_STOP
            mpBoundaryConditionsContainer->MergeFromArchive(archive, mpMesh);
        }
    }
    BccType p_default_bcc;
    archive >> p_default_bcc;
    if (p_default_bcc)
    {
        // This always holds, so we never need to load the default BCs
        assert(p_bcc == p_default_bcc);
    }

    // Are we a bidomain problem?
    BidomainProblem<ELEMENT_DIM>* p_bidomain_problem = dynamic_cast<BidomainProblem<ELEMENT_DIM>*>(this);
    if (p_bidomain_problem)
    {
        assert(ELEMENT_DIM == SPACE_DIM);
        p_bidomain_problem->LoadExtraArchiveForBidomain(archive, version);
    }
}

namespace boost {
namespace serialization {
/**
 * Specify a version number for archive backwards compatibility.
 *
 * This is how to do BOOST_CLASS_VERSION(AbstractCardiacProblem, 1)
 * with a templated class.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM,  unsigned PROBLEM_DIM>
struct version<AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> >
{
    ///Macro to set the version number of templated archive in known versions of Boost
    CHASTE_VERSION_CONTENT(4);
};
} // namespace serialization
} // namespace boost
#endif /*ABSTRACTCARDIACPROBLEM_HPP_*/
