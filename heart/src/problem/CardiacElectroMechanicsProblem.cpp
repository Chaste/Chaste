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

#include "CardiacElectroMechanicsProblem.hpp"

#include "OutputFileHandler.hpp"
#include "ReplicatableVector.hpp"
#include "HeartConfig.hpp"
#include "LogFile.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "TimeStepper.hpp"
#include "TrianglesMeshWriter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "PetscTools.hpp"
#include "ImplicitCardiacMechanicsSolver.hpp"
#include "ExplicitCardiacMechanicsSolver.hpp"
#include "CmguiDeformedSolutionsWriter.hpp"
#include "VoltageInterpolaterOntoMechanicsMesh.hpp"
#include "BidomainProblem.hpp"


template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::DetermineWatchedNodes()
{
    assert(mIsWatchedLocation);

    // find the nearest electrics mesh node
    double min_dist = DBL_MAX;
    unsigned node_index = UNSIGNED_UNSET;
    for (unsigned i=0; i<mpElectricsMesh->GetNumNodes(); i++)
    {
        double dist = norm_2(mWatchedLocation - mpElectricsMesh->GetNode(i)->rGetLocation());
        if (dist < min_dist)
        {
            min_dist = dist;
            node_index = i;
        }
    }

    // set up watched node, if close enough
    assert(node_index != UNSIGNED_UNSET); // should def have found something
    c_vector<double,DIM> pos = mpElectricsMesh->GetNode(node_index)->rGetLocation();

    if (min_dist > 1e-8)
    {
        // LCOV_EXCL_START
        std::cout << "ERROR: Could not find an electrics node very close to requested watched location - "
                  << "min distance was " << min_dist << " for node " << node_index
                  << " at location " << pos << std::flush;;

        //// the following causes a seg fault for some reason (!!???!!!)
        //EXCEPTION("Could not find an electrics node very close to requested watched location");
        NEVER_REACHED;
        // LCOV_EXCL_STOP
    }
    else
    {
        LOG(2,"Chosen electrics node "<<node_index<<" at location " << pos << " to be watched");
        mWatchedElectricsNodeIndex = node_index;
    }

    // find nearest mechanics mesh
    min_dist = DBL_MAX;
    node_index = UNSIGNED_UNSET;
    c_vector<double,DIM> pos_at_min;

    for (unsigned i=0; i<mpMechanicsMesh->GetNumNodes(); i++)
    {
        c_vector<double,DIM> position = mpMechanicsMesh->GetNode(i)->rGetLocation();

        double dist = norm_2(position-mWatchedLocation);

        if (dist < min_dist)
        {
            min_dist = dist;
            node_index = i;
            pos_at_min = position;
        }
    }

    // set up watched node, if close enough
    assert(node_index != UNSIGNED_UNSET); // should def have found something

    if (min_dist > 1e-8)
    {
        // LCOV_EXCL_START
        std::cout << "ERROR: Could not find a mechanics node very close to requested watched location - "
                  << "min distance was " << min_dist << " for node " << node_index
                  << " at location " << pos_at_min;

        //// the following causes a seg fault for some reason (!!???!!!)
        //EXCEPTION("Could not find a mechanics node very close to requested watched location");
        NEVER_REACHED;
        // LCOV_EXCL_STOP
    }
    else
    {
        LOG(2,"Chosen mechanics node "<<node_index<<" at location " << pos << " to be watched");
        mWatchedMechanicsNodeIndex = node_index;
    }

    OutputFileHandler handler(mOutputDirectory,false);
    mpWatchedLocationFile = handler.OpenOutputFile("watched.txt");
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::WriteWatchedLocationData(double time, Vec voltage)
{
    assert(mIsWatchedLocation);

    std::vector<c_vector<double,DIM> >& deformed_position = mpMechanicsSolver->rGetDeformedPosition();

    ///\todo Improve efficiency of this method?
    ReplicatableVector voltage_replicated(voltage);
    double V = voltage_replicated[mWatchedElectricsNodeIndex];

    //// Removed the following which also took this nodes calcium concentration and printing, because (it isn't that
    //// important and) it won't work in parallel and has the hardcoded index issue described below.
    //     // Metadata is currently being added to CellML models and then this will be avoided by asking for Calcium.
    //    double Ca = mpElectricsProblem->GetMonodomainTissue()->GetCardiacCell(mWatchedElectricsNodeIndex)->GetIntracellularCalciumConcentration();

    *mpWatchedLocationFile << time << " ";
    for (unsigned i=0; i<DIM; i++)
    {
        *mpWatchedLocationFile << deformed_position[mWatchedMechanicsNodeIndex](i) << " ";
    }
    *mpWatchedLocationFile << V << "\n";
    mpWatchedLocationFile->flush();
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
c_matrix<double,DIM,DIM>& CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::rCalculateModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,DIM,DIM>& rOriginalConductivity, unsigned domainIndex)
{

    // first get the deformation gradient for this electrics element
    unsigned containing_mechanics_elem = mpMeshPair->rGetCoarseElementsForFineElementCentroids()[elementIndex];
    c_matrix<double,DIM,DIM>& r_deformation_gradient = mDeformationGradientsForEachMechanicsElement[containing_mechanics_elem];

    // compute sigma_def = F^{-1} sigma_undef F^{-T}
    c_matrix<double,DIM,DIM> inv_F = Inverse(r_deformation_gradient);
    mModifiedConductivityTensor = prod(inv_F, rOriginalConductivity);
    mModifiedConductivityTensor = prod(mModifiedConductivityTensor, trans(inv_F));

    return mModifiedConductivityTensor;
}


/**
 * Helper "function" for the constructor, to create the electrics sub-problem without dynamic_cast.
 */
template<unsigned DIM, unsigned PROBLEM_DIM>
class CreateElectricsProblem
{
public:
    /**
     * The actual creation method.
     * @param problemType
     * @param pCellFactory
     * @return newly created cardiac problem
     */
    static AbstractCardiacProblem<DIM, DIM, PROBLEM_DIM>* Create(ElectricsProblemType problemType,
                                                                 AbstractCardiacCellFactory<DIM>* pCellFactory);
};

/**
 * Specialization for monodomain problems
 */
template<unsigned DIM>
class CreateElectricsProblem<DIM, 1u>
{
public:
    /**
     * The actual creation method.
     * @param problemType
     * @param pCellFactory
     * @return newly created cardiac problem
     */
    static AbstractCardiacProblem<DIM, DIM, 1u>* Create(ElectricsProblemType problemType,
                                                        AbstractCardiacCellFactory<DIM>* pCellFactory)
    {
        if (problemType == MONODOMAIN)
        {
            return new MonodomainProblem<DIM>(pCellFactory);
        }
        EXCEPTION("The second template parameter should be 2 when a bidomain problem is chosen");
    }
};

/**
 * Specialization for bidomain problems
 */
template<unsigned DIM>
class CreateElectricsProblem<DIM, 2u>
{
public:
    /**
     * The actual creation method.
     * @param problemType
     * @param pCellFactory
     * @return newly created cardiac problem
     */
    static AbstractCardiacProblem<DIM, DIM, 2u>* Create(ElectricsProblemType problemType,
                                                        AbstractCardiacCellFactory<DIM>* pCellFactory)
    {
        if (problemType == BIDOMAIN)
        {
            return new BidomainProblem<DIM>(pCellFactory, false);//false-> no bath
        }
        if (problemType == BIDOMAIN_WITH_BATH)
        {
            return new BidomainProblem<DIM>(pCellFactory, true);// true-> bath
        }
        EXCEPTION("The second template parameter should be 1 when a monodomain problem is chosen");
    }
};


template<unsigned DIM, unsigned ELEC_PROB_DIM>
CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::CardiacElectroMechanicsProblem(
            CompressibilityType compressibilityType,
            ElectricsProblemType electricsProblemType,
            TetrahedralMesh<DIM,DIM>* pElectricsMesh,
            QuadraticMesh<DIM>* pMechanicsMesh,
            AbstractCardiacCellFactory<DIM>* pCellFactory,
            ElectroMechanicsProblemDefinition<DIM>* pProblemDefinition,
            std::string outputDirectory)
      : mCompressibilityType(compressibilityType),
        mpCardiacMechSolver(NULL),
        mpMechanicsSolver(NULL),
        mpElectricsMesh(pElectricsMesh),
        mpMechanicsMesh(pMechanicsMesh),
        mpProblemDefinition(pProblemDefinition),
        mHasBath(false),
        mpMeshPair(NULL),
        mNoElectricsOutput(false),
        mIsWatchedLocation(false),
        mWatchedElectricsNodeIndex(UNSIGNED_UNSET),
        mWatchedMechanicsNodeIndex(UNSIGNED_UNSET),
        mNumTimestepsToOutputDeformationGradientsAndStress(UNSIGNED_UNSET)
{
    // Do some initial set up...
    // However, NOTE, we don't use either the passed in meshes or the problem_definition.
    // These pointers are allowed to be NULL, in case a child constructor wants to set
    // them up (eg CardiacElectroMechProbRegularGeom).
    // The meshes and problem_defn are used for the first time in Initialise().


    // Start-up mechanics event handler..
    MechanicsEventHandler::Reset();
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ALL);
    // disable the electric event handler, because we use a problem class but
    // don't call Solve, so we would have to worry about starting and ending any
    // events in AbstractCardiacProblem::Solve() (esp. calling EndEvent(EVERYTHING))
    // if we didn't disable it.
    HeartEventHandler::Disable();

    assert(HeartConfig::Instance()->GetSimulationDuration()>0.0);
    assert(HeartConfig::Instance()->GetPdeTimeStep()>0.0);

    // Create the monodomain problem.
    // **NOTE** WE ONLY USE THIS TO: set up the cells, get an initial condition
    // (voltage) vector, and get a solver. We won't ever call Solve on the cardiac problem class
    assert(pCellFactory != NULL);
    mpElectricsProblem = CreateElectricsProblem<DIM,ELEC_PROB_DIM>::Create(electricsProblemType, pCellFactory);

    if (electricsProblemType == BIDOMAIN_WITH_BATH)
    {
        mHasBath = true;
    }
    // check whether output is required
    mWriteOutput = (outputDirectory!="");
    if (mWriteOutput)
    {
        mOutputDirectory = outputDirectory;
        // create the directory
        OutputFileHandler handler(mOutputDirectory);
        mDeformationOutputDirectory = mOutputDirectory + "/deformation";
        HeartConfig::Instance()->SetOutputDirectory(mOutputDirectory + "/electrics");
        HeartConfig::Instance()->SetOutputFilenamePrefix("voltage");
    }
    else
    {
        mDeformationOutputDirectory = "";
    }

//    mpImpactRegion=NULL;
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::~CardiacElectroMechanicsProblem()
{
    // NOTE if SetWatchedLocation but not Initialise has been called, mpWatchedLocationFile
    // will be uninitialised and using it will cause a seg fault. Hence the mpMechanicsMesh!=NULL
    // it is true if Initialise has been called.
    if (mIsWatchedLocation && mpMechanicsMesh)
    {
        mpWatchedLocationFile->close();
    }

    delete mpElectricsProblem;
    delete mpCardiacMechSolver;
    delete mpMeshPair;

    LogFile::Close();
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::Initialise()
{
    assert(mpElectricsMesh!=NULL);
    assert(mpMechanicsMesh!=NULL);
    assert(mpProblemDefinition!=NULL);
    assert(mpCardiacMechSolver==NULL);

    ///\todo This is fragile: check how the TimeStepper does it, and possibly refactor the behaviour there
    /// into a static helper method if it isn't already.
    mNumElecTimestepsPerMechTimestep = (unsigned) floor((mpProblemDefinition->GetMechanicsSolveTimestep()/HeartConfig::Instance()->GetPdeTimeStep())+0.5);
    if (fabs(mNumElecTimestepsPerMechTimestep*HeartConfig::Instance()->GetPdeTimeStep() - mpProblemDefinition->GetMechanicsSolveTimestep()) > 1e-6)
    {
        EXCEPTION("Electrics PDE timestep does not divide mechanics solve timestep");
    }

    // Create the Logfile (note we have to do this after the output dir has been
    // created, else the log file might get cleaned away
    std::string log_dir = mOutputDirectory; // just the TESTOUTPUT dir if mOutputDir="";
    LogFile::Instance()->Set(2, mOutputDirectory);
    LogFile::Instance()->WriteHeader("Electromechanics");
    LOG(2, DIM << "d Implicit CardiacElectroMechanics Simulation:");
    LOG(2, "End time = " << HeartConfig::Instance()->GetSimulationDuration() << ", electrics time step = " << HeartConfig::Instance()->GetPdeTimeStep() << ", mechanics timestep = " << mpProblemDefinition->GetMechanicsSolveTimestep() << "\n");
    LOG(2, "Contraction model ode timestep " << mpProblemDefinition->GetContractionModelOdeTimestep());
    LOG(2, "Output is written to " << mOutputDirectory << "/[deformation/electrics]");

    LOG(2, "Electrics mesh has " << mpElectricsMesh->GetNumNodes() << " nodes");
    LOG(2, "Mechanics mesh has " << mpMechanicsMesh->GetNumNodes() << " nodes");

    LOG(2, "Initialising..");


    if (mIsWatchedLocation)
    {
        DetermineWatchedNodes();
    }

    // initialise electrics problem
    mpElectricsProblem->SetMesh(mpElectricsMesh);
    mpElectricsProblem->Initialise();

    if (mCompressibilityType==INCOMPRESSIBLE)
    {
        switch(mpProblemDefinition->GetSolverType())
        {
            case EXPLICIT:
                mpCardiacMechSolver = new ExplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<DIM>,DIM>(
                                        *mpMechanicsMesh,*mpProblemDefinition,mDeformationOutputDirectory);
                break;
            case IMPLICIT:
                mpCardiacMechSolver = new ImplicitCardiacMechanicsSolver<IncompressibleNonlinearElasticitySolver<DIM>,DIM>(
                                        *mpMechanicsMesh,*mpProblemDefinition,mDeformationOutputDirectory);
                break;
            default:
                NEVER_REACHED;
        }
    }
    else
    {
        // repeat above with Compressible solver rather than incompressible -
        // not the neatest but avoids having to template this class.
        switch(mpProblemDefinition->GetSolverType())
        {
            case EXPLICIT:
                mpCardiacMechSolver = new ExplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<DIM>,DIM>(
                                        *mpMechanicsMesh,*mpProblemDefinition,mDeformationOutputDirectory);
                break;
            case IMPLICIT:
                mpCardiacMechSolver = new ImplicitCardiacMechanicsSolver<CompressibleNonlinearElasticitySolver<DIM>,DIM>(
                                        *mpMechanicsMesh,*mpProblemDefinition,mDeformationOutputDirectory);
                break;
            default:
                NEVER_REACHED;
        }
    }


    mpMechanicsSolver = dynamic_cast<AbstractNonlinearElasticitySolver<DIM>*>(mpCardiacMechSolver);
    assert(mpMechanicsSolver);

    // set up mesh pair and determine the fine mesh elements and corresponding weights for each
    // quadrature point in the coarse mesh
    mpMeshPair = new FineCoarseMeshPair<DIM>(*mpElectricsMesh, *mpMechanicsMesh);
    mpMeshPair->SetUpBoxesOnFineMesh();
    mpMeshPair->ComputeFineElementsAndWeightsForCoarseQuadPoints(*(mpCardiacMechSolver->GetQuadratureRule()), false);
    mpMeshPair->DeleteFineBoxCollection();

    mpCardiacMechSolver->SetFineCoarseMeshPair(mpMeshPair);
    mpCardiacMechSolver->Initialise();

    unsigned num_quad_points = mpCardiacMechSolver->GetTotalNumQuadPoints();
    mInterpolatedCalciumConcs.assign(num_quad_points, 0.0);
    mInterpolatedVoltages.assign(num_quad_points, 0.0);

    if (mpProblemDefinition->ReadFibreSheetDirectionsFromFile())
    {
       mpCardiacMechSolver->SetVariableFibreSheetDirections(mpProblemDefinition->GetFibreSheetDirectionsFile(),
                                                            mpProblemDefinition->GetFibreSheetDirectionsDefinedPerQuadraturePoint());
    }


    if (mpProblemDefinition->GetDeformationAffectsConductivity() || mpProblemDefinition->GetDeformationAffectsCellModels())
    {
        mpMeshPair->SetUpBoxesOnCoarseMesh();
    }


    if (mpProblemDefinition->GetDeformationAffectsCellModels() || mpProblemDefinition->GetDeformationAffectsConductivity())
    {
        // initialise the stretches saved for each mechanics element
        mStretchesForEachMechanicsElement.resize(mpMechanicsMesh->GetNumElements(), 1.0);

        // initialise the store of the F in each mechanics element (one constant value of F) in each
        mDeformationGradientsForEachMechanicsElement.resize(mpMechanicsMesh->GetNumElements(),identity_matrix<double>(DIM));
    }


    if (mpProblemDefinition->GetDeformationAffectsCellModels())
    {
        // compute the coarse elements which contain each fine node -- for transferring stretch from
        // mechanics solve electrics cell models
        mpMeshPair->ComputeCoarseElementsForFineNodes(false);

    }

    if (mpProblemDefinition->GetDeformationAffectsConductivity())
    {
        // compute the coarse elements which contain each fine element centroid -- for transferring F from
        // mechanics solve to electrics mesh elements
        mpMeshPair->ComputeCoarseElementsForFineElementCentroids(false);

        // tell the abstract tissue class that the conductivities need to be modified, passing in this class
        // (which is of type AbstractConductivityModifier)
        mpElectricsProblem->GetTissue()->SetConductivityModifier(this);
    }

    if (mWriteOutput)
    {
        TrianglesMeshWriter<DIM,DIM> mesh_writer(mOutputDirectory,"electrics_mesh",false);
        mesh_writer.WriteFilesUsingMesh(*mpElectricsMesh);
    }
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::Solve()
{
    // initialise the meshes and mechanics solver
    if (mpCardiacMechSolver==NULL)
    {
        Initialise();
    }

    bool verbose_during_solve = (    mpProblemDefinition->GetVerboseDuringSolve()
                                  || CommandLineArguments::Instance()->OptionExists("-mech_verbose")
                                  || CommandLineArguments::Instance()->OptionExists("-mech_very_verbose"));


    mpProblemDefinition->Validate();

    boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,ELEC_PROB_DIM> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,ELEC_PROB_DIM>);
    p_bcc->DefineZeroNeumannOnMeshBoundary(mpElectricsMesh, 0);
    mpElectricsProblem->SetBoundaryConditionsContainer(p_bcc);

    // get an electrics solver from the problem. Note that we don't call
    // Solve() on the CardiacProblem class, we do the looping here.
    AbstractDynamicLinearPdeSolver<DIM,DIM,ELEC_PROB_DIM>* p_electrics_solver
       = mpElectricsProblem->CreateSolver();

    // set up initial voltage etc
    Vec electrics_solution=NULL; //This will be set and used later
    Vec calcium_data= mpElectricsMesh->GetDistributedVectorFactory()->CreateVec();
    Vec initial_voltage = mpElectricsProblem->CreateInitialCondition();

    // write the initial position
    unsigned counter = 0;

    TimeStepper stepper(0.0, HeartConfig::Instance()->GetSimulationDuration(), mpProblemDefinition->GetMechanicsSolveTimestep());

    CmguiDeformedSolutionsWriter<DIM>* p_cmgui_writer = NULL;

    std::vector<std::string> variable_names;

    if (mWriteOutput)
    {
        mpMechanicsSolver->SetWriteOutput();
        mpMechanicsSolver->WriteCurrentSpatialSolution("undeformed","nodes");

        p_cmgui_writer = new CmguiDeformedSolutionsWriter<DIM>(mOutputDirectory+"/deformation/cmgui",
                                                               "solution",
                                                               *(this->mpMechanicsMesh),
                                                               WRITE_QUADRATIC_MESH);
        variable_names.push_back("V");
        if (ELEC_PROB_DIM==2)
        {
            variable_names.push_back("Phi_e");
            if (mHasBath==true)
            {
                std::vector<std::string> regions;
                regions.push_back("tissue");
                regions.push_back("bath");
                p_cmgui_writer->SetRegionNames(regions);
            }
        }
        p_cmgui_writer->SetAdditionalFieldNames(variable_names);
        p_cmgui_writer->WriteInitialMesh("undeformed");
    }

    /////////////////////////////////////////////////////////////////
    ////
    ////  Solve a static problem which might involve finding
    ////  equilibrium state given initial internal pressures etc
    ////
    ////////////////////////////////////////////////////////////////

    LOG(2, "\nSolving for initial deformation");
    // LCOV_EXCL_START
    if (verbose_during_solve)
    {
        std::cout << "\n\n ** Solving for initial deformation\n";
    }
    // LCOV_EXCL_STOP

    mpMechanicsSolver->SetWriteOutput(false);

    mpMechanicsSolver->SetCurrentTime(0.0);
    MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ALL_MECH);

    mpMechanicsSolver->SetIncludeActiveTension(false);
    if (mNumTimestepsToOutputDeformationGradientsAndStress!=UNSIGNED_UNSET)
    {
        mpMechanicsSolver->SetComputeAverageStressPerElementDuringSolve(true);
    }

    unsigned total_newton_iters = 0;
    for (unsigned index=1; index<=mpProblemDefinition->GetNumIncrementsForInitialDeformation(); index++)
    {
        // LCOV_EXCL_START
        if (verbose_during_solve)
        {
            std::cout << "    Increment " << index << " of " << mpProblemDefinition->GetNumIncrementsForInitialDeformation() << "\n";
        }
        // LCOV_EXCL_STOP

        if (mpProblemDefinition->GetTractionBoundaryConditionType()==PRESSURE_ON_DEFORMED)
        {
            mpProblemDefinition->SetPressureScaling(((double)index)/mpProblemDefinition->GetNumIncrementsForInitialDeformation());
        }
        mpMechanicsSolver->Solve();

        total_newton_iters += mpMechanicsSolver->GetNumNewtonIterations();
    }

    mpMechanicsSolver->SetIncludeActiveTension(true);
    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ALL_MECH);
    LOG(2, "    Number of newton iterations = " << total_newton_iters);


    unsigned mech_writer_counter = 0;

    if (mWriteOutput)
    {
        LOG(2, "  Writing output");
        mpMechanicsSolver->SetWriteOutput();
        mpMechanicsSolver->WriteCurrentSpatialSolution("solution","nodes",mech_writer_counter);
        p_cmgui_writer->WriteDeformationPositions(rGetDeformedPosition(), mech_writer_counter);

        if (!mNoElectricsOutput)
        {
            // the writer inside monodomain problem uses the printing timestep
            // inside HeartConfig to estimate total number of timesteps, so make
            // sure this is set to what we will use.
            HeartConfig::Instance()->SetPrintingTimeStep(mpProblemDefinition->GetMechanicsSolveTimestep());

            // When we consider archiving these simulations we will need to get a bool back from the following
            // command to decide whether or not the file is being extended, and hence whether to write the
            // initial conditions into the .h5 file.
            mpElectricsProblem->InitialiseWriter();
            mpElectricsProblem->WriteOneStep(stepper.GetTime(), initial_voltage);
        }

        if (mIsWatchedLocation)
        {
            WriteWatchedLocationData(stepper.GetTime(), initial_voltage);
        }

        if (mNumTimestepsToOutputDeformationGradientsAndStress!=UNSIGNED_UNSET)
        {
            mpMechanicsSolver->WriteCurrentStrains(DEFORMATION_GRADIENT_F,"deformation_gradient",mech_writer_counter);
            mpMechanicsSolver->WriteCurrentAverageElementStresses("second_PK",mech_writer_counter);
        }
    }


    PrepareForSolve();

//// For attempting to improve Newton convergence by quadratically extrapolating from
//// last two solutions to guess next solution. Needs further investigation - appears
//// to improve convergence but lead to decreased robustness
//    std::vector<double> current_solution_previous_time_step = mpMechanicsSolver->rGetCurrentSolution();
//    std::vector<double> current_solution_second_last_time_step = mpMechanicsSolver->rGetCurrentSolution();
//    bool first_step = true;


    // reset this to false, may be reset again below
    mpMechanicsSolver->SetComputeAverageStressPerElementDuringSolve(false);

    while (!stepper.IsTimeAtEnd())
    {
        LOG(2, "\nCurrent time = " << stepper.GetTime());
        // LCOV_EXCL_START
        if (verbose_during_solve)
        {
            // also output time to screen as newton solve information will be output
            std::cout << "\n\n ** Current time = " << stepper.GetTime() << "\n";
        }
        // LCOV_EXCL_STOP

        /////////////////////////////////////////////////////////////////////////////////////
        ////
        ////  Pass the current deformation information back to the electro-physiology
        ////  solver (MEF)
        ////
        //////////////////////////////////////////////////////////////////////////////////////
        if (mpProblemDefinition->GetDeformationAffectsCellModels() || mpProblemDefinition->GetDeformationAffectsConductivity())
        {
            //  Determine the stretch and deformation gradient on each element.
            //
            //  If mpProblemDefinition->GetDeformationAffectsCellModels()==true:
            //  Stretch will be passed to the cell models.
            //
            //  If mpProblemDefinition->GetDeformationAffectsConductivity()==true:
            //  The deformation gradient needs to be set up but does not need to be passed to the tissue
            //  so that F is used to compute the conductivity. Instead this is
            //  done through the line "mpElectricsProblem->GetMonodomainTissue()->SetConductivityModifier(this);" line above, which means
            //  rCalculateModifiedConductivityTensor() will be called on this class by the tissue, which then uses the F

            mpCardiacMechSolver->ComputeDeformationGradientAndStretchInEachElement(mDeformationGradientsForEachMechanicsElement, mStretchesForEachMechanicsElement);
        }

        if (mpProblemDefinition->GetDeformationAffectsCellModels())
        {
            //  Set the stretches on each of the cell models
            for (unsigned global_index = mpElectricsMesh->GetDistributedVectorFactory()->GetLow();
                         global_index < mpElectricsMesh->GetDistributedVectorFactory()->GetHigh();
                         global_index++)
            {
                unsigned containing_elem = mpMeshPair->rGetCoarseElementsForFineNodes()[global_index];
                double stretch = mStretchesForEachMechanicsElement[containing_elem];
                mpElectricsProblem->GetTissue()->GetCardiacCell(global_index)->SetStretch(stretch);
            }
        }

        p_electrics_solver->SetTimeStep(HeartConfig::Instance()->GetPdeTimeStep());

        /////////////////////////////////////////////////////////////////////////
        ////
        ////  Solve for the electrical activity
        ////
        /////////////////////////////////////////////////////////////////////////
        LOG(2, "  Solving electrics");
        MechanicsEventHandler::BeginEvent(MechanicsEventHandler::NON_MECH);
        for (unsigned i=0; i<mNumElecTimestepsPerMechTimestep; i++)
        {
            double current_time = stepper.GetTime() + i*HeartConfig::Instance()->GetPdeTimeStep();
            double next_time = stepper.GetTime() + (i+1)*HeartConfig::Instance()->GetPdeTimeStep();

            // solve the electrics
            p_electrics_solver->SetTimes(current_time, next_time);
            p_electrics_solver->SetInitialCondition( initial_voltage );

            electrics_solution = p_electrics_solver->Solve();

            PetscReal min_voltage, max_voltage;
            VecMax(electrics_solution,PETSC_NULL,&max_voltage); //the second param is where the index would be returned
            VecMin(electrics_solution,PETSC_NULL,&min_voltage);
            if (i==0)
            {
                LOG(2, "  minimum and maximum voltage is " << min_voltage <<", "<<max_voltage);
            }
            else if (i==1)
            {
                LOG(2, "  ..");
            }

            PetscTools::Destroy(initial_voltage);
            initial_voltage = electrics_solution;
        }

        if (mpProblemDefinition->GetDeformationAffectsConductivity())
        {
            p_electrics_solver->SetMatrixIsNotAssembled();
        }


        /////////////////////////////////////////////////////////////////////////
        ////
        ////  Pass the electrical information to the contraction models:
        ////
        /////////////////////////////////////////////////////////////////////////

        // compute Ca_I at each quad point (by interpolation, using the info on which
        // electrics element the quad point is in. Then set Ca_I on the mechanics solver
        LOG(2, "  Interpolating Ca_I and voltage");

        //Collect the distributed calcium data into one Vec to be later replicated
        for (unsigned node_index = 0; node_index<mpElectricsMesh->GetNumNodes(); node_index++)
        {
            if (mpElectricsMesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(node_index))
            {
                double calcium_value = mpElectricsProblem->GetTissue()->GetCardiacCell(node_index)->GetIntracellularCalciumConcentration();
                VecSetValue(calcium_data, node_index ,calcium_value, INSERT_VALUES);
            }
        }
        PetscTools::Barrier();//not sure this is needed

        //Replicate electrics solution and calcium (replication is inside this constructor of ReplicatableVector)
        ReplicatableVector electrics_solution_repl(electrics_solution);//size=(number of electrics nodes)*ELEC_PROB_DIM
        ReplicatableVector calcium_repl(calcium_data);//size = number of electrics nodes

        //interpolate values onto mechanics mesh
        for (unsigned i=0; i<mpMeshPair->rGetElementsAndWeights().size(); i++)
        {
            double interpolated_CaI = 0;
            double interpolated_voltage = 0;

            Element<DIM,DIM>& element = *(mpElectricsMesh->GetElement(mpMeshPair->rGetElementsAndWeights()[i].ElementNum));

            for (unsigned node_index = 0; node_index<element.GetNumNodes(); node_index++)
            {
                unsigned global_index = element.GetNodeGlobalIndex(node_index);
                double CaI_at_node =  calcium_repl[global_index];
                interpolated_CaI += CaI_at_node*mpMeshPair->rGetElementsAndWeights()[i].Weights(node_index);
                //the following line assumes interleaved solution for ELEC_PROB_DIM>1 (e.g, [Vm_0, phi_e_0, Vm1, phi_e_1...])
                interpolated_voltage += electrics_solution_repl[global_index*ELEC_PROB_DIM]*mpMeshPair->rGetElementsAndWeights()[i].Weights(node_index);
            }

            assert(i<mInterpolatedCalciumConcs.size());
            assert(i<mInterpolatedVoltages.size());
            mInterpolatedCalciumConcs[i] = interpolated_CaI;
            mInterpolatedVoltages[i] = interpolated_voltage;
        }

        LOG(2, "  Setting Ca_I. max value = " << Max(mInterpolatedCalciumConcs));

        // NOTE IF NHS: HERE WE SHOULD PERHAPS CHECK WHETHER THE CELL MODELS HAVE Ca_Trop
        // AND UPDATE FROM NHS TO CELL_MODEL, BUT NOT SURE HOW TO DO THIS.. (esp for implicit)

        // set [Ca], V, t
        mpCardiacMechSolver->SetCalciumAndVoltage(mInterpolatedCalciumConcs, mInterpolatedVoltages);
        MechanicsEventHandler::EndEvent(MechanicsEventHandler::NON_MECH);


        /////////////////////////////////////////////////////////////////////////
        ////
        ////  Solve for the deformation
        ////
        /////////////////////////////////////////////////////////////////////////
        LOG(2, "  Solving mechanics ");
        mpMechanicsSolver->SetWriteOutput(false);

        // make sure the mechanics solver knows the current time (in case
        // the traction say is time-dependent).
        mpMechanicsSolver->SetCurrentTime(stepper.GetTime());

        // see if we will need to output stresses at the end of this timestep
        if (mNumTimestepsToOutputDeformationGradientsAndStress!=UNSIGNED_UNSET
            && (counter+1)%mNumTimestepsToOutputDeformationGradientsAndStress == 0)
        {
            mpMechanicsSolver->SetComputeAverageStressPerElementDuringSolve(true);
        }

//// For attempting to improve Newton convergence by quadratically extrapolating from
//// last two solutions to guess next solution. See comments above
//        for (unsigned i=0; i<mpMechanicsSolver->rGetCurrentSolution().size(); i++)
//        {
//            double current = mpMechanicsSolver->rGetCurrentSolution()[i];
//            double previous = current_solution_previous_time_step[i];
//            double second_last = current_solution_second_last_time_step[i];
//            //double guess = 2*current - previous;
//            double guess = 3*current - 3*previous + second_last;
//
//            if (!first_step)
//            {
//                current_solution_second_last_time_step[i] = current_solution_previous_time_step[i];
//            }
//            current_solution_previous_time_step[i] = mpMechanicsSolver->rGetCurrentSolution()[i];
//            mpMechanicsSolver->rGetCurrentSolution()[i] = guess;
//        }
//        first_step = false;

        MechanicsEventHandler::BeginEvent(MechanicsEventHandler::ALL_MECH);
        mpCardiacMechSolver->Solve(stepper.GetTime(), stepper.GetNextTime(), mpProblemDefinition->GetContractionModelOdeTimestep());
        MechanicsEventHandler::EndEvent(MechanicsEventHandler::ALL_MECH);

        LOG(2, "    Number of newton iterations = " << mpMechanicsSolver->GetNumNewtonIterations());

         // update the current time
        stepper.AdvanceOneTimeStep();
        counter++;


        /////////////////////////////////////////////////////////////////////////
        ////
        ////  Output the results
        ////
        /////////////////////////////////////////////////////////////////////////
        MechanicsEventHandler::BeginEvent(MechanicsEventHandler::OUTPUT);
        if (mWriteOutput && (counter%WRITE_EVERY_NTH_TIME==0))
        {
            LOG(2, "  Writing output");
            // write deformed position
            mech_writer_counter++;
            mpMechanicsSolver->SetWriteOutput();
            mpMechanicsSolver->WriteCurrentSpatialSolution("solution","nodes",mech_writer_counter);

            p_cmgui_writer->WriteDeformationPositions(rGetDeformedPosition(), counter);

            if (!mNoElectricsOutput)
            {
                mpElectricsProblem->mpWriter->AdvanceAlongUnlimitedDimension();
                mpElectricsProblem->WriteOneStep(stepper.GetTime(), electrics_solution);
            }

            if (mIsWatchedLocation)
            {
                WriteWatchedLocationData(stepper.GetTime(), electrics_solution);
            }
            OnEndOfTimeStep(counter);

            if (mNumTimestepsToOutputDeformationGradientsAndStress!=UNSIGNED_UNSET && counter%mNumTimestepsToOutputDeformationGradientsAndStress==0)
            {
                mpMechanicsSolver->WriteCurrentStrains(DEFORMATION_GRADIENT_F,"deformation_gradient",mech_writer_counter);
                mpMechanicsSolver->WriteCurrentAverageElementStresses("second_PK",mech_writer_counter);
            }
            mpMechanicsSolver->SetComputeAverageStressPerElementDuringSolve(false);
        }
        MechanicsEventHandler::EndEvent(MechanicsEventHandler::OUTPUT);

        // write the total elapsed time..
        LogFile::Instance()->WriteElapsedTime("  ");
    }

    if ((mWriteOutput) && (!mNoElectricsOutput))
    {
        HeartConfig::Instance()->Reset();
        mpElectricsProblem->mpWriter->Close();
        delete mpElectricsProblem->mpWriter;

        std::string input_dir = mOutputDirectory+"/electrics";

        // Convert simulation data to meshalyzer format - commented
        std::string config_directory = HeartConfig::Instance()->GetOutputDirectory();
        HeartConfig::Instance()->SetOutputDirectory(input_dir);

//      Hdf5ToMeshalyzerConverter<DIM,DIM> meshalyzer_converter(FileFinder(input_dir, RelativeTo::ChasteTestOutput),
//                                                              "voltage", mpElectricsMesh,
//                                                                HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(),
//                                                                HeartConfig::Instance()->GetVisualizerOutputPrecision() );

        // convert output to CMGUI format
        Hdf5ToCmguiConverter<DIM,DIM> cmgui_converter(FileFinder(input_dir, RelativeTo::ChasteTestOutput),
                                                      "voltage", mpElectricsMesh, mHasBath,
                                                      HeartConfig::Instance()->GetVisualizerOutputPrecision());

        // Write mesh in a suitable form for meshalyzer
        //std::string output_directory =  mOutputDirectory + "/electrics/output";
        // Write the mesh
        //MeshalyzerMeshWriter<DIM,DIM> mesh_writer(output_directory, "mesh", false);
        //mesh_writer.WriteFilesUsingMesh(*mpElectricsMesh);
        // Write the parameters out
        //HeartConfig::Instance()->Write();

        // interpolate the electrical data onto the mechanics mesh nodes and write CMGUI...
        // Note: this calculates the data on ALL nodes of the mechanics mesh (incl internal,
        // non-vertex ones), which won't be used if linear CMGUI visualisation
        // of the mechanics solution is used.
        VoltageInterpolaterOntoMechanicsMesh<DIM> converter(*mpElectricsMesh,*mpMechanicsMesh,variable_names,input_dir,"voltage");

        // reset to the default value
        HeartConfig::Instance()->SetOutputDirectory(config_directory);
    }

    if (p_cmgui_writer)
    {
        if (mNoElectricsOutput)
        {
            p_cmgui_writer->WriteCmguiScript("","undeformed");
        }
        else
        {
            p_cmgui_writer->WriteCmguiScript("../../electrics/cmgui_output/voltage_mechanics_mesh","undeformed");
        }
        delete p_cmgui_writer;
    }
    PetscTools::Destroy(electrics_solution);
    PetscTools::Destroy(calcium_data);
    delete p_electrics_solver;

    MechanicsEventHandler::EndEvent(MechanicsEventHandler::ALL);
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
double CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::Max(std::vector<double>& vec)
{
    double max = -1e200;
    for (unsigned i=0; i<vec.size(); i++)
    {
        if (vec[i]>max) max=vec[i];
    }
    return max;
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::SetNoElectricsOutput()
{
    mNoElectricsOutput = true;
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::SetWatchedPosition(c_vector<double,DIM> watchedLocation)
{
    mIsWatchedLocation = true;
    mWatchedLocation = watchedLocation;
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
void CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::SetOutputDeformationGradientsAndStress(double timeStep)
{
    mNumTimestepsToOutputDeformationGradientsAndStress = (unsigned) floor((timeStep/mpProblemDefinition->GetMechanicsSolveTimestep())+0.5);
    if (fabs(mNumTimestepsToOutputDeformationGradientsAndStress*mpProblemDefinition->GetMechanicsSolveTimestep() - timeStep) > 1e-6)
    {
        EXCEPTION("Timestep provided for SetOutputDeformationGradientsAndStress() is not a multiple of mechanics solve timestep");
    }
}

template<unsigned DIM, unsigned ELEC_PROB_DIM>
std::vector<c_vector<double,DIM> >& CardiacElectroMechanicsProblem<DIM,ELEC_PROB_DIM>::rGetDeformedPosition()
{
    return mpMechanicsSolver->rGetDeformedPosition();
}

// Explicit instantiation

//note: 1d incompressible material doesn't make sense
template class CardiacElectroMechanicsProblem<2,1>;
template class CardiacElectroMechanicsProblem<3,1>;
template class CardiacElectroMechanicsProblem<2,2>;
template class CardiacElectroMechanicsProblem<3,2>;
