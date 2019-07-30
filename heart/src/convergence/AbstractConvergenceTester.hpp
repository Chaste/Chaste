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


#ifndef ABSTRACTCONVERGENCETESTER_HPP_
#define ABSTRACTCONVERGENCETESTER_HPP_

#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"

#include <petscvec.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include "AbstractTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "Hdf5DataReader.hpp"
#include "GeneralPlaneStimulusCellFactory.hpp"
#include "CuboidMeshConstructor.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "ConstBoundaryCondition.hpp"
#include "StimulusBoundaryCondition.hpp"

#include "Warnings.hpp"
#include "Timer.hpp"

typedef enum StimulusType_
{
    PLANE=0,
    QUARTER,
    NEUMANN
} StimulusType;

/**
 * RampedQuarterStimulusCellFactory stimulates a quarter of a mesh of width mMeshWidth
 * ie all the cells in 0 < x <= mMeshWidth/4
 */
template <class CELL, unsigned DIM>
class RampedQuarterStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    /** define a new vector of stimuli - one for each step in x-direction*/
    std::vector< boost::shared_ptr<SimpleStimulus> > mpStimuli;
    /** Width (x-width) of mesh*/
    double mMeshWidth;
    /** Step size of mesh is derived from the width and the number of elements across*/
    double mStepSize;
    /** The number of stimulated levels covering 0 <= x < mMeshWidth/4.
     * Note that the nodes on the quarter level are not included and are unstimulated
     */
    unsigned mLevels;
public:

    /**
     * Constructor.
     * @param meshWidth x-width of mesh
     * @param numElemAcross this allows us to deduce the mesh step size.
     * @param fullStim  the maximum stimulus level
     */
    RampedQuarterStimulusCellFactory(double meshWidth, unsigned numElemAcross, double fullStim=-1e7)
        : AbstractCardiacCellFactory<DIM>(),
          mMeshWidth(meshWidth),
          mStepSize(meshWidth/numElemAcross),
          mLevels(numElemAcross/4)
    {
        assert(numElemAcross%4 == 0); //numElemAcross is supposed to be a multiple of 4

        for (unsigned level=0; level<mLevels; level++)
        {
            double this_stim = fullStim - (level*fullStim)/mLevels;
            //this_stim is full_stim at the zero level and would be zero at level=mLevels
            mpStimuli.push_back((boost::shared_ptr<SimpleStimulus>)new SimpleStimulus(this_stim, 0.5));
        }
    }


    /**
     * @return a newly createdcell model
     *
     * @param pNode  Pointer to node object
     */
    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        double x = pNode->GetPoint()[0];
        double d_level = x/mStepSize;
        unsigned level = (unsigned) d_level;
        assert(fabs(level-d_level) < DBL_MAX); //x ought to really be a multiple of the step size

        if (level < mLevels)
        {
            return new CELL(this->mpSolver, this->mpStimuli[level]);
        }
        else
        {
            return new CELL(this->mpSolver, this->mpZeroStimulus);
        }
    }
};


/**
 * AbstractUntemplatedConvergenceTester
 * contains core functionality used in more convergence testers
 */
class AbstractUntemplatedConvergenceTester
{
protected:
    /** Mesh width (for cuboid mesh)*/
    double mMeshWidth;
public:
    /** OdeTimeStep to be varied in OdeConvergenceTester etc*/
    double OdeTimeStep;
    /** PdeTimeStep to be varied in PdeConvergenceTester etc*/
    double PdeTimeStep;
    /** Mesh number - mesh 0 has 4 elements in each space dimension  0.05cm on a 0.2cm mesh
     *              - mesh 1 has 8 (0.025cm)
     *              - mesh 2 has 16 (0.0125cm)
     */
    unsigned MeshNum;
    double RelativeConvergenceCriterion; /**< Main convergence test is LastDifference < RelativeConvergenceCriterion */
    double LastDifference; /**< Used to store and retrieve the difference between success runs (relative 2-norm of Vm at 3rd quarter node).*/
    double Apd90FirstQn; /**< Used to store and retrieve the APD90 of a node in the first quarter x-value*/
    double Apd90ThirdQn; /**< Used to store and retrieve the APD90 of a node in the third quarter x-value*/
    double ConductionVelocity;  /**< Used to store and retrieve the conduction velocity between the first & third quarter nodes*/
    /** Set to true once the result of a test is known (either read from file
     * or produced by the coarsest run of the tester).
     */
    bool PopulatedResult;
    /** true if converging to a known standard result
     *  used in StimulusConverger in projects/jmpf
     */
    bool FixedResult;

    /**
     * A value to be used with a plane stimulus (not scaled by space-step)
     *  used in StimulusConverger in projects/jmpf
     */
    double AbsoluteStimulus;
    bool SimulateFullActionPotential; /**< Set it true in order to run simulations for long enough to get a whole AP (and thus get convergence history for APD90).  Note that this slackens the relative L2 norm criterion, since we are in plateau phase for longer*/
    bool Converged; /**< Set to true when convergence has been reached */
    StimulusType Stimulus; /**< The type of stimulus: PLANE (x=0), REGION (first quarter in x) or NEUMANN (monodomain only)*/
    double NeumannStimulus; /**< Quantity of face stimulus to use in the Neumann case */

    AbstractUntemplatedConvergenceTester();

    /**
     * Run the same test at different levels of refinement until
     * some convergence criterion is met.
     * @param nameOfTest The name of the convergence test (typically the name in the suite) for use in naming files.
     */
    virtual void Converge(std::string nameOfTest)=0;

    virtual ~AbstractUntemplatedConvergenceTester();
};

/**
 * Use template specialization to set the appropriate conductivities for the problem.
 */
template<class CARDIAC_PROBLEM, unsigned DIM>
void SetConductivities(CARDIAC_PROBLEM& rCardiacProblem);

template<unsigned DIM>
void SetConductivities(BidomainProblem<DIM>& rProblem)
{
    c_vector<double, DIM> conductivities;
    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 1.75;
    }
    HeartConfig::Instance()->SetIntracellularConductivities(conductivities);

    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 7.0;
    }
    HeartConfig::Instance()->SetExtracellularConductivities(conductivities);
}

template<unsigned DIM>
void SetConductivities(MonodomainProblem<DIM>& rProblem)
{
    c_vector<double, DIM> conductivities;
    for (unsigned i=0; i<DIM; i++)
    {
        conductivities[i] = 1.75;
    }
    HeartConfig::Instance()->SetIntracellularConductivities(conductivities);
}

/**
 * AbstractConvergenceTester
 * Run convergence for a particular cell type, mono/bidomain and dimension
 */
template<class CELL, class CARDIAC_PROBLEM, unsigned DIM, unsigned PROBLEM_DIM>
class AbstractConvergenceTester : public AbstractUntemplatedConvergenceTester
{
public:
    /**
     * Run the same test at different levels of refinement until some convergence criterion is met.
     * @param nameOfTest The name of the convergence test (typically the name in the suite) for use in naming files.
     * \todo This is a scarily long method; could do with some parts extracted?
     */
    void Converge(std::string nameOfTest)
    {
        // Create the meshes on which the test will be based
        const std::string mesh_dir = "ConvergenceMesh";
        OutputFileHandler output_file_handler(mesh_dir);
        ReplicatableVector voltage_replicated;

        unsigned file_num=0;

        // Create a file for storing conduction velocity and AP data and write the header
        OutputFileHandler conv_info_handler("ConvergencePlots"+nameOfTest, false);
        out_stream p_conv_info_file;
        if (PetscTools::AmMaster())
        {
            std::cout << "=========================== Beginning Test...==================================\n";
            p_conv_info_file = conv_info_handler.OpenOutputFile(nameOfTest+"_info.csv");
            (*p_conv_info_file) << "#Abcisa\t"
                                << "l2-norm-full\t"
                                << "l2-norm-onset\t"
                                << "Max absolute err\t"
                                << "APD90_1st_quad\t"
                                << "APD90_3rd_quad\t"
                                << "Conduction velocity (relative diffs)" << std::endl;
        }
        SetInitialConvergenceParameters();

        double prev_apd90_first_qn=0.0;
        double prev_apd90_third_qn=0.0;
        double prev_cond_velocity=0.0;
        std::vector<double> prev_voltage;
        std::vector<double> prev_times;
        PopulateStandardResult(prev_voltage, prev_times);

        do
        {
            CuboidMeshConstructor<DIM> constructor;

            //If the printing time step is too fine, then simulations become I/O bound without much improvement in accuracy
            double printing_step = this->PdeTimeStep;
// LCOV_EXCL_START
            while (printing_step < 1.0e-4)
            {
                printing_step *= 2.0;
                std::cout<<"Warning: PrintingTimeStep increased\n";
            }
// LCOV_EXCL_STOP
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(this->OdeTimeStep, this->PdeTimeStep, printing_step);
// LCOV_EXCL_START
            if (SimulateFullActionPotential)
            {
                HeartConfig::Instance()->SetSimulationDuration(400.0);
            }
            else
            {
                HeartConfig::Instance()->SetSimulationDuration(8.0);
            }
// LCOV_EXCL_STOP
            HeartConfig::Instance()->SetOutputDirectory("Convergence"+nameOfTest);
            HeartConfig::Instance()->SetOutputFilenamePrefix("Results");

            DistributedTetrahedralMesh<DIM, DIM> mesh;
            constructor.Construct(mesh, this->MeshNum, mMeshWidth);

            unsigned num_ele_across = SmallPow(2u, this->MeshNum+2); // number of elements in each dimension

            AbstractCardiacCellFactory<DIM>* p_cell_factory=NULL;

            switch (this->Stimulus)
            {
                case NEUMANN:
                {
                    p_cell_factory = new ZeroStimulusCellFactory<CELL, DIM>();
                    break;
                }
                case PLANE:
                {
                    p_cell_factory = new GeneralPlaneStimulusCellFactory<CELL, DIM>(num_ele_across, constructor.GetWidth(), this->AbsoluteStimulus);
                    break;
                }
                case QUARTER:
                {
                    ///\todo consider reducing all stimuli to match this one.
                    p_cell_factory = new RampedQuarterStimulusCellFactory<CELL, DIM>(constructor.GetWidth(), num_ele_across, this->AbsoluteStimulus/10.0);
                    break;
                }
                default:
                    NEVER_REACHED;
            }


            CARDIAC_PROBLEM cardiac_problem(p_cell_factory);
            cardiac_problem.SetMesh(&mesh);

            // Calculate positions of nodes 1/4 and 3/4 through the mesh
            unsigned third_quadrant_node;
            unsigned first_quadrant_node;
            switch(DIM)
            {
                case 1:
                {
                    first_quadrant_node = (unsigned) (0.25*constructor.GetNumElements());
                    third_quadrant_node = (unsigned) (0.75*constructor.GetNumElements());
                    break;
                }
                case 2:
                {
                    unsigned n= SmallPow (2u, this->MeshNum+2);
                    first_quadrant_node =   (n+1)*(n/2)+  n/4 ;
                    third_quadrant_node =   (n+1)*(n/2)+3*n/4 ;
                    break;
                }
                case 3:
                {
                    const unsigned first_quadrant_nodes_3d[5]={61, 362, 2452, 17960, 137296};
                    const unsigned third_quadrant_nodes_3d[5]={63, 366, 2460, 17976, 137328};
                    assert(this->PdeTimeStep<5);
                    first_quadrant_node = first_quadrant_nodes_3d[this->MeshNum];
                    third_quadrant_node = third_quadrant_nodes_3d[this->MeshNum];
                    break;
                }

                default:
                    NEVER_REACHED;
            }

            double mesh_width=constructor.GetWidth();

            // We only need the output of these two nodes
            std::vector<unsigned> nodes_to_be_output;
            nodes_to_be_output.push_back(first_quadrant_node);
            nodes_to_be_output.push_back(third_quadrant_node);
            cardiac_problem.SetOutputNodes(nodes_to_be_output);

            // The results of the tests were originally obtained with the following conductivity
            // values. After implementing fibre orientation the defaults changed. Here we set
            // the former ones to be used.
            SetConductivities(cardiac_problem);

            cardiac_problem.Initialise();

            boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM>);
            SimpleStimulus stim(NeumannStimulus, 0.5);
            if (Stimulus==NEUMANN)
            {

                StimulusBoundaryCondition<DIM>* p_bc_stim = new StimulusBoundaryCondition<DIM>(&stim);

                // get mesh
                AbstractTetrahedralMesh<DIM, DIM> &r_mesh = cardiac_problem.rGetMesh();
                // loop over boundary elements
                typename AbstractTetrahedralMesh<DIM, DIM>::BoundaryElementIterator iter;
                iter = r_mesh.GetBoundaryElementIteratorBegin();
                while (iter != r_mesh.GetBoundaryElementIteratorEnd())
                {
                    double x = ((*iter)->CalculateCentroid())[0];
                    ///\todo remove magic number? (#1884)
                    if (x*x<=1e-10)
                    {
                        p_bcc->AddNeumannBoundaryCondition(*iter, p_bc_stim);
                    }
                    iter++;
                }
                // pass the bcc to the problem
                cardiac_problem.SetBoundaryConditionsContainer(p_bcc);
            }

            DisplayRun();
            Timer::Reset();
            //// use this to get some info printed out
            //cardiac_problem.SetWriteInfo();

            try
            {
                cardiac_problem.Solve();
            }
            catch (Exception e)
            {
                WARNING("This run threw an exception.  Check convergence results\n");
                std::cout << e.GetMessage() << std::endl;
            }

            if (PetscTools::AmMaster())
            {
                std::cout << "Time to solve = " << Timer::GetElapsedTime() << " seconds\n";
            }

            OutputFileHandler results_handler("Convergence"+nameOfTest, false);
            Hdf5DataReader results_reader = cardiac_problem.GetDataReader();

            {
                std::vector<double> transmembrane_potential = results_reader.GetVariableOverTime("V", third_quadrant_node);
                std::vector<double> time_series = results_reader.GetUnlimitedDimensionValues();

                OutputFileHandler plot_file_handler("ConvergencePlots"+nameOfTest, false);
                if (PetscTools::AmMaster())
                {
                    // Write out the time series for the node at third quadrant
                    {
                        std::stringstream plot_file_name_stream;
                        plot_file_name_stream<< nameOfTest << "_Third_quadrant_node_run_"<< file_num << ".csv";
                        out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                        for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                        {
                            (*p_plot_file) << time_series[data_point] << "\t" << transmembrane_potential[data_point] << "\n";
                        }
                        p_plot_file->close();
                    }
                    // Write time series for first quadrant node
                    {
                        std::vector<double> transmembrane_potential_1qd=results_reader.GetVariableOverTime("V", first_quadrant_node);
                        std::vector<double> time_series_1qd = results_reader.GetUnlimitedDimensionValues();
                        std::stringstream plot_file_name_stream;
                        plot_file_name_stream<< nameOfTest << "_First_quadrant_node_run_"<< file_num << ".csv";
                        out_stream p_plot_file = plot_file_handler.OpenOutputFile(plot_file_name_stream.str());
                        for (unsigned data_point = 0; data_point<time_series.size(); data_point++)
                        {
                            (*p_plot_file) << time_series_1qd[data_point] << "\t" << transmembrane_potential_1qd[data_point] << "\n";
                        }
                        p_plot_file->close();
                    }
                }
                // calculate conduction velocity and APD90 error
                PropagationPropertiesCalculator ppc(&results_reader);

                try
                {
                    // LCOV_EXCL_START
                    if (SimulateFullActionPotential)
                    {
                        Apd90FirstQn = ppc.CalculateActionPotentialDuration(90.0, first_quadrant_node);
                        Apd90ThirdQn = ppc.CalculateActionPotentialDuration(90.0, third_quadrant_node);
                    }
                    // LCOV_EXCL_STOP
                    ConductionVelocity  = ppc.CalculateConductionVelocity(first_quadrant_node,third_quadrant_node,0.5*mesh_width);
                }
                // LCOV_EXCL_START
                catch (Exception e)
                {
                    std::cout << "Warning - this run threw an exception in calculating propagation.  Check convergence results\n";
                    std::cout << e.GetMessage() << std::endl;
                }
                // LCOV_EXCL_STOP
                double cond_velocity_error = 1e10;
                double apd90_first_qn_error = 1e10;
                double apd90_third_qn_error = 1e10;

                if (this->PopulatedResult)
                {
                    if (prev_cond_velocity != 0.0)
                    {
                        cond_velocity_error = fabs(ConductionVelocity - prev_cond_velocity) / prev_cond_velocity;
                    }
// LCOV_EXCL_START
                    if (prev_apd90_first_qn != 0.0)
                    {
                        apd90_first_qn_error = fabs(Apd90FirstQn - prev_apd90_first_qn) / prev_apd90_first_qn;
                    }
                    if (prev_apd90_third_qn != 0.0)
                    {
                        apd90_third_qn_error = fabs(Apd90ThirdQn - prev_apd90_third_qn) / prev_apd90_third_qn;
                    }
                    if (apd90_first_qn_error == 0.0)
                    {
                        apd90_first_qn_error = DBL_EPSILON; //Avoid log zero on plot
                    }
                    if (apd90_third_qn_error == 0.0)
                    {
                        apd90_third_qn_error = DBL_EPSILON; //Avoid log zero on plot
                    }
// LCOV_EXCL_STOP
                    if (cond_velocity_error == 0.0)
                    {
                        cond_velocity_error = DBL_EPSILON; //Avoid log zero on plot
                    }
                }

                prev_cond_velocity = ConductionVelocity;
                prev_apd90_first_qn = Apd90FirstQn;
                prev_apd90_third_qn = Apd90ThirdQn;

                // calculate l2norm
                double max_abs_error = 0;
                double sum_sq_abs_error =0;
                double sum_sq_prev_voltage = 0;
                double sum_sq_abs_error_full =0;
                double sum_sq_prev_voltage_full = 0;
                if (this->PopulatedResult)
                {
                    //If the PDE step is varying then we'll have twice as much data now as we use to have

                    unsigned time_factor=(time_series.size()-1) / (prev_times.size()-1);
                    assert (time_factor == 1 || time_factor == 2 || time_factor == 8);
                    //Iterate over the shorter time series data
                    for (unsigned data_point = 0; data_point<prev_times.size(); data_point++)
                    {
                        unsigned this_data_point=time_factor*data_point;

                        assert(time_series[this_data_point] == prev_times[data_point]);
                        double abs_error = fabs(transmembrane_potential[this_data_point]-prev_voltage[data_point]);
                        max_abs_error = (abs_error > max_abs_error) ? abs_error : max_abs_error;
                        //Only do resolve the upstroke...
                        sum_sq_abs_error_full += abs_error*abs_error;
                        sum_sq_prev_voltage_full += prev_voltage[data_point] * prev_voltage[data_point];
                        if (time_series[this_data_point] <= 8.0)
                        {
                            //In most regular cases we always do this, since the simulation stops at ms
                            sum_sq_abs_error += abs_error*abs_error;
                            sum_sq_prev_voltage += prev_voltage[data_point] * prev_voltage[data_point];
                        }
                    }

                }
                if (!this->PopulatedResult || !FixedResult)
                {
                    prev_voltage = transmembrane_potential;
                    prev_times = time_series;
                }

                if (this->PopulatedResult)
                {
                    double l2_norm_upstroke = sqrt(sum_sq_abs_error/sum_sq_prev_voltage);
                    double l2_norm_full = sqrt(sum_sq_abs_error_full/sum_sq_prev_voltage_full);

                    if (PetscTools::AmMaster())
                    {
                        (*p_conv_info_file) << std::setprecision(8)
                                            << Abscissa() << "\t"
                                            << l2_norm_full << "\t"
                                            << l2_norm_upstroke << "\t"
                                            << max_abs_error << "\t"
                                            << Apd90FirstQn <<" "<< apd90_first_qn_error <<""<< "\t"
                                            << Apd90ThirdQn <<" "<< apd90_third_qn_error <<""<< "\t"
                                            << ConductionVelocity <<" "<< cond_velocity_error  <<""<< std::endl;
                    }
                    // convergence criterion
                    this->Converged = l2_norm_full < this->RelativeConvergenceCriterion;
                    this->LastDifference = l2_norm_full;
// LCOV_EXCL_START
                    assert (time_series.size() != 1u);
                    if (time_series.back() == 0.0)
                    {
                        std::cout << "Failed after successful convergence - give up this convergence test\n";
                        break;
                    }
// LCOV_EXCL_STOP

                }

                if (time_series.back() != 0.0)
                {
                    // Simulation ran to completion
                    this->PopulatedResult=true;

                }
            }

            // Get ready for the next test by halving the time step
            if (!this->Converged)
            {
                UpdateConvergenceParameters();
                file_num++;
            }
            delete p_cell_factory;
        }
        while (!GiveUpConvergence() && !this->Converged);


        if (PetscTools::AmMaster())
        {
            p_conv_info_file->close();

            std::cout << "Results: " << std::endl;
            FileFinder info_finder = conv_info_handler.FindFile(nameOfTest + "_info.csv");
            std::ifstream info_file(info_finder.GetAbsolutePath().c_str());
            if (info_file)
            {
                std::cout << info_file.rdbuf();
                info_file.close();
            }
        }
    }

    /**
     * Show convergence details relevant to this run in std::cout
     */
    void DisplayRun()
    {
        if (!PetscTools::AmMaster())
        {
            return; //Only master displays this
        }
        unsigned num_ele_across = SmallPow(2u, this->MeshNum+2);// number of elements in each dimension
        double scaling = mMeshWidth/(double) num_ele_across;

        std::cout<<"================================================================================"<<std::endl;
        std::cout<<"Solving in " << DIM << "D\t";
        std::cout<<"Space step " << scaling << " cm (mesh " << this->MeshNum << ")" << "\n";
        std::cout<<"PDE step " << this->PdeTimeStep << " ms" << "\t";
        std::cout<<"ODE step " << this->OdeTimeStep << " ms" << "\t";
        if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
        {
            std::cout<<"KSP absolute "<<HeartConfig::Instance()->GetAbsoluteTolerance()<<"\t";
        }
        else
        {
            std::cout<<"KSP relative "<<HeartConfig::Instance()->GetRelativeTolerance()<<"\t";
        }
        switch (this->Stimulus)
        {
            case PLANE:
            std::cout<<"Stimulus = Plane\n";
            break;

            case QUARTER:
            std::cout<<"Stimulus = Ramped first quarter\n";
            break;

            case NEUMANN:
            std::cout<<"Stimulus = Neumann\n";
            break;

        }
        // Keep track of what Nightly things are doing
        std::time_t rawtime;
        std::time( &rawtime );
        std::cout << std::ctime(&rawtime);
        std::cout << std::flush;
        //HeartEventHandler::Headings();
        //HeartEventHandler::Report();

    }

public:
    virtual ~AbstractConvergenceTester() {}

    /** Initial values of parameters at the beginning of the convergence test (the parameter to be varied will be larger than the expected value at convergence)*/
    virtual void SetInitialConvergenceParameters()=0;
    /** Update the parameter which is being varied*/
    virtual void UpdateConvergenceParameters()=0;
    /** @return whether to abort the convergence test (convergence is unlikely to happen).*/
    virtual bool GiveUpConvergence()=0;
    /** @return value of the parameter which is being varied*/
    virtual double Abscissa()=0;
    /** This is currently used as stub for convergence testers which need to converge
     * to a known standardised result (the StimulusConvergence tester in projects/jmpf).
     *
     * @param result a standard vector to be sized and filled with V_m values by this method (in subclass)
     * @param times a standard vector to be sized and filled with times values by this method (in subclass)
     */
    virtual void PopulateStandardResult(std::vector<double> &result, std::vector<double> &times)
    {
        assert(this->PopulatedResult==false);
        assert(result.size()==0);
        assert(times.size()==0);
    }

    /**
     * @return when the convergence criterion is met
     */
    bool IsConverged()
    {
        return Converged;
    }

    /**
     * @param meshWidth set the dimension of the cuboid mesh (default value is 0.2cm)
     */
    void SetMeshWidth(double meshWidth)
    {
        mMeshWidth=meshWidth;
    }
};

#endif /*ABSTRACTCONVERGENCETESTER_HPP_*/
