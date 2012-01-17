/*

Copyright (C) University of Oxford, 2005-2012

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


#ifndef EXTENDEDBIDOMAINPROBLEM_HPP_
#define EXTENDEDBIDOMAINPROBLEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCardiacProblem.hpp"
#include "AbstractCardiacTissue.hpp"
#include "AbstractExtendedBidomainSolver.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "Electrodes.hpp"
#include "ExtendedBidomainTissue.hpp"

#include "HeartRegionCodes.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractStimulusFactory.hpp"
#include "ElectrodesStimulusFactory.hpp"

/**
 * Class which specifies and solves an extended bidomain problem.
 *
 * See Buist ML, Poh YC. An Extended Bidomain Framework Incorporating Multiple Cell Types. Biophysical Journal, Volume 99, Issue 1, 13-18, 7 July 2010.
 *
 * Briefly, the problems consits of 3 equations, 2 parabolic and one elliptic. The space is divided in 3 compartments:
 * - First cell
 * - Second cell
 * - Extracellular space.
 *
 * The three unkowns are the the intracellular potential of the two cells and the extracellular potential.
 * This class allows the user to specify the parameters specific for the these simulations and also different intracellular conductivities for the two cells.
 *
 * The solution vector used for calculations is arranged as a striped vector with this order:
 *
 *  - Intracellular potential of the first cell
 *  - Intracellular potential of the second cell
 *  - Extracellular potential
 *
 * However, the OUTPUT is rearranged as follows. The rearrangament occurs AFTER the solution is calculated:
 *
 *  - Transmembrane potential of the first cell (Intracellular potential of the first cell - Extracellular potential)
 *  - Transmembrane potential of the second cell (Intracellular potential of the second cell - Extracellular potential)
 *  - Extracellular potential
 *
 *  Unlike a bidomain problem, a node-wise extracellular stimulus can be set up in extended bidomain problems in absence of a bath.
 *  This is done by setting a stimulus factory via the method SetExtracellularStimulusFactory.
 *  See documentation of AbstractStimulusFactory and ElectrodesStimulusFactory for more details.
 *  Note that compatibility conditions will be taken care of by this class by calling specific methods within the stimulus factory class.
 *
 */
template<unsigned DIM>
class ExtendedBidomainProblem : public AbstractCardiacProblem<DIM,DIM, 3>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;

    friend class TestArchivingExtendedBidomain;//for testing

    /**
     * Save the member variables to an archive.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractCardiacProblem<DIM, DIM, 3> >(*this);
        archive & mpExtendedBidomainTissue;
        archive & mFixedExtracellularPotentialNodes;
        //archive & mIntracellularConductivitiesSecondCell; not allowed in some versions of boost
        archive & mVariablesIDs;
        archive & mUserSpecifiedSecondCellConductivities;
        archive & mUserHasSetBidomainValuesExplicitly;
        archive & mAmFirstCell;
        archive & mAmSecondCell;
        archive & mAmGap;
        archive & mCmFirstCell;
        archive & mCmSecondCell;
        archive & mGGap;
        archive & mGgapHeterogeneityRegions;
        archive & mGgapHeterogenousValues;
        archive & mRowForAverageOfPhiZeroed;
        archive & mApplyAveragePhieZeroConstraintAfterSolving;
        archive & mUserSuppliedExtracellularStimulus;
        archive & mHasBath;
        //archive & mpSolver;

        //archive the values for the conductivies of the second cell
        for (unsigned i = 0; i < DIM; i++)
        {
            double conductivity = mIntracellularConductivitiesSecondCell(i);
            archive & conductivity;
        }

        bool has_solution = (this->mSolution != NULL);
        archive & has_solution;
        if (has_solution)
        {
            std::string filename = ArchiveLocationInfo::GetArchiveDirectory() + "AbstractCardiacProblem_mSolution.vec";

            Hdf5DataWriter writer(*this->mpMesh->GetDistributedVectorFactory(), ArchiveLocationInfo::GetArchiveRelativePath(), "AbstractCardiacProblem_mSolution", false);
            writer.DefineFixedDimension(this->mpMesh->GetDistributedVectorFactory()->GetProblemSize());
            writer.DefineUnlimitedDimension("Time", "msec", 1);

            // Make sure the file does not take more disc space than really needed (#1200)
            writer.SetFixedChunkSize(1);

            ///\todo #1369
            assert(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering() == false );

            int V = writer.DefineVariable("V","mV");
            int V_2 = writer.DefineVariable("V_2","mV");
            int phie = writer.DefineVariable("Phi_e","mV");
            std::vector<int> variable_ids;
            variable_ids.push_back(V);
            variable_ids.push_back(V_2);
            variable_ids.push_back(phie);
            writer.EndDefineMode();
            writer.PutUnlimitedVariable(0.0);

            //re-arrange to write out voltages...
            Vec voltages_to_be_written =  this->mpMesh->GetDistributedVectorFactory()->CreateVec(3);
            DistributedVector wrapped_voltages_to_be_written = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(voltages_to_be_written);

            DistributedVector distr_solution = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(this->mSolution);
            DistributedVector::Stripe phi_i_first_cell_stripe(distr_solution,0);
            DistributedVector::Stripe phi_i_second_cell_stripe(distr_solution,1);
            DistributedVector::Stripe phi_e_stripe(distr_solution,2);


            DistributedVector::Stripe wrapped_voltages_to_be_written_first_stripe(wrapped_voltages_to_be_written,0);
            DistributedVector::Stripe wrapped_voltages_to_be_written_second_stripe(wrapped_voltages_to_be_written,1);
            DistributedVector::Stripe wrapped_voltages_to_be_written_third_stripe(wrapped_voltages_to_be_written,2);

            for (DistributedVector::Iterator index = distr_solution.Begin();
                 index != distr_solution.End();
                 ++index)
            {
                wrapped_voltages_to_be_written_first_stripe[index] = phi_i_first_cell_stripe[index] - phi_e_stripe[index];
                wrapped_voltages_to_be_written_second_stripe[index] = phi_i_second_cell_stripe[index] - phi_e_stripe[index];
                wrapped_voltages_to_be_written_third_stripe[index] = phi_e_stripe[index];
            }
            distr_solution.Restore();
            wrapped_voltages_to_be_written.Restore();

            writer.PutStripedVector(variable_ids, voltages_to_be_written);
            PetscTools::Destroy(voltages_to_be_written);
        }
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
        archive & boost::serialization::base_object<AbstractCardiacProblem<DIM, DIM, 3> >(*this);
        archive & mpExtendedBidomainTissue;
        archive & mFixedExtracellularPotentialNodes;
        //archive & mIntracellularConductivitiesSecondCell; not allowed in some versions of boost
        archive & mVariablesIDs;
        archive & mUserSpecifiedSecondCellConductivities;
        archive & mUserHasSetBidomainValuesExplicitly;
        archive & mAmFirstCell;
        archive & mAmSecondCell;
        archive & mAmGap;
        archive & mCmFirstCell;
        archive & mCmSecondCell;
        archive & mGGap;
        archive & mGgapHeterogeneityRegions;
        archive & mGgapHeterogenousValues;
        archive & mRowForAverageOfPhiZeroed;
        archive & mApplyAveragePhieZeroConstraintAfterSolving;
        archive & mUserSuppliedExtracellularStimulus;
        archive & mHasBath;
        //archive & mpSolver;

        //load the values for the conductivies of the second cell
        for (unsigned i = 0; i < DIM; i++)
        {
            double conductivity;
            archive & conductivity;
            mIntracellularConductivitiesSecondCell(i) = conductivity;
        }

        bool has_solution;
        archive & has_solution;

        if (has_solution)
        {
            std::string filename = ArchiveLocationInfo::GetArchiveDirectory() + "AbstractCardiacProblem_mSolution.vec";

            this->mSolution = this->mpMesh->GetDistributedVectorFactory()->CreateVec(3);
            DistributedVector mSolution_distri = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(this->mSolution);

            std::string archive_dir = ArchiveLocationInfo::GetArchiveRelativePath();
            Hdf5DataReader reader(archive_dir, "AbstractCardiacProblem_mSolution", !FileFinder::IsAbsolutePath(archive_dir));

            Vec V = this->mpMesh->GetDistributedVectorFactory()->CreateVec();
            Vec V_2 = this->mpMesh->GetDistributedVectorFactory()->CreateVec();
            Vec phie = this->mpMesh->GetDistributedVectorFactory()->CreateVec();

            reader.GetVariableOverNodes(V, "V", 0);
            reader.GetVariableOverNodes(V_2, "V_2", 0);
            reader.GetVariableOverNodes(phie, "Phi_e", 0);

            //from transmembrane voltages back to phi_i now...
            DistributedVector vm_first_cell_distri = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(V);
            DistributedVector vm_second_cell_distri = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(V_2);
            DistributedVector phie_distri = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(phie);

            DistributedVector::Stripe mSolution_phi_1(mSolution_distri,0);
            DistributedVector::Stripe mSolution_phi_2(mSolution_distri,1);
            DistributedVector::Stripe mSolution_phie(mSolution_distri,2);

            for (DistributedVector::Iterator index = mSolution_distri.Begin();
                 index != mSolution_distri.End();
                 ++index)
            {
                mSolution_phi_1[index] = vm_first_cell_distri[index] + phie_distri[index];//phi_i = Vm + phi_e
                mSolution_phi_2[index] = vm_second_cell_distri[index] + phie_distri[index];
                mSolution_phie[index] = phie_distri[index];
            }
            PetscTools::Destroy(V);
            PetscTools::Destroy(V_2);
            PetscTools::Destroy(phie);

            mSolution_distri.Restore();
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()


protected:

    /** The cell factory creates the cells for each node */
    AbstractCardiacCellFactory<DIM,DIM>* mpSecondCellFactory;

    /** The bidomain tissue object for the extended problem*/
    ExtendedBidomainTissue<DIM>* mpExtendedBidomainTissue;

    /** Nodes at which the extracellular voltage is fixed to zero (replicated) */
    std::vector<unsigned> mFixedExtracellularPotentialNodes;

    /**stores the values of the conductivities for the second cell.*/
    c_vector<double, DIM>  mIntracellularConductivitiesSecondCell;

    /** Used by the writer, transmembrane potential of the first cell*/
    unsigned mVoltageColumnId_Vm1;
    /** Used by the writer, transmembrane potential of the second cell*/
    unsigned mVoltageColumnId_Vm2;
    /** Used by the writer, extracellular potential*/
    unsigned mVoltageColumnId_Phie;
    /** Used by the writer, stores the variable output names*/
    std::vector<signed int> mVariablesIDs;

    /** Flag for checking that the user specified conductivities for the second cell. Get method available for this variable.*/
    bool mUserSpecifiedSecondCellConductivities;

    /** flag used to check whwther the user wanted specific, out-of-heartconfig  values*/
    bool mUserHasSetBidomainValuesExplicitly;

    /**Am for the first cell, set by the user and, if so,  set into the PDE*/
    double mAmFirstCell;
    /**Am for the second cell, set by the user and, if so,  set into the PDE*/
    double mAmSecondCell;
    /**Am for the gap junctions, set by the user and, if so,  set into the PDE*/
    double mAmGap;
    /**Cm for the first cell, set by the user and, if so,  set into the PDE*/
    double mCmFirstCell;
    /**Cm for the second cell, set by the user and, if so,  set into the PDE*/
    double mCmSecondCell;
    /**Conductance, in mS of the gap junction conductance, set by the user and, if so,  set into the PDE (otherwise its default value is 0 */
    double mGGap;

    /**
     * Vector of Ggap heterogeneity regions. Set by the user, it is passed on to the tissue object.
     * If empty (i.e., the user did not specify any heterogeneities) then the empty vector will still be passed to the tissue object, which,
     * seeing the empty vector, will apply mGgap everywhere
     */
    std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > > mGgapHeterogeneityRegions;

    /**Corresponding vector of values of Ggap for every heterogeneous region in mGgapHeterogeneityRegions*/
    std::vector<double> mGgapHeterogenousValues;


    /**Factory to generate the stimulus*/
    AbstractStimulusFactory<DIM>* mpExtracellularStimulusFactory;

    /**
     * Another method of resolving the singularity in the bidomain equations.
     * Specifies a row of the matrix at which to impose an extra condition.
     */
    int mRowForAverageOfPhiZeroed;

    /**
     * An alternative method for constraining the solution and resolving singularity.
     * Singular system is solved and THEN the constraint is applied on the phi_e (after calculating.
     *
     * This variable is checked within the specialization of "OnEndOfTimestep" method.
     * False by default as initialised in the constructor.
     */
    unsigned mApplyAveragePhieZeroConstraintAfterSolving;

    /** keeps track of whether the user set an extracellular stimulus. True if the user did.*/
    bool mUserSuppliedExtracellularStimulus;

    /** Whether the mesh has a bath, ie whether this is a bath simulation */
    bool mHasBath;

    /**
     *  Create normal initial condition but overwrite V to zero for bath nodes, if
     *  there are any.
     */
    Vec CreateInitialCondition();

    /**
     * Annotate bath nodes with the correct region code, if a bath is present.
     * Will throw if #mHasBath is set but no bath is present in the mesh.
     */
    void AnalyseMeshForBath();

    /**
     * Small helper tidy-up method that incorporates everything
     * related to the initialisation of the extracellular stimulus.
     *
     * It checks whether there is any extracellular stimulus, If not, it creates a one
     * (with default implementation to zero stimulus).
     * In any case, it passes the mesh into the stimulus, calls the SetCompatibleExtracellularStimulus method
     * and also checks if there are any grounded nodes. If so, it calls  SetFixedExtracellularPotentialNodes.
     */
    void ProcessExtracellularStimulus();

    /**
     * We need to save the assembler that is being used to be able to switch
     * off the electrodes (by adding default boundary conditions to the
     * assembler)
     */
    AbstractExtendedBidomainSolver<DIM,DIM>* mpSolver;

    /** Create our bidomain PDE object */
    virtual AbstractCardiacTissue<DIM> *CreateCardiacTissue();

    /** Create a suitable bidomain solver */
    virtual AbstractDynamicLinearPdeSolver<DIM,DIM,3>* CreateSolver();

public:
    /**
     * Constructor
     * @param pCellFactory User defined cell factory which shows how the pde should
     *   create cells.
     *  @param pSecondCellFactory User defined cell factory which shows how the pde should
     *   create cells.
     *  @param hasBath Whether the simulation has a bath (if this is true, all elements with
     *   attribute = 1 will be set to be bath elements (the rest should have
     *   attribute = 0)).
     */
    ExtendedBidomainProblem(AbstractCardiacCellFactory<DIM>* pCellFactory, AbstractCardiacCellFactory<DIM>* pSecondCellFactory, bool hasBath = false);

    /**
     * Archiving constructor
     */
    ExtendedBidomainProblem();

    /**
     * Destructor
     */
    ~ExtendedBidomainProblem();

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
     * Allow the user to specify different values for the conductivities of the second cell type
     *
     * @param conductivities : the intracellular conductivities of the second cell that you wish to set
     */
    void SetIntracellularConductivitiesForSecondCell(c_vector<double, DIM> conductivities);

    /**
     * Allow the user to specify different values for the bidomain parameters in the extened bidomain framework.
     * This method switches the mUserHasSetBidomainValuesExplicitly to true and tells the problem to overwrite the corresponding parameters in HeartConfig.
     * It needs to be called before Initialise() to have any effect.
     *
     * @param Am1 value of the surface-to-volume ratio for the first cell
     * @param Am2 value of the surface-to-volume ratio for the second cell
     * @param AmGap value of the surface-to-volume ratio for the gap junction
     * @param Cm1 value of the capacitance  for the first cell
     * @param Cm2 value of the capacitance for the second cell
     * @param Ggap value, in mS of the gap junction conductance
     */
    void SetExtendedBidomainParameters(double Am1, double Am2, double AmGap, double Cm1, double Cm2, double Ggap);

    /**
     *  Set the values of mCellHeterogeneityRegions and mGgapValues for the heterogeneities of Ggap.
     *  It just sets the member variables here that will later be passed on to the tissue object.
     *  It also checks that the two have the same size. Throws otherwise.
     *
     *  @param rGgapHeterogeneityRegions a vector of heterogeneity regions for gap junctions
     *  @param rGgapValues a vector (of the same size as rGgapHeterogeneityRegions) with the respective values of Ggap for every region.
     */
    void SetGgapHeterogeneities ( std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rGgapHeterogeneityRegions, std::vector<double>& rGgapValues);

    /**
     * Sets a stimulus factory as the extracellular one.
     *
     * @param pFactory the factory to be set.
     */
    void SetExtracellularStimulusFactory( AbstractStimulusFactory<DIM>* pFactory);


    /**
     * Set which row of the linear system should be used to enforce the
     * condition that the average of phi_e is zero.  If not called, this
     * condition will not be used.
     *
     * @param node  the mesh node index giving the row at which to impose the constraint
     */
    void SetNodeForAverageOfPhiZeroed(unsigned node);


    /**
     *  Access the tissue object. Can only be called after Initialise()
     */
    ExtendedBidomainTissue<DIM>* GetExtendedBidomainTissue();

    /**
     *  Print out time and max/min of the intracellular potential of the two cells and phi_e values at current time.
     *
     *  @param time  current time.
     */
    void WriteInfo(double time);

    /**
     * Define what variables are written to the primary results file.
     * Hardcoded variable names are "V", "V_2" and "Phi_e" for the three variables.
     * If you request any extra variable (via HeartConfig), this method will also define
     * the extra variables by calling a method in the parent class.
     *
     * @param extending  whether we are extending an existing results file
     */
    virtual void DefineWriterColumns(bool extending);

    /**
     * Write one timestep of output data to the primary results file.
     * The solution of the problem is in term of Phi_i_1, Phi_i_2 and Phi_e
     * (i.e., intracellular potential of the two cells plus extracellular potential).
     * This method works out the transmembrane potential of the two cells (V and V_2)
     * calculated as V = Phi_i_1 - Phi_e and V_2 = Phi_i_2-Phi_e.
     * It then writes to file V, V_2 and Phi_e.
     *
     * It also writes any extra variable (defined by HeartConfig).
     * Note that it does the job by calling the method in the parent class.
     * This implies that the extra variable must be in the first cell (not the second)
     * because the generic method to write extra variables only looks into the first cell factory.
     * TODO: write a specific method for extended problems that looks in both cell factories
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
     * Return whether this is a problem with bath or not
     */
    bool GetHasBath();

    /**
     * Set whether there is a bath or not.
     * Useful for covering exceptions only (for example if bath problems are not supported and there is no method to analyse the mesh for a bath).
     *
     * @param hasBath whether we want the problem to have bath or not.
     */
    void SetHasBath(bool hasBath);

};

#include "SerializationExportWrapper.hpp" // Must be last
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtendedBidomainProblem)

#endif /*EXTENDEDBIDOMAINPROBLEM_HPP_*/
