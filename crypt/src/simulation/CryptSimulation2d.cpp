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

#include "CryptSimulation2d.hpp"
#include "WntConcentration.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "SmartPointers.hpp"

CryptSimulation2d::CryptSimulation2d(AbstractCellPopulation<2>& rCellPopulation,
                                     bool deleteCellPopulationInDestructor,
                                     bool initialiseCells)
    : OffLatticeSimulation<2>(rCellPopulation,
                             deleteCellPopulationInDestructor,
                             initialiseCells),
      mWriteBetaCatenin(false),
      mUsingMeshBasedCellPopulation(false)
{
    /* Throw an exception message if not using a  MeshBasedCellPopulation or a VertexBasedCellPopulation.
     * This is to catch NodeBasedCellPopulations as AbstactOnLatticeBasedCellPopulations are caught in
     * the OffLatticeSimulation constructor.
     */
    if ( (dynamic_cast<VertexBasedCellPopulation<2>*>(&rCellPopulation) == NULL)
         &&(dynamic_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation) == NULL) )
    {
        EXCEPTION("CryptSimulation2d is to be used with MeshBasedCellPopulation or VertexBasedCellPopulation (or subclasses) only");
    }

    if (dynamic_cast<MeshBasedCellPopulation<2>*>(&mrCellPopulation))
    {
        mUsingMeshBasedCellPopulation = true;
    }

    /*
     * To check if beta-catenin results will be written to file, we test if the first
     * cell has a cell-cycle model that is a subclass of AbstractVanLeeuwen2009WntSwatCellCycleModel.
     * In doing so, we assume that all cells in the simulation have the same cell cycle
     * model.
     */
    if (dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(mrCellPopulation.Begin()->GetCellCycleModel()))
    {
        mWriteBetaCatenin = true;
    }

    if (!mDeleteCellPopulationInDestructor)
    {
        // Pass a CryptSimulationBoundaryCondition object into mBoundaryConditions
        MAKE_PTR_ARGS(CryptSimulationBoundaryCondition<2>, p_bc, (&rCellPopulation));
        AddCellPopulationBoundaryCondition(p_bc);
    }
}

CryptSimulation2d::~CryptSimulation2d()
{
}

c_vector<double, 2> CryptSimulation2d::CalculateCellDivisionVector(CellPtr pParentCell)
{
    if (mUsingMeshBasedCellPopulation)
    {
        // Location of parent and daughter cells
        c_vector<double, 2> parent_coords = mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, 2> daughter_coords;

        // Get separation parameter
        double separation =
            static_cast<MeshBasedCellPopulation<2>*>(&mrCellPopulation)->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, 2> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        double random_angle = RandomNumberGenerator::Instance()->ranf();
        random_angle *= 2.0*M_PI;

        random_vector(0) = 0.5*separation*cos(random_angle);
        random_vector(1) = 0.5*separation*sin(random_angle);

        c_vector<double, 2> proposed_new_parent_coords = parent_coords - random_vector;
        c_vector<double, 2> proposed_new_daughter_coords = parent_coords + random_vector;

        if ((proposed_new_parent_coords(1) >= 0.0) && (proposed_new_daughter_coords(1) >= 0.0))
        {
            // We are not too close to the bottom of the cell population, so move parent
            parent_coords = proposed_new_parent_coords;
            daughter_coords = proposed_new_daughter_coords;
        }
        else
        {
            proposed_new_daughter_coords = parent_coords + 2.0*random_vector;
            while (proposed_new_daughter_coords(1) < 0.0)
            {
                random_angle = RandomNumberGenerator::Instance()->ranf();
                random_angle *= 2.0*M_PI;

                random_vector(0) = separation*cos(random_angle);
                random_vector(1) = separation*sin(random_angle);
                proposed_new_daughter_coords = parent_coords + random_vector;
            }
            daughter_coords = proposed_new_daughter_coords;
        }

        assert(daughter_coords(1) >= 0.0); // to make sure dividing cells stay in the cell population
        assert(parent_coords(1) >= 0.0);   // to make sure dividing cells stay in the cell population

        // Set the parent to use this location
        ChastePoint<2> parent_coords_point(parent_coords);

        unsigned node_index = mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else // using a VertexBasedCellPopulation
    {
        c_vector<double, 2> axis_of_division = zero_vector<double>(2);

        // We don't need to prescribe how 'stem' cells divide if Wnt is present
        bool is_wnt_included = WntConcentration<2>::Instance()->IsWntSetUp();
        if (!is_wnt_included)
        {
            WntConcentration<2>::Destroy();
            if (pParentCell->GetCellCycleModel()->GetCellProliferativeType() == STEM)
            {
                axis_of_division(0) = 1.0;
                axis_of_division(1) = 0.0;
            }
        }
        return axis_of_division;
    }
}

void CryptSimulation2d::SetupWriteBetaCatenin()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);
    mVizBetaCateninResultsFile = output_file_handler.OpenOutputFile("results.vizbetacatenin");
    *mpVizSetupFile << "BetaCatenin\n";
}

void CryptSimulation2d::WriteBetaCatenin(double time)
{
    *mVizBetaCateninResultsFile <<  time << "\t";

    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        unsigned global_index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        double x = mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[0];
        double y = mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

        // We should only be calling this code block if mWriteBetaCatenin has been set to true in the constructor
        assert(mWriteBetaCatenin);

        AbstractVanLeeuwen2009WntSwatCellCycleModel* p_model = dynamic_cast<AbstractVanLeeuwen2009WntSwatCellCycleModel*>(cell_iter->GetCellCycleModel());
        double b_cat_membrane = p_model->GetMembraneBoundBetaCateninLevel();
        double b_cat_cytoplasm = p_model->GetCytoplasmicBetaCateninLevel();
        double b_cat_nuclear = p_model->GetNuclearBetaCateninLevel();

        *mVizBetaCateninResultsFile << global_index << " " << x << " " << y << " " << b_cat_membrane << " " << b_cat_cytoplasm << " " << b_cat_nuclear << " ";
    }

    *mVizBetaCateninResultsFile << "\n";
}

void CryptSimulation2d::SetupSolve()
{
    // First call method on base class
    OffLatticeSimulation<2>::SetupSolve();

    /*
     * If there are any cells in the simulation, and mWriteBetaCatenin has been set
     * to true in the constructor, then set up the beta-catenin results file and
     * write the initial conditions to file.
     */
    bool any_cells_present = (mrCellPopulation.Begin() != mrCellPopulation.End());
    if (any_cells_present && mWriteBetaCatenin)
    {
        SetupWriteBetaCatenin();
        double current_time = SimulationTime::Instance()->GetTime();
        WriteBetaCatenin(current_time);
    }
}

void CryptSimulation2d::PostSolve()
{
    // First call method on base class
    OffLatticeSimulation<2>::PostSolve();

    SimulationTime* p_time = SimulationTime::Instance();

    if ((p_time->GetTimeStepsElapsed()+1)%mSamplingTimestepMultiple == 0)
    {
        /*
         * If there are any cells in the simulation, and mWriteBetaCatenin has been set
         * to true in the constructor, then set up the beta-catenin results file and
         * write the initial conditions to file.
         */
        bool any_cells_present = (mrCellPopulation.Begin() != mrCellPopulation.End());
        if (any_cells_present && mWriteBetaCatenin)
        {
            double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
            WriteBetaCatenin(time_next_step);
        }
    }
}

void CryptSimulation2d::AfterSolve()
{
    // First call method on base class
    OffLatticeSimulation<2>::AfterSolve();

    /*
     * If there are any cells in the simulation, and mWriteBetaCatenin has been set
     * to true in the constructor, then close the beta-catenin results file.
     */
    bool any_cells_present = (mrCellPopulation.Begin() != mrCellPopulation.End());
    if (any_cells_present && mWriteBetaCatenin)
    {
        mVizBetaCateninResultsFile->close();
    }
}

void CryptSimulation2d::UseJiggledBottomCells()
{
    // The CryptSimulationBoundaryCondition object is the first element of mBoundaryConditions
    boost::static_pointer_cast<CryptSimulationBoundaryCondition<2> >(mBoundaryConditions[0])->SetUseJiggledBottomCells(true);
}

void CryptSimulation2d::SetBottomCellAncestors()
{
    /*
     * We use a different height threshold depending on which type of cell
     * population we are using, a MeshBasedCellPopulationWithGhostNodes or
     * a VertexBasedCellPopulation.
     */
    double threshold_height = 1.0;
    if (mUsingMeshBasedCellPopulation)
    {
        threshold_height = 0.5;
    }

    unsigned index = 0;
    for (AbstractCellPopulation<2>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        if (mrCellPopulation.GetLocationOfCellCentre(*cell_iter)[1] < threshold_height)
        {
            cell_iter->SetAncestor(index++);
        }
    }
}

void CryptSimulation2d::OutputSimulationParameters(out_stream& rParamsFile)
{
    double width = mrCellPopulation.GetWidth(0);
    bool use_jiggled_bottom_cells = boost::static_pointer_cast<CryptSimulationBoundaryCondition<2> >(mBoundaryConditions[0])->GetUseJiggledBottomCells();

    *rParamsFile << "\t\t<CryptCircumference>" << width << "</CryptCircumference>\n";
    *rParamsFile << "\t\t<UseJiggledBottomCells>" << use_jiggled_bottom_cells << "</UseJiggledBottomCells>\n";

    // Call method on direct parent class
    OffLatticeSimulation<2>::OutputSimulationParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CryptSimulation2d)
