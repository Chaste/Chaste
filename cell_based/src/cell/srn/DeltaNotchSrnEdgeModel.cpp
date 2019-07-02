#include "DeltaNotchSrnEdgeModel.hpp"

DeltaNotchSrnEdgeModel::DeltaNotchSrnEdgeModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractOdeSrnModel(2, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchSrnEdgeModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchSrnEdgeModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

DeltaNotchSrnEdgeModel::DeltaNotchSrnEdgeModel(const DeltaNotchSrnEdgeModel& rModel)
        : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */

    assert(rModel.GetOdeSystem());
    SetOdeSystem(new DeltaNotchOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* DeltaNotchSrnEdgeModel::CreateSrnModel()
{
    return new DeltaNotchSrnEdgeModel(*this);
}

void DeltaNotchSrnEdgeModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateDeltaNotch();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void DeltaNotchSrnEdgeModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new DeltaNotchOdeSystem);
}

void DeltaNotchSrnEdgeModel::UpdateDeltaNotch()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double mean_delta = mpCell->GetCellEdgeData()->GetItem("mean delta")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("Mean Delta", mean_delta);
}

double DeltaNotchSrnEdgeModel::GetNotch()
{
    assert(mpOdeSystem != nullptr);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

void DeltaNotchSrnEdgeModel::SetNotch(double value) {
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[0] = value;
}

double DeltaNotchSrnEdgeModel::GetDelta()
{
    assert(mpOdeSystem != nullptr);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

void DeltaNotchSrnEdgeModel::SetDelta(double value) {
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[1] = value;

}

double DeltaNotchSrnEdgeModel::GetMeanNeighbouringDelta()
{
    assert(mpOdeSystem != nullptr);
    double mean_neighbouring_delta = mpOdeSystem->GetParameter("Mean Delta");
    return mean_neighbouring_delta;
}

void DeltaNotchSrnEdgeModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}





// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchSrnEdgeModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchSrnEdgeModel)