#include "DeltaNotchOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"
#include "Debug.hpp"

DeltaNotchOdeSystem::DeltaNotchOdeSystem(double meanDelta, std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<DeltaNotchOdeSystem>);

    /**
     * The state variables are as follows:
     * 
     * 0 - Notch concentration for this cell
     * 1 - Delta concentration for this cell
     * 2 - average Delta concentration for this cell's immediate neighbours
     * 
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1.0); // soon overwritten
    SetDefaultInitialCondition(1, 1.0); // soon overwritten
    SetDefaultInitialCondition(2, 0.5);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

DeltaNotchOdeSystem::~DeltaNotchOdeSystem()
{
}

void DeltaNotchOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double notch = rY[0];
    double delta = rY[1];
    double mean_delta = rY[2];


    // The next two lines define the ODE system by Collier et al. (1996)
    double dx1 = mean_delta*mean_delta/(0.01 + mean_delta*mean_delta) - notch;
    double dx2 = 1/(1 + 100*notch*notch) - delta;

    rDY[0] = dx1;
    rDY[1] = dx2;
    rDY[2] = 0.0; // don't change the mean Delta level over the course of the mechanics time step
}

template<>
void CellwiseOdeSystemInformation<DeltaNotchOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Mean Delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchOdeSystem)
