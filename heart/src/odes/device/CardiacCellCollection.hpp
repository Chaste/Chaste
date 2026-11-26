/*

Copyright (c) 2005-2026, University of Oxford.
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


#ifndef CARDIACCELLCOLLECTION_HPP_
#define CARDIACCELLCOLLECTION_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>

#include "CardiacCellDeviceData.hpp"
#include "DeviceTools.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "TimeStepper.hpp"

/**
 * A set of cardiac cells of one model which are solved together on a GPU.
 *
 * CellType is a cell class generated from CellML with GPU support, and it determines both the equations and
 * the numerical method: for example, CellLuoRudy1991FromCellML gives the Luo-Rudy 1991 model solved with
 * forward Euler. Any solver given to the cells themselves is not used. ValueType is the floating point type
 * used on the GPU, float or double.
 *
 * Each solve copies the cells' state variables, parameters and stimuli to the GPU, advances them there, and
 * copies the new state variables back, so the cells can be read and changed between solves as usual. The
 * stimulus is evaluated at the start of each solve and held constant over it, so a solve should be no longer
 * than the time over which the stimulus can be taken as constant, e.g. one PDE time step.
 */
template <typename CellType, typename ValueType>
class CardiacCellCollection
{
private:
    /** The cells. */
    std::vector<boost::shared_ptr<CellType> > mCells;

    /** The number of state variables of a cell. */
    unsigned mNumStateVariables;

    /** The number of parameters of a cell. */
    unsigned mNumParameters;

    /** The time step to solve with (ms). Defaults to the ODE time step in HeartConfig. */
    double mDt;

    /** State variables of all the cells, cell after cell, in the type used on the GPU. */
    std::vector<ValueType> mStates;

    /** Parameters of all the cells, cell after cell. */
    std::vector<ValueType> mParameters;

    /** The stimulus of each cell. */
    std::vector<ValueType> mStimuli;

    /** The voltage of each cell. */
    std::vector<ValueType> mFixedVoltages;

    /** Device memory for #mStates. */
    ValueType* mpDeviceStates;

    /** Device memory for #mParameters. */
    ValueType* mpDeviceParameters;

    /** Device memory for #mStimuli. */
    ValueType* mpDeviceStimuli;

    /** Device memory for #mFixedVoltages. */
    ValueType* mpDeviceFixedVoltages;

    /**
     * Allocate device memory for a host array.
     *
     * @param rHost  the host array
     * @return a pointer to device memory of the same size
     */
    ValueType* AllocateOnDevice(const std::vector<ValueType>& rHost)
    {
        return static_cast<ValueType*>(DeviceTools::Allocate(rHost.size() * sizeof(ValueType)));
    }

    /**
     * Copy a host array to the device.
     *
     * @param pDevice  the device memory for the array
     * @param rHost  the host array
     */
    void CopyToDevice(ValueType* pDevice, const std::vector<ValueType>& rHost)
    {
        DeviceTools::CopyToDevice(pDevice, rHost.data(), rHost.size() * sizeof(ValueType));
    }

    /**
     * Advance all the cells on the GPU.
     *
     * @param tStart  beginning of the time interval to solve over (ms)
     * @param tEnd  end of the time interval to solve over (ms)
     * @param setVoltageDerivativeToZero  whether to keep the voltage of each cell fixed
     */
    void Solve(double tStart, double tEnd, bool setVoltageDerivativeToZero);

public:
    /**
     * Create a collection and allocate GPU memory for it.
     *
     * @param rCells  the cells to solve together, at least one
     */
    CardiacCellCollection(const std::vector<boost::shared_ptr<CellType> >& rCells);

    /**
     * Release the GPU memory. The cells are not affected.
     */
    ~CardiacCellCollection();

    /** Collections own GPU memory, so they cannot be copied. */
    CardiacCellCollection(const CardiacCellCollection&) = delete;

    /** Collections own GPU memory, so they cannot be assigned. */
    CardiacCellCollection& operator=(const CardiacCellCollection&) = delete;

    /**
     * @return the number of cells in the collection
     */
    unsigned GetNumberOfCells() const
    {
        return mCells.size();
    }

    /**
     * Set the time step to solve with.
     *
     * @param dt  the time step (ms)
     */
    void SetTimestep(double dt)
    {
        mDt = dt;
    }

    /**
     * Simulate all the cells between the time interval [tStart, tEnd], updating their state variables.
     *
     * @param tStart  beginning of the time interval to simulate (ms)
     * @param tEnd  end of the time interval to simulate (ms)
     */
    void SolveAndUpdateState(double tStart, double tEnd)
    {
        Solve(tStart, tEnd, false);
    }

    /**
     * Simulate all the cells between the time interval [tStart, tEnd], but without updating their voltages.
     *
     * @param tStart  beginning of the time interval to simulate (ms)
     * @param tEnd  end of the time interval to simulate (ms)
     */
    void ComputeExceptVoltage(double tStart, double tEnd)
    {
        Solve(tStart, tEnd, true);
    }
};

template <typename CellType, typename ValueType>
CardiacCellCollection<CellType, ValueType>::CardiacCellCollection(const std::vector<boost::shared_ptr<CellType> >& rCells)
    : mCells(rCells),
      mDt(HeartConfig::Instance()->GetOdeTimeStep()),
      mpDeviceStates(NULL),
      mpDeviceParameters(NULL),
      mpDeviceStimuli(NULL),
      mpDeviceFixedVoltages(NULL)
{
    if (mCells.empty())
    {
        EXCEPTION("A CardiacCellCollection needs at least one cell");
    }
    mNumStateVariables = mCells[0]->GetNumberOfStateVariables();
    mNumParameters = mCells[0]->GetNumberOfParameters();

    mStates.resize(mCells.size() * mNumStateVariables);
    mParameters.resize(mCells.size() * mNumParameters);
    mStimuli.resize(mCells.size());
    mFixedVoltages.resize(mCells.size());

    mpDeviceStates = AllocateOnDevice(mStates);
    mpDeviceParameters = AllocateOnDevice(mParameters);
    mpDeviceStimuli = AllocateOnDevice(mStimuli);
    mpDeviceFixedVoltages = AllocateOnDevice(mFixedVoltages);
}

template <typename CellType, typename ValueType>
CardiacCellCollection<CellType, ValueType>::~CardiacCellCollection()
{
    DeviceTools::Free(mpDeviceStates);
    DeviceTools::Free(mpDeviceParameters);
    DeviceTools::Free(mpDeviceStimuli);
    DeviceTools::Free(mpDeviceFixedVoltages);
}

template <typename CellType, typename ValueType>
void CardiacCellCollection<CellType, ValueType>::Solve(double tStart, double tEnd, bool setVoltageDerivativeToZero)
{
    // Take the same time steps as the CPU solvers: steps of mDt, with a shorter last step if needed
    TimeStepper stepper(tStart, tEnd, mDt);
    unsigned num_steps = 0;
    double last_time_step = 0.0;
    while (!stepper.IsTimeAtEnd())
    {
        last_time_step = stepper.GetNextTimeStep();
        stepper.AdvanceOneTimeStep();
        num_steps++;
    }

    const unsigned num_cells = mCells.size();
    for (unsigned cell = 0; cell < num_cells; cell++)
    {
        const std::vector<double>& r_states = mCells[cell]->rGetStateVariables();
        for (unsigned i = 0; i < mNumStateVariables; i++)
        {
            mStates[cell * mNumStateVariables + i] = r_states[i];
        }
        for (unsigned i = 0; i < mNumParameters; i++)
        {
            mParameters[cell * mNumParameters + i] = mCells[cell]->GetParameter(i);
        }
        mStimuli[cell] = mCells[cell]->GetIntracellularAreaStimulus(tStart);
        mFixedVoltages[cell] = mCells[cell]->GetVoltage();
    }
    CopyToDevice(mpDeviceStates, mStates);
    CopyToDevice(mpDeviceParameters, mParameters);
    CopyToDevice(mpDeviceStimuli, mStimuli);
    CopyToDevice(mpDeviceFixedVoltages, mFixedVoltages);

    CardiacCellDeviceData<ValueType> data;
    data.NumCells = num_cells;
    data.pStates = mpDeviceStates;
    data.pParameters = mpDeviceParameters;
    data.pStimuli = mpDeviceStimuli;
    data.pFixedVoltages = mpDeviceFixedVoltages;
    data.Capacitance = HeartConfig::Instance()->GetCapacitance();
    data.SetVoltageDerivativeToZero = setVoltageDerivativeToZero;
    data.StartTime = tStart;
    data.TimeStep = mDt;
    data.NumSteps = num_steps;
    data.LastTimeStep = last_time_step;
    CellType::SolveOnDevice(data);

    DeviceTools::CopyToHost(mStates.data(), mpDeviceStates, mStates.size() * sizeof(ValueType));
    for (unsigned cell = 0; cell < num_cells; cell++)
    {
        const double voltage = mCells[cell]->GetVoltage();
        std::vector<double>& r_states = mCells[cell]->rGetStateVariables();
        for (unsigned i = 0; i < mNumStateVariables; i++)
        {
            r_states[i] = mStates[cell * mNumStateVariables + i];
        }
        if (setVoltageDerivativeToZero)
        {
            // The voltage has not changed, but may have lost precision on its way through the GPU
            mCells[cell]->SetVoltage(voltage);
#ifndef NDEBUG
            // As in AbstractCardiacCell::ComputeExceptVoltage
            mCells[cell]->VerifyStateVariables();
#endif // NDEBUG
        }
    }
}

#endif // CARDIACCELLCOLLECTION_HPP_
