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

#ifndef FORWARDEULERDEVICESOLVER_HPP_
#define FORWARDEULERDEVICESOLVER_HPP_

// This header contains a kernel, so it can only be included in files compiled by a GPU compiler (.cu files)

#include "CardiacCellDeviceData.hpp"
#include "DeviceTools.hpp"

/**
 * Kernel that advances each cell with the forward Euler method, using one thread per cell.
 *
 * DeviceModel is the generated struct of device functions for a cell model, e.g. CellLuoRudy1991FromCellMLDevice.
 *
 * @param data  the cells to advance and the time interval to solve over
 */
template <typename DeviceModel, typename ValueType>
__global__ void ForwardEulerKernel(CardiacCellDeviceData<ValueType> data)
{
    const unsigned cell = blockIdx.x * blockDim.x + threadIdx.x;
    if (cell >= data.NumCells)
    {
        return;
    }

    const unsigned num_states = DeviceModel::NUM_STATE_VARIABLES;
    ValueType* p_cell_states = data.pStates + cell * num_states;
    const ValueType* p_cell_parameters = data.pParameters + cell * DeviceModel::NUM_PARAMETERS;
    const ValueType stimulus = data.pStimuli[cell];
    const ValueType fixed_voltage = data.pFixedVoltages[cell];

    ValueType y[DeviceModel::NUM_STATE_VARIABLES];
    ValueType dy[DeviceModel::NUM_STATE_VARIABLES];
    for (unsigned i = 0; i < num_states; i++)
    {
        y[i] = p_cell_states[i];
    }

    for (unsigned step = 0; step < data.NumSteps; step++)
    {
        const ValueType time = data.StartTime + ValueType(step) * data.TimeStep;
        const ValueType dt = (step + 1u == data.NumSteps) ? data.LastTimeStep : data.TimeStep;
        DeviceModel::EvaluateYDerivatives(time, y, dy, p_cell_parameters, stimulus, data.Capacitance,
                                          data.SetVoltageDerivativeToZero, fixed_voltage);
        for (unsigned i = 0; i < num_states; i++)
        {
            y[i] += dt * dy[i];
        }
    }

    for (unsigned i = 0; i < num_states; i++)
    {
        p_cell_states[i] = y[i];
    }
}

/**
 * Advance a set of cells on the GPU with the forward Euler method, returning when the solve has finished.
 *
 * DeviceModel is the generated struct of device functions for a cell model, e.g. CellLuoRudy1991FromCellMLDevice.
 *
 * @param rData  the cells to advance and the time interval to solve over
 */
template <typename DeviceModel, typename ValueType>
void SolveCellsForwardEulerOnDevice(const CardiacCellDeviceData<ValueType>& rData)
{
    const unsigned threads_per_block = 128u;
    const unsigned num_blocks = (rData.NumCells + threads_per_block - 1u) / threads_per_block;
    ForwardEulerKernel<DeviceModel, ValueType><<<num_blocks, threads_per_block>>>(rData);
    DeviceTools::WaitForKernel();
}

#endif // FORWARDEULERDEVICESOLVER_HPP_
