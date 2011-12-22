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

#include "AbstractConvergenceTester.hpp"
#include "Exception.hpp"


AbstractUntemplatedConvergenceTester::AbstractUntemplatedConvergenceTester()
    : mMeshWidth(0.2),//cm
      OdeTimeStep(0.0025),//Justification from 1D test with this->PdeTimeStep held at 0.01 (allowing two hits at convergence)
      PdeTimeStep(0.005),//Justification from 1D test with this->OdeTimeStep held at 0.0025
      MeshNum(5u),//Justification from 1D test
      RelativeConvergenceCriterion(1e-4),
      LastDifference(1.0),
      Apd90FirstQn(0.0),
      Apd90ThirdQn(0.0),
      ConductionVelocity(0.0),
      PopulatedResult(false),
      FixedResult(false),
      UseAbsoluteStimulus(false),
      AbsoluteStimulus(-1e7),
      SimulateFullActionPotential(false),
      Converged(false),
      Stimulus(PLANE),
      NeumannStimulus(4000)
{
}





AbstractUntemplatedConvergenceTester::~AbstractUntemplatedConvergenceTester()
{
}

