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

#include "Nash2004ContractionModel.hpp"

const double Nash2004ContractionModel::e0ByT0 = 1.0/100; //See class doxygen for discussion on where this comes from
const double Nash2004ContractionModel::kTa = 47.9;

template<>
void OdeSystemInformation<Nash2004ContractionModel>::Initialise()
{
    // this is definitely (currently) covered [put an assert(0) in here and run
    // TestContractionModels and it will fail] [also see the call of
    //   this->mpSystemInfo = OdeSystemInformation<Nash2004ContractionModel>::Instance()
    // in this class' constructor], but for some reason gcov doesn't see it.
    #define COVERAGE_IGNORE
    this->mVariableNames.push_back("Ta");
    this->mVariableUnits.push_back("kPa");
    this->mInitialised = true;
    #undef COVERAGE_IGNORE
}
