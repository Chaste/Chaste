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


#include "LobePropertiesCalculator.hpp"

#ifdef CHASTE_VTK

#include "vtkMassProperties.h"
#include "vtkSTLReader.h"

LobePropertiesCalculator::LobePropertiesCalculator()
{
}

LobePropertiesCalculator::~LobePropertiesCalculator()
{
}

void LobePropertiesCalculator::AddLobe(const std::string& rFileName, LungLocation lungLocation, const std::string& rName)
{
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(rFileName.c_str());
    reader->Update();

    AddLobe(reader->GetOutput(), lungLocation, rName);
}


void LobePropertiesCalculator::AddLobe(vtkSmartPointer<vtkPolyData> pLobeSurface, LungLocation lungLocation, const std::string& rName)
{
    mLobesMap[rName] = std::make_pair(pLobeSurface, lungLocation);
}

double LobePropertiesCalculator::GetTotalVolume()
{
    double total_volume = 0.0;

    std::map< std::string, std::pair< vtkSmartPointer<vtkPolyData>, LungLocation> >::iterator iter;
    for (iter = mLobesMap.begin(); iter != mLobesMap.end(); ++iter)
    {
        total_volume += GetLobeVolume(iter->second.first);
    }

    return total_volume;
}

double LobePropertiesCalculator::GetLobeVolume(const std::string& rName)
{
    return GetLobeVolume(mLobesMap[rName].first);
}

double LobePropertiesCalculator::GetLobeVolume(vtkSmartPointer<vtkPolyData> pLobeSurface)
{
    vtkSmartPointer<vtkMassProperties> mass_properties = vtkSmartPointer<vtkMassProperties>::New();
#if VTK_MAJOR_VERSION >= 6
        mass_properties->SetInputData(pLobeSurface);
#else
        mass_properties->SetInput(pLobeSurface);
#endif

    return mass_properties->GetVolume();
}

double LobePropertiesCalculator::GetLobeVolumeFraction(const std::string& rName)
{
    return GetLobeVolume(rName)/GetTotalVolume();
}

double LobePropertiesCalculator::GetLobeVolumeFraction(vtkSmartPointer<vtkPolyData> pLobeSurface)
{
    return GetLobeVolume(pLobeSurface)/GetTotalVolume();
}

#endif // CHASTE_VTK


