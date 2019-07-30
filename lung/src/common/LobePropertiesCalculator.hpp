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


#ifndef LOBE_PROPERTIES_CALCULATOR_HPP_
#define LOBE_PROPERTIES_CALCULATOR_HPP_

#include <map>
#include "LungTools.hpp"

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkSTLReader.h"

/**
 * Calculates a number of morphological properties for a set of lung lobes
 */
class LobePropertiesCalculator
{
public:
    /**
     * Constructor
     *
     */
    LobePropertiesCalculator();

    /**
     * Destructor
     */
    ~LobePropertiesCalculator();

    /**
     * Allows the user to add a lung lobe from a file
     *
     * This method assumes that the file is a .stl file.
     *
     * @param rFileName The name of the STL file representing the lobe surface
     * @param lungLocation Is this a right or a left lung lobe?
     * @param rName A string identifying the lobe
     */
    void AddLobe(const std::string& rFileName, LungLocation lungLocation, const std::string& rName);

    /**
     * Adds a lobe to the mesh generator
     *
     * @param pLobeSurface A VTK poly data containing the surface of the lobe
     * @param lungLocation Is this a right or a left lung lobe?
     * @param rName A string identifying the lobe
     */
    void AddLobe(vtkSmartPointer<vtkPolyData> pLobeSurface, LungLocation lungLocation, const std::string& rName);

    /**
     * @return The total volume of all the registered lobes
     */
    double GetTotalVolume();

    /**
     * @param rName The name of the lobe to get the volume for
     * @return The absolute volume of the given lobe
     */
    double GetLobeVolume(const std::string& rName);

    /**
     * @param pLobeSurface A VTK poly data containing the surface of the lobe
     * @return The absolute volume of the given lobe
     */
    double GetLobeVolume(vtkSmartPointer<vtkPolyData> pLobeSurface);


    /**
     * @param rName The name of the lobe to get the fraction of total volume for
     * @return The volume of the given lobe as a fraction of the total
     */
    double GetLobeVolumeFraction(const std::string& rName);

    /**
     * @param pLobeSurface A VTK poly data containing the surface of the lobe
     * @return The volume of the given lobe as a fraction of the total
     */
    double GetLobeVolumeFraction(vtkSmartPointer<vtkPolyData> pLobeSurface);

private:
    /** A map containing the lobar surface definitions */
    std::map< std::string, std::pair< vtkSmartPointer<vtkPolyData>, LungLocation> > mLobesMap;
};

#endif //CHASTE_VTK

#endif // LOBE_PROPERTIES_CALCULATOR
