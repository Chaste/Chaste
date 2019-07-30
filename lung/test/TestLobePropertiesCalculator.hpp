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

#ifndef _TESTLOBEPROPERTIESCALCULATOR_HPP_
#define _TESTLOBEPROPERTIESCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "LobePropertiesCalculator.hpp"

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkSTLReader.h"
#include "vtkTriangleFilter.h"


#endif //CHASTE_VTK

class TestLobePropertiesCalculator : public CxxTest::TestSuite
{
public:
    void TestVolumes()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

        LobePropertiesCalculator calculator;
        calculator.AddLobe(CreateSphere(1.0,0.0,-1.0,1.0), LEFT, "lll");
        calculator.AddLobe(CreateSphere(1.0,0.0,1.0,1.5), LEFT, "lul");
        calculator.AddLobe(CreateSphere(-1.0,0.0,-1.0,0.5), LEFT, "rll");
        calculator.AddLobe(CreateSphere(-1.0,0.0,1.0,2.0), LEFT, "rul");

        TS_ASSERT_DELTA(calculator.GetTotalVolume(), 4*M_PI*(1.0*1.0*1.0 + 1.5*1.5*1.5 + 0.5*0.5*0.5 + 2.0*2.0*2.0)/3, 0.5); //Resolution of the sphere isn't great, hence low tolerance

        TS_ASSERT_DELTA(calculator.GetLobeVolume("lll"), 4*M_PI/3, 2e-1);
        TS_ASSERT_DELTA(calculator.GetLobeVolume("lul"), 4*M_PI/3*1.5*1.5*1.5, 2e-1);
        TS_ASSERT_DELTA(calculator.GetLobeVolume("rll"), 4*M_PI/3*0.5*0.5*0.5, 2e-1);
        TS_ASSERT_DELTA(calculator.GetLobeVolume("rul"), 4*M_PI/3*2.0*2.0*2.0, 2e-1);

        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("lll"), 0.08, 1e-6);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("lul"), 0.27, 1e-6);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("rll"), 0.01, 1e-6);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("rul"), 0.64, 1e-6);

        //Just do some independent checks
        TS_ASSERT_DELTA(calculator.GetLobeVolume(CreateSphere(-1.0,0.0,1.0,1.0)), 4*M_PI/3, 2e-1);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction(CreateSphere(-1.0,0.0,1.0,1.0)), 0.08, 2e-1);
#endif
    }

    void TestRealVolumeFractions()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

        LobePropertiesCalculator calculator;
        calculator.AddLobe("lung/test/data/lll.stl", LEFT, "lll");
        calculator.AddLobe("lung/test/data/lul.stl", LEFT, "lul");
        calculator.AddLobe("lung/test/data/rll.stl", RIGHT, "rll");
        calculator.AddLobe("lung/test/data/rml.stl", RIGHT, "rml");
        calculator.AddLobe("lung/test/data/rul.stl", RIGHT, "rul");

        //Confirmed manually to be correct for this subject
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("lll"), 0.1686, 1e-3);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("lul"), 0.2820, 1e-3);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("rll"), 0.2264, 1e-3);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("rml"), 0.0987, 1e-3);
        TS_ASSERT_DELTA(calculator.GetLobeVolumeFraction("rul"), 0.2241, 1e-3);
#endif
    }

private:
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
    vtkSmartPointer<vtkPolyData> CreateSphere(double XCentre, double YCentre, double ZCentre, double radius)
    {
        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetCenter(XCentre, YCentre, ZCentre);
        sphere->SetRadius(radius);
        sphere->SetThetaResolution(50);
        sphere->SetPhiResolution(50);
        sphere->Update();

        vtkSmartPointer<vtkPolyData> sphere_data = sphere->GetOutput();
        return sphere_data;
    }
#endif
};

#endif /*_TESTLOBEPROPERTIESCALCULATOR_HPP_*/
