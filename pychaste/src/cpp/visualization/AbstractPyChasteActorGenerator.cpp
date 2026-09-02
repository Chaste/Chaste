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

#include "AbstractPyChasteActorGenerator.hpp"

#include <vtkColorTransferFunction.h>
#include <vtkCoordinate.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>

#include "PyChasteColormaps.hpp"
#include "UblasIncludes.hpp"
#include "UblasVectorInclude.hpp"

template <unsigned DIM>
AbstractPyChasteActorGenerator<DIM>::AbstractPyChasteActorGenerator()
        : mpContinuousColorTransferFunction(vtkSmartPointer<vtkColorTransferFunction>::New()),
          mpDiscreteColorTransferFunction(vtkSmartPointer<vtkColorTransferFunction>::New()),
          mEdgeColor(zero_vector<double>(3)),
          mPointColor(unit_vector<double>(3, 2)),
          mVolumeColor(unit_vector<double>(3, 0)),
          mpScaleBar(vtkSmartPointer<vtkScalarBarActor>::New())
{
    // The default point and volume colours are initialised as unit RGB vectors
    // (blue and red respectively); scale them to the 0-255 range used elsewhere.
    mPointColor *= 255.0;
    mVolumeColor *= 255.0;

    // Build the continuous colour transfer function from the viridis colourmap,
    // sampling the unit interval [0, 1] across the 256 entries.
    for (unsigned idx = 0; idx < 256; idx++)
    {
        mpContinuousColorTransferFunction->AddRGBPoint(double(idx) / 255.0,
                                                       viridis_colors[idx][0],
                                                       viridis_colors[idx][1],
                                                       viridis_colors[idx][2]);
    }

    // Build the discrete colour transfer function from the accent colourmap.
    for (unsigned idx = 0; idx < 256; idx++)
    {
        mpDiscreteColorTransferFunction->AddRGBPoint(double(idx) / 255.0,
                                                     accent_colors[idx][0],
                                                     accent_colors[idx][1],
                                                     accent_colors[idx][2]);
    }

    // Configure the scale bar: a horizontal bar centred near the top of the
    // viewport, with plain (non-italic, non-bold) black title and labels.
    mpScaleBar->SetOrientationToHorizontal();
    mpScaleBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    mpScaleBar->GetPositionCoordinate()->SetValue(0.25, 0.84);
    mpScaleBar->SetWidth(0.5);
    mpScaleBar->SetHeight(0.1);
    mpScaleBar->GetTitleTextProperty()->ItalicOff();
    mpScaleBar->GetLabelTextProperty()->ItalicOff();
    mpScaleBar->GetTitleTextProperty()->BoldOff();
    mpScaleBar->GetLabelTextProperty()->BoldOff();
    mpScaleBar->SetLabelFormat("%.2g");
    mpScaleBar->GetTitleTextProperty()->SetFontSize(5.0);
    mpScaleBar->GetLabelTextProperty()->SetFontSize(5.0);
    mpScaleBar->GetTitleTextProperty()->SetColor(0.0, 0.0, 0.0);
    mpScaleBar->GetLabelTextProperty()->SetColor(0.0, 0.0, 0.0);
}

template <unsigned DIM>
vtkSmartPointer<vtkColorTransferFunction> AbstractPyChasteActorGenerator<DIM>::GetContinuousColorTransferFunction() const
{
    return mpContinuousColorTransferFunction;
}

template <unsigned DIM>
vtkSmartPointer<vtkColorTransferFunction> AbstractPyChasteActorGenerator<DIM>::GetDiscreteColorTransferFunction() const
{
    return mpDiscreteColorTransferFunction;
}

template <unsigned DIM>
vtkSmartPointer<vtkScalarBarActor> AbstractPyChasteActorGenerator<DIM>::GetScaleBar() const
{
    return mpScaleBar;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetDataLabel(const std::string& rLabel)
{
    mDataLabel = rLabel;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetEdgeColor(const c_vector<double, 3>& rColor)
{
    mEdgeColor = rColor;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetEdgeSize(double size)
{
    mEdgeSize = size;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetPointColor(const c_vector<double, 3>& rColor)
{
    mPointColor = rColor;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetPointSize(double size)
{
    mPointSize = size;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetShowEdges(bool show)
{
    mShowEdges = show;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetShowPoints(bool show)
{
    mShowPoints = show;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetShowScaleBar(double show)
{
    mShowScaleBar = show;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetShowVolume(bool show)
{
    mShowVolume = show;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetVolumeColor(const c_vector<double, 3>& rColor)
{
    mVolumeColor = rColor;
}

template <unsigned DIM>
void AbstractPyChasteActorGenerator<DIM>::SetVolumeOpacity(double opacity)
{
    mVolumeOpacity = opacity;
}

template class AbstractPyChasteActorGenerator<2>;
template class AbstractPyChasteActorGenerator<3>;
