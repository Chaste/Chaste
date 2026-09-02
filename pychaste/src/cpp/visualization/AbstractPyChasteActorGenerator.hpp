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

#ifndef ABSTRACTPYCHASTEACTORGENERATOR_HPP_
#define ABSTRACTPYCHASTEACTORGENERATOR_HPP_

#include <vector>

#include <vtkColorTransferFunction.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkTextProperty.h>

#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"

/**
 * Abstract base for the PyChaste actor generators.
 *
 * It does not build any actors itself - AddActor() is pure virtual and
 * implemented by each concrete generator. Instead it holds the rendering
 * configuration they share: the continuous and discrete colour transfer
 * functions, the edge/point/volume colours, sizes and opacity, the data label
 * to colour by, and the scale bar. The constructor seeds these with sensible
 * defaults, and the Set... methods adjust them.
 */
template <unsigned DIM>
class AbstractPyChasteActorGenerator
{

protected:
    /**
     * The color lookup for continuous entities
     */
    vtkSmartPointer<vtkColorTransferFunction> mpContinuousColorTransferFunction;

    /**
     * The color lookup for discrete entities
     */
    vtkSmartPointer<vtkColorTransferFunction> mpDiscreteColorTransferFunction;

    /**
     * Show the edges, using a tube filter
     */
    bool mShowEdges = true;

    /**
     * Show the points, using a glyph filter
     */
    bool mShowPoints = false;

    /**
     * Show the volume
     */
    bool mShowVolume = true;

    /**
     * The edge color in RGB
     */
    c_vector<double, 3> mEdgeColor;

    /**
     * The point color in RGB
     */
    c_vector<double, 3> mPointColor;

    /**
     * The volume color in RGB
     */
    c_vector<double, 3> mVolumeColor;

    /**
     * The volume opacity
     */
    double mVolumeOpacity = 0.9;

    /**
     * The default size for points
     */
    double mPointSize = 0.5;

    /**
     * The default size for edges
     */
    double mEdgeSize = 0.01;

    /**
     * The name of the data field to colour by, when colouring by data
     */
    std::string mDataLabel;

    /**
     * The scale bar
     */
    vtkSmartPointer<vtkScalarBarActor> mpScaleBar;

    /**
     * Show the scale bar
     */
    bool mShowScaleBar = false;

public:

    /**
     * Constructor. Initialises the continuous (viridis) and discrete (accent)
     * colour transfer functions and the scale bar with their default settings.
     */
    AbstractPyChasteActorGenerator();

    /**
     * Destructor
     */
    virtual ~AbstractPyChasteActorGenerator() = default;

    /**
     * Build the generated actors and add them to the renderer. Pure virtual:
     * each concrete generator adds the actors appropriate to what it renders.
     * @param pRenderer the renderer to add the actors to
     */
    virtual void AddActor(vtkSmartPointer<vtkRenderer> pRenderer) = 0;

    /**
     * @return the continuous (data) colour transfer function
     */
    vtkSmartPointer<vtkColorTransferFunction> GetContinuousColorTransferFunction() const;

    /**
     * @return the discrete (cell id/type) colour transfer function
     */
    vtkSmartPointer<vtkColorTransferFunction> GetDiscreteColorTransferFunction() const;

    /**
     * @return the scale bar
     */
    vtkSmartPointer<vtkScalarBarActor> GetScaleBar() const;

    /**
     * Set the name of the data field to colour by (used when colouring by data)
     * @param rLabel the data field name
     */
    void SetDataLabel(const std::string& rLabel);

    /**
     * Set the edge color in RGB (e.g. (255,255,255) is white)
     * @param rColor the edge color
     */
    void SetEdgeColor(const c_vector<double, 3>& rColor);

    /**
     * Set the default edge size
     * @param size the default edge size
     */
    void SetEdgeSize(double size);

    /**
     * Set the point color in RGB (e.g. (255,255,255) is white)
     * @param rColor the point color
     */
    void SetPointColor(const c_vector<double, 3>& rColor);

    /**
     * Set the default point size
     * @param size the default point size
     */
    void SetPointSize(double size);

    /**
     * Set whether to show the edges
     * @param show whether to show the edges
     */
    void SetShowEdges(bool show);

    /**
     * Set whether to show the points
     * @param show whether to show the points
     */
    void SetShowPoints(bool show);

    /**
     * Set whether to show the scale bar
     * @param show whether to show the scale bar
     */
    void SetShowScaleBar(double show);

    /**
     * Set whether to show the volume
     * @param show whether to show the volume
     */
    void SetShowVolume(bool show);

    /**
     * Set the volume color in RGB (e.g. (255,255,255) is white)
     * @param rColor the volume color
     */
    void SetVolumeColor(const c_vector<double, 3>& rColor);

    /**
     * Set the opacity for the volume
     * @param opacity the opacity for the volume
     */
    void SetVolumeOpacity(double opacity);
};

#endif // ABSTRACTPYCHASTEACTORGENERATOR_HPP_
