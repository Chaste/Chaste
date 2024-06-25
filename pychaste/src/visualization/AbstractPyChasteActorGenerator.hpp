/*

Copyright (c) 2005-2024, University of Oxford.
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
#include "SmartPointers.hpp"
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkColorTransferFunction.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include "UblasVectorInclude.hpp"

/**
 * This class generates VTK actors for geometric features.
 */
template<unsigned DIM>
class AbstractPyChasteActorGenerator
{

protected:

    /**
     * The color lookup for continuous entities
     */
    vtkSmartPointer<vtkColorTransferFunction> mpColorTransferFunction;

    /**
     * The color lookup for discrete entities
     */
    vtkSmartPointer<vtkColorTransferFunction> mpDiscreteColorTransferFunction;

    /**
     * Show the edges, using a tube filter
     */
    bool mShowEdges;

    /**
     * Show the points, using a glyph filter
     */
    bool mShowPoints;

    /**
     * Show the volume
     */
    bool mShowVolume;

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
    double mVolumeOpacity;

    /**
     * The default size for points
     */
    double mPointSize;

    /**
     * The default size for edges
     */
    double mEdgeSize;

    /**
     * The label for contouring on data
     */
    std::string mDataLabel;

    /**
     * The scale bar
     */
    vtkSmartPointer<vtkScalarBarActor> mpScaleBar;

    /**
     * Show the scale bar
     */
    bool mShowScaleBar;

public:

    /**
     * Constructor
     */
    AbstractPyChasteActorGenerator();

    /**
     * Destructor
     */
    virtual ~AbstractPyChasteActorGenerator();

    /**
     * @return return the color transfer function
     */
    vtkSmartPointer<vtkColorTransferFunction> GetColorTransferFunction();

    /**
     * @return return the discrete color transfer function
     */
    vtkSmartPointer<vtkColorTransferFunction> GetDiscreteColorTransferFunction();

    /**
     * @return return the scale bar
     */
    vtkSmartPointer<vtkScalarBarActor> GetScaleBar();

    /**
     * Add the Abstract actor to the renderer
     * @param pRenderer the current renderer
     */
    virtual void AddActor(vtkSmartPointer<vtkRenderer> pRenderer) = 0;

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
     * Set whether to show the volume
     * @param show whether to show the volumes
     */
    void SetShowVolume(bool show);

    /**
     * Set the edge color in RGB (e.g. (255,255,255) is white)
     * @param rColor the edge color
     */
    void SetEdgeColor(const c_vector<double, 3>& rColor);

    /**
     * Set the point color in RGB (e.g. (255,255,255) is white)
     * @param rColor the point color
     */
    void SetPointColor(const c_vector<double, 3>& rColor);

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

    /**
     * Set the default point size
     * @param size the default point size
     */
    void SetPointSize(double size);

    /**
     * Set the default edge size
     * @param size the default edge size
     */
    void SetEdgeSize(double size);

    /**
     * Set the label for contouring data
     * @param rLabel the label for contouring data
     */
    void SetDataLabel(const std::string& rLabel);

    /**
     * Set show scale bar
     * @param show show scale bar
     */
    void SetShowScaleBar(double show);

};

#endif /* ABSTRACTPYCHASTEACTORGENERATOR_HPP_*/
