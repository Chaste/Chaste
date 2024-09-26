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

#ifndef VTKSCENE_HPP_
#define VTKSCENE_HPP_

#include <vector>

#include <vtkAutoInit.h>
#include <vtkLookupTable.h>
#include <vtkOggTheoraWriter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVersion.h>
#include <vtkWindowToImageFilter.h>

#if VTK_MAJOR_VERSION <= 6
VTK_MODULE_INIT(vtkRenderingOpenGL);
#else
VTK_MODULE_INIT(vtkRenderingOpenGL2);
#endif
VTK_MODULE_INIT(vtkRenderingFreeType);

#include "AbstractCellPopulation.hpp"
#include "CellPopulationPyChasteActorGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

/**
 * A simple VTK renderer for cell populations
 */
template <unsigned DIM>
class VtkScene
{
    /**
     * The vtk renderer
     */
    vtkSmartPointer<vtkRenderer> mpRenderer;

    /**
     * The vtk render window
     */
    vtkSmartPointer<vtkRenderWindow> mpRenderWindow;

    /**
     * The vtk render window interactor
     */
    vtkSmartPointer<vtkRenderWindowInteractor> mpRenderWindowInteractor;

    /**
     * The path for output
     */
    std::string mOutputFilePath;

    /**
     * The animation writer
     */
    vtkSmartPointer<vtkOggTheoraWriter> mAnimationWriter;

    /**
     * The image to window filter
     */
    vtkSmartPointer<vtkWindowToImageFilter> mWindowToImageFilter;

    /**
     * Is the rendering interactive
     */
    bool mIsInteractive;

    /**
     * Save as an animation
     */
    bool mSaveAsAnimation;

    /**
     * Save as an image
     */
    bool mSaveAsImages;

    /**
     * Has the renderer started
     */
    bool mHasStarted;

    /**
     * Add annotation
     */
    bool mAddAnnotations;

    /**
     * How often to update the renderer during a simulation
     */
    unsigned mOutputFrequency;

    /**
     * The cell population
     */
    boost::shared_ptr<CellPopulationPyChasteActorGenerator<DIM> > mpCellPopulationGenerator;

public:
    /**
     * Constructor
     */
    VtkScene();

    /**
     * Destructor
     */
    virtual ~VtkScene();

    /**
     * Shut down the scene and close the animation
     */
    void End();

    /**
     * Render the current scene and return it as a char array that can be passed
     * into a Python buffer for display.
     * @return the scene as a char array
     */
    vtkSmartPointer<vtkUnsignedCharArray> GetSceneAsCharBuffer();

    /**
     * Return the renderer
     * @return the vtk renderer
     */
    vtkSmartPointer<vtkRenderer> GetRenderer();

    /**
     * Get the cell population actor generator
     * @return the cell population actor generator
     */
    boost::shared_ptr<CellPopulationPyChasteActorGenerator<DIM> > GetCellPopulationActorGenerator();

    /**
     * Update the renderer, this will update the population actor and write output images
     * @param timeStep the curren time step, for annotating output files
     */
    virtual void ResetRenderer(unsigned timeStep = 0);

    /**
     * Render the scene
     */
    void Start();

    /**
     * Set the cell population
     * @param pCellPopulation the cell population for rendering
     */
    void SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation);

    /**
     * Set the path for output
     * @param rPath the path for output
     */
    void SetOutputFilePath(const std::string& rPath);

    /**
     * Set run as an interactive window
     * @param isInteractive run as an interactive window
     */
    void SetIsInteractive(bool isInteractive);

    /**
     * Whether to save as an animation
     * @param saveAsAnimation save as an animation
     */
    void SetSaveAsAnimation(bool saveAsAnimation);

    /**
     * Whether to save as images (default)
     * @param saveAsImages save as images
     */
    void SetSaveAsImages(bool saveAsImages);

    /**
     * Start the event handler for window interaction
     */
    void StartInteractiveEventHandler();
};

#endif // VTKSCENE_HPP_
