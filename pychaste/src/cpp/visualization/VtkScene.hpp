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

#ifndef VTKSCENE_HPP_
#define VTKSCENE_HPP_

#include <vector>

#include <vtkLookupTable.h>
#include <vtkOggTheoraWriter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkWindowToImageFilter.h>

#include "AbstractCellPopulation.hpp"
#include "CellPopulationPyChasteActorGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "SmartPointers.hpp"

class vtkPNGWriter;

/**
 * Renders a Chaste cell population with VTK for visualisation.
 *
 * A population is supplied via SetCellPopulation(), and its actors are built
 * by a CellPopulationPyChasteActorGenerator. The scene can be captured to an
 * in-memory PNG (GetSceneAsCharBuffer()) or, during a simulation, driven by a
 * VtkSceneModifier that calls RenderFrame() at each output step to write a
 * per-step PNG and/or an Ogg/Theora animation. The output mode is chosen with
 * SetSaveAsImages() / SetSaveAsAnimation() / SetIsInteractive().
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
    bool mIsInteractive = false;

    /**
     * Save as an animation
     */
    bool mSaveAsAnimation = false;

    /**
     * Save as an image
     */
    bool mSaveAsImages = false;

    /**
     * Has the renderer started
     */
    bool mHasStarted = false;

    /**
     * Has the animation writer been started but not yet ended
     */
    bool mAnimationWriterStarted = false;

    /**
     * The cell population
     */
    boost::shared_ptr<CellPopulationPyChasteActorGenerator<DIM> > mpCellPopulationGenerator;

public:

    /**
     * Constructor. Creates the renderer, render window and interactor, with a
     * white background and a trackball-camera interactor style.
     */
    VtkScene();

    /**
     * Destructor. Finalises any in-progress animation output.
     */
    virtual ~VtkScene();

    /**
     * Finalise and close the animation file, if one is being written.
     */
    void End();

    /**
     * Get the cell population actor generator
     * @return the cell population actor generator
     */
    boost::shared_ptr<CellPopulationPyChasteActorGenerator<DIM> > GetCellPopulationActorGenerator() const;

    /**
     * Return the render window (e.g. to set its size before capturing a frame)
     * @return the vtk render window
     */
    vtkSmartPointer<vtkRenderWindow> GetRenderWindow() const;

    /**
     * Return the renderer
     * @return the vtk renderer
     */
    vtkSmartPointer<vtkRenderer> GetRenderer() const;

    /**
     * Rebuild and render the scene off-screen and return it as an in-memory PNG,
     * for inline display (e.g. in a Jupyter notebook). Writes no output files.
     * @return the rendered scene as PNG bytes
     */
    vtkSmartPointer<vtkUnsignedCharArray> GetSceneAsCharBuffer();

    /**
     * Update the scene and write output: refreshes the actors via UpdateScene(),
     * then renders and writes an image/animation frame if enabled.
     * @param timeStep the current time step, for annotating output files
     */
    virtual void RenderFrame(unsigned timeStep = 0);

    /**
     * Set the cell population
     * @param pCellPopulation the cell population for rendering
     */
    void SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation);

    /**
     * Set whether to render to an interactive on-screen window (rather than
     * off-screen).
     * @param isInteractive whether to render interactively
     */
    void SetIsInteractive(bool isInteractive);

    /**
     * Set the output file path prefix for saved images/animations.
     * @param rPath the output path prefix
     */
    void SetOutputFilePath(const std::string& rPath);

    /**
     * Set whether to write the scene to an Ogg/Theora animation, one frame per
     * RenderFrame() call (the output path sets the .ogg file name).
     * @param saveAsAnimation whether to save an animation
     */
    void SetSaveAsAnimation(bool saveAsAnimation);

    /**
     * Set whether to write a PNG image per RenderFrame() call, named from the
     * output path and time step.
     * @param saveAsImages whether to save per-step images
     */
    void SetSaveAsImages(bool saveAsImages);

    /**
     * Perform the one-time scene setup: build the actors, select off-screen
     * rendering, and (if enabled) initialise the image/animation output pipeline.
     * Called automatically by RenderFrame() / GetSceneAsCharBuffer() on first use.
     */
    void Start();

    /**
     * Start the interactor's event loop to handle window interaction; blocks
     * until the window is closed. Only meaningful in interactive mode.
     */
    void StartInteractiveEventHandler();

    /**
     * Rebuild the scene's actors to reflect the current cell population, and
     * reset the camera. Does not render or write any output (see RenderFrame).
     */
    void UpdateScene();

private:

    /**
     * Create a PNG writer connected to the window-to-image filter, configured to
     * write the rendered frame to an in-memory buffer.
     * @return the configured writer
     */
    vtkSmartPointer<vtkPNGWriter> MakePngWriter();

    /**
     * Render the current frame off-screen and flag the window-to-image filter,
     * leaving it ready to be captured by a PNG/animation writer.
     */
    void RenderForCapture();
};

#endif // VTKSCENE_HPP_
