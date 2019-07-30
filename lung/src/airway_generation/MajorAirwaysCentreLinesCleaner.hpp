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


#ifndef MAJOR_AIRWAYS_CENTRE_LINES_CLEANER_HPP_
#define MAJOR_AIRWAYS_CENTRE_LINES_CLEANER_HPP_

#include "MutableMesh.hpp"
#include "AirwayTreeWalker.hpp"
#include "AirwayPropertiesCalculator.hpp"

/**
 * Class to trim back centrelines of major airways.
 *
 * Segmentations of the airways decrease in quality towards the most distal segmented branches.
 * Furthermore, the airway generator always starts generating with a bifurcation. This may not
 * be appropriate for a partially segmented airway. The MajorAirwaysCentreLinesCleaner allows
 * the removal of airways in a major airways centerline below a given order. The model will
 * be trimmed back to the bifurcation below which either child is of the given order.
 */
class MajorAirwaysCentreLinesCleaner
{
public:
    /**
     * Constructor
     *
     * @param rMesh The centrelines mesh to clean. Note that this mesh will be edited in place.
     * @param rootIndex The root node index corresponding to the trachea.
     */
    MajorAirwaysCentreLinesCleaner(MutableMesh<1,3>& rMesh,
                                   unsigned rootIndex);

    /**
     * Trims back the mesh to the given order.
     *
     * All branches below each bifurcation with one or more children of the given order will be removed.
     *
     * @param order Horsfield order to remove branches from
     */
    void CleanUsingHorsfieldOrder(unsigned order);

    /**
     * Lengthens or shortens terminal airways to the expected length.
     *
     * Also smoothes radii.
     */
    void CleanTerminalsHueristic();

    /**
     * Removes any nodes from an airway tree that are not associated with elements.
     */
    void CleanIsolatedNodes();

private:
    /**
     * A mesh containing the major airways
     */
    MutableMesh<1,3>& mrMesh;

    /**
     * The index of the root of the airway tree (trachea)
     */
    unsigned mOutletNodeIndex;

    /**
     * Used to navigate through the airways mesh.
     */
    AirwayTreeWalker mWalker;

    /**
     * Used to calculate order on the airways mesh.
     */
    AirwayPropertiesCalculator mCalculator;

    /**
     * Maximum order below which elements will be discarded
     */
    unsigned mMaxOrder;

    /**
     * Recursively processes the given element
     *
     * @param pElement The element to recursively process
     * @param delete_me A flag to determine if this element should be deleted or not *after* processing its children
     */
    void CleanElementUsingHorsfieldOrder(Element<1,3>* pElement, bool delete_me);
};

#endif // MAJOR_AIRWAYS_CENTRE_LINES_CLEANER_HPP_
