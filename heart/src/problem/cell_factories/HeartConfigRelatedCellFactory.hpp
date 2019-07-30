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

#ifndef HEARTCONFIGRELATEDCELLFACTORY_HPP_
#define HEARTCONFIGRELATEDCELLFACTORY_HPP_

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#include "AbstractCardiacCellFactory.hpp"
#include "HeartConfig.hpp"

#include "AbstractStimulusFunction.hpp"
#include "MultiStimulus.hpp" // Included here for archiving - see below.
#include "SimpleStimulus.hpp"

#include "ChasteCuboid.hpp"
#include "AbstractChasteRegion.hpp"

#include "DynamicCellModelLoader.hpp"

/*
 * Even though these classes are only used in the .cpp file, they need to be
 * included here for serialization to work - the archiving code needs to see
 * the CHASTE_CLASS_EXPORT incantations.
 */
#include "DiFrancescoNoble1985.hpp"
#include "FaberRudy2000.hpp"
#include "FaberRudy2000Opt.hpp"
#include "FoxModel2002.hpp"
#include "FoxModel2002BackwardEuler.hpp"
#include "HodgkinHuxley1952.hpp"
#include "LuoRudy1991.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "Mahajan2008.hpp"
#include "Mahajan2008BackwardEuler.hpp"
#include "Maleckar2008.hpp"
#include "TenTusscher2006Epi.hpp"
#include "TenTusscher2006EpiBackwardEuler.hpp"

/*
 * The following cell models aren't used by this class, but do exist in Chaste.
 * So users of the source code could see issues with archiving if they aren't included.
 */
#include "FoxModel2002Opt.hpp"
#include "LuoRudy1991Opt.hpp"
#include "Mahajan2008Opt.hpp"
#include "Maleckar2008Opt.hpp"
#include "NobleVargheseKohlNoble1998a.hpp"
#include "NobleVargheseKohlNoble1998aOpt.hpp"
#include "NobleVargheseKohlNoble1998aBackwardEuler.hpp"
#include "Shannon2004.hpp"
#include "TenTusscher2006EpiOpt.hpp"

/**
 * This is a cardiac cell factory which uses the settings from HeartConfig to set up the cells.
 *
 * It thus supports such features as heterogeneities (in both the type of cell, and for some cells,
 * cell model parameters), and the ability to dynamically convert CellML files into C++ code,
 * compile this, and load the resulting model.
 */
template<unsigned SPACE_DIM>
class HeartConfigRelatedCellFactory : public AbstractCardiacCellFactory<SPACE_DIM>
{
private:
    /** Default cardiac cell model to be used in all tissue (except heterogeneous regions)*/
    cp::ionic_model_selection_type mDefaultIonicModel;
    /** List of axis-aligned box regions which contain heterogeneous cardiac ionic model types*/
    std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > mIonicModelRegions;
    /** List of ionic model (size matches that of mIonicModelRegions)*/
    std::vector<cp::ionic_model_selection_type> mIonicModelsDefined;

    /** List of axis-aligned box or ellipsoid regions which represent areas to stimulate*/
    std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > mStimulatedAreas;
    /** List of intracellular current stimuli to apply (size matches that of mStimulatedAreas)*/
    std::vector<boost::shared_ptr<AbstractStimulusFunction> > mStimuliApplied;

    /**
     *  List of regions which represent areas in which to give parametric heterogeneity (scaling gating parameters)
     *  This vector will be filled in by the HeartConfig::GetCellHeterogeneity method if the user requested
     *  to specify the heterogeneity areas by cuboids, or, alternatively, by the FillInCellularTransmuralAreas method
     *  in the cell factory called by the problem class AFTER setting the mesh
     *  (which is needed for the calculations of the distance maps for the calculations of heterogeneities).
     *
     *  When creating a cardiac cell for each node (CreateCardiacCellForTissueNode) the code will check whether
     *  that node is contained in the heterogeneity area or not.
     *
     */
    std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > mCellHeterogeneityAreas;
    /** List of scale factors for Gks scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorGks;
    /** List of scale factors for Ito scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorIto;
    /** List of scale factors for Gkr scaling in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<double> mScaleFactorGkr;

    /** Named parameters to be set in each region (size of list matches that of mCellHeterogeneityAreas)*/
    std::vector<std::map<std::string, double> > mParameterSettings;

    /**
     * Called by the constructor to convert any CellML files used as dynamically loaded cell models to shared libraries.
     * This is necessary since the conversion process must be done collectively.
     *
     * @note Must be called collectively.
     */
    void PreconvertCellmlFiles();

    /**
     * @return a loader for the given (dynamically loadable) cell model.
     * @param rModel  model to load
     * @param isCollective  whether we are being called collectively
     */
    DynamicCellModelLoaderPtr LoadDynamicModel(const cp::ionic_model_selection_type& rModel,
                                               bool isCollective);

public:
    /**
     * Constructor reads settings from the configuration file.
     *
     * @note Must be called collectively.
     */
    HeartConfigRelatedCellFactory();

    /** Destructor*/
    ~HeartConfigRelatedCellFactory();

    /**
     * @return a newly created correct tissue cell for a given region in the mesh
     * @param intracellularStimulus is computed in CreateCardiacCellForTissueNode determined by the list of stimulation regions
     * @param nodeIndex is the global index within the mesh
     */
    AbstractCardiacCellInterface* CreateCellWithIntracellularStimulus(
            boost::shared_ptr<AbstractStimulusFunction> intracellularStimulus,
            unsigned nodeIndex);

    /**
     * @return a newly created correct stimulated tissue cell for a given region in the mesh
     * The stimulus is determined in this method (using the list of stimulation regions).
     * The cardiac cell type (and parameters) are
     * determined in the CreateCellWithIntracellularStimulus method
     * @param pNode pointer to the node.
     */
    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<SPACE_DIM>* pNode);

    /**
     * Helper method to calculate and fill in the heterogeneities areas (mCellHeterogeneityAreas)
     */
    void FillInCellularTransmuralAreas();

    /**
     * Set parameters of a cell based on heterogeneity definitions in HeartConfig.
     *
     * @param pCell  the cell to set parameters of
     * @param nodeIndex  the index of the node corresponding to this cell in the mesh
     */
    void SetCellParameters(AbstractCardiacCellInterface* pCell, unsigned nodeIndex);

    /**
     * Set the intracellular stimulus for a cell based on the definitions in HeartConfig.
     *
     * @param pCell  the cell to set stimulus of
     * @param nodeIndex  the index of the node corresponding to this cell in the mesh
     */
    void SetCellIntracellularStimulus(AbstractCardiacCellInterface* pCell, unsigned nodeIndex);
};


#endif /*HEARTCONFIGRELATEDCELLFACTORY_HPP_*/
