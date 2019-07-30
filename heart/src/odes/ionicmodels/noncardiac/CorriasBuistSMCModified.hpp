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

#ifndef CorriasBuistSMCModified_HPP_
#define CorriasBuistSMCModified_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
/**
 * This class is a modified version of the model of a gastric
 * Smooth Muscle Cell.
 *
 * Reference publication is:
 *
 * Corrias A, Buist ML.
 * "A quantitative model of gastric smooth muscle cellular activation."
 * Ann Biomed Eng. 2007 Sep;35(9):1595-607. Epub 2007 May 8.
 *
 * Modifications include:
 * - ability to include/exclude built-in fake ICC stimulus
 * - ability to set K+ channels-affecting CO concentrations
 */
class CorriasBuistSMCModified : public AbstractCardiacCell
{
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacCell >(*this);
    }

private:

    /**
     * Scale factor for CO-affected currents
     * Note that this the number that multiply the currents, hence it is not [CO],
     * but a function of [CO] (for example, 2.8*[CO] - 0.1)
     */
    double mScaleFactorCarbonMonoxide;

    /**
     * True if the fake built-in ICC stimulus is present
     */
    bool mFakeIccStimulusPresent;

    double Cm;/**< membrane capacitance, pF*/

    double Asurf_in_cm_square;/**< Surface area in cm^2*/
    double Asurf;/**< surface area (mm^2)*/

    double VolCell;/**< cell volume (mm^3)*/
    double hCa;/**< conc for half inactivation of fCa */
    double sCa;/**< lope factor for inactivation of fCa */

    /* concentrations */
    double Ki;       /**< intra K conc (mM)*/
    double Nai;        /**< intra Na conc (mM)*/
    double ACh;       /**< acetylcholine conc (mM)*/
    double CaiRest;     /**< baseline Ca conc (mM)*/

    /* maximum conductances*/
    double gLVA_max; /**< max conductance of ILVA*/                 // (0.18 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gCaL_max; /**< max conductance of ICaL*/                 // (65.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gBK_max; /**< max conductance of IBK)*/                 // (45.7 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gKb_max; /**< max conductance of IKb*/                  // (0.0144 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gKA_max; /**< max conductance of IKA*/                  // (9.0  nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gKr_max; /**< max conductance of IKr*/                  // (35.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gNa_max;  /**< max conductance of INa*/                  // (3.0  nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gnsCC_max; /**< max conductance of InsCC*/                // (50.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double gcouple;   /**<  coupling conductance bewteen fake ICC and SMC*/        // 1.3 nS * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
    double JCaExt_max;     /**< max flux of CaSR (mM/ms)*/

    /* Temperature corrections */
    double Q10Ca;         /**< (dim)*/
    double Q10K;         /**< (dim)*/ //1.365
    double Q10Na;        /**< (dim)*/
    double Texp;       /**< (degK)*/

    double T_correct_Ca ;/**< temperature correction for Ca (dim)*/
    double T_correct_K ;  /**< temperature correction for K (dim)*/
    double T_correct_Na;/**< temperature correction for Na (dim)*/
    double T_correct_gBK; /**< temperature correction for gBK*/  // (nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2

    /* Nernst potentials */
    double EK;                  /**< Nernst potential for K (mV)*/
    double ENa ;               /**< Nernst potential for Na (mV)*/
    double EnsCC;                          /**< Nernst potential for nsCC (mV)*/

    double Ca_o;         /**<  mM */
    double K_o;          /**< mM */
    double Na_o;         /**< mM */

    /* Nernst parameters */
    double R;    /**<  pJ/nmol/K*/
    double T;       /**<  degK*/
    double F;     /**<  nC/nmol*/
    double FoRT;     /**<  1/mV*/
    double RToF;    /**<  mV*/

public:

    /**
     * Constructor
     *
     * @param pSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    CorriasBuistSMCModified(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~CorriasBuistSMCModified();

    /**
     * Now empty
     */
    void VerifyStateVariables();

    /**
     * Calculates the ionic current
     *
     * @param pStateVariables the state variables of this model
     * @return the total ionic current
     */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);

    /**
     * Compute the RHS of the FitHugh-Nagumo system of ODEs
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);


    /**
     * Set whether we want the fake ICC stimulus or not.
     * It changes the member variable mFakeIccStimulusPresent (which is true by default).
     *
     * @param present - true if we want the fake ICC stimulus, false otherwise
     */
    void SetFakeIccStimulusPresent(bool present);

    /**
     * @return true if the fake ICC stimulus is present
     */
    bool GetFakeIccStimulusPresent();

    /**
     * @return the Carbon Monoxide scale for
     */
    double SetCarbonMonoxideScaleFactor();

    /**
     * Set the carbon monoxide scale factor.
     * This will multiply the following currents: I_kr, I_Ka, Ibk
     *
     * @param scaleFactor the scale factor that multiply the currents.
     */
    void SetCarbonMonoxideScaleFactor(double scaleFactor);

    /**
     * @return the Carbon Monoxide scale factor
     */
    double GetCarbonMonoxideScaleFactor();
};

// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CorriasBuistSMCModified)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const CorriasBuistSMCModified * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, CorriasBuistSMCModified * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)CorriasBuistSMCModified(p_solver, p_stimulus);
        }
    }
}

#endif // CorriasBuistSMCModified_HPP_
