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

#ifndef CorriasBuistICCModified_HPP_
#define CorriasBuistICCModified_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"

/**
 * This class is a modified version of the  model of a gastric
 * Interstitial Cell of Cajal.
 *
 * Reference publication is:
 *
 * Corrias A, Buist ML.
 * "Quantitative cellular description of gastric slow wave activity."
 * Am J Physiol Gastrointest Liver Physiol. 2008 Apr;294(4):G989-95. Epub 2008 Feb 14.
 *
 * Modifications include:
 * - simplified mitochondria dynamics (assumed mitochondrial potential is almost constant
 * - ability to set K+ channels-affecting CO concentrations
 * - ability to deflect a fraction of VDDR channels into the pacemaker unit.
 */
class CorriasBuistICCModified : public AbstractCardiacCell
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
        archive &  mFractionOfVDDRInPU;
        archive &  mIP3Concentration;
        archive &  mScaleFactorSerca;
        archive &  mScaleFactorCarbonMonoxide;
    }

private:

    /**fraction of VDDR channel in the PU, initialised to zero*/
    double mFractionOfVDDRInPU;
    /**the IP3 concentration, defaults to 0.0006 mM*/
    double mIP3Concentration;
    /**scales the flux through the SERCA pump.1.0-> control. 0.0-> blocked*/
    double mScaleFactorSerca;
    /**
     * Scale factor for CO-affected currents
     * Note that this the number that multiply the currents, hence it is not [CO],
     * but a function of [CO] (for example, 2.8*[CO] - 0.1)
     */
    double mScaleFactorCarbonMonoxide;

    /* Concentrations */
    double Ca_o;/**<  mM */
    double Cl_o;/**<  mM */
    double K_o;/**<  mM */
    double Na_o;/**<  mM */

    /* Nernst parameters */
    double R;/**<  pJ/nmol/K */
    double T;/**<  degK */
    double F;/**<  nC/nmol */
    double FoRT;/**<  1/mV */
    double RToF;/**<  mV */

    double  Cm ;/**<  pF */
    double Asurf_in_cm_square;/**< cm2 */
    double  Asurf ;/**<  mm2 */
    double  Cl_i  ;/**<  mM */
    double  K_i   ;/**<  mM */
    double  Na_i    ;/**<  mM */
    double  P_cyto;/**<  dim */
    double  Vol  ;/**<  mm3 */
    double  fc ;/**<  dim */
    double  fe  ;/**<  dim */
    double  fm ;/**<  dim */
    double  Q10Ca ;/**<  dim */
    double  Q10K ;/**<  dim */
    double  Q10Na  ;/**<  dim */
    double  T_exp  ;/**<  degK */

    double  G_max_BK  ;/**<  nS */
    double  G_max_CaCl  ;/**<  nS */
    double  G_max_ERG  ;/**<  nS */
    double  G_max_Ltype ;/**<  nS */
    double  G_max_NSCC  ;/**<  nS */
    double  G_max_Na ;/**<  nS */
    double  G_max_VDDR ;/**<  nS */
    double  G_max_bk    ;/**<  nS */
    double  G_max_kv11 ;/**<  nS */


    double  J_max_PMCA  ;/**<  mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms */
    double  J_max_PMCA_PU ;/**<  mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms */
    double  J_ERleak   ;/**<  1/ms    (1/s) * 1/1000 (ms/s) = 1/ms */
    double  J_max_leak ;/**<  1/ms    (1/s) * 1/1000 (ms/s) = 1/ms */
    double  Jmax_IP3  ;/**<  1/ms    (1/s) * 1/1000 (ms/s) = 1/ms */
    double  Jmax_NaCa ;/**<  mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms */
    double  Jmax_serca   ;/**<  mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms */
    double  Jmax_uni    ;/**<  1/ms    (1/s) * 1/1000 (ms/s) = 1/ms */

    double  NaPerm_o_Kperm   ;/**<  dim */
    double  L    ;/**<  dim */
    double  P_ER  ;/**<  dim */
    double  P_PU   ;/**<  dim */
    double  P_mito ;/**<  dim */
    double  b ;/**<  dim */
    double  na ;/**<  dim */

    double  K_Ca   ;/**<  mM */
    double  K_Na ;/**<  mM */
    double  K_act ;/**<  mM */
    double  K_trans   ;/**<  mM */
    double  k_serca ;/**<  mM */
    double  conc   ;/**<  mM */
    double  d_ACT   ;/**<  mM */
    double  d_IP3 ;/**<  mM */
    double  d_INH  ;/**<  mM */

    double  tau_d_CaCl;/**<  ms(s) * 1000 (ms/s) = ms */
    double  tau_d_NSCC ;/**<  ms(s) * 1000 (ms/s) = ms */
    double  tauh;/**<  ms(s ) * 1000 (ms/s) = ms */

    double  deltaPsi_B;/**<  mV */
    double  deltaPsi_star;/**<  mV */
    double  deltaPsi;/**<  mV */


     /////////////////////
     //Calculated constants
     ////////////////////
     /* Volumes */
     double V_cyto;              /**<  mm3 */
     double V_ER;                /**<  mm3 */
     double V_MITO;              /**<  mm3 */
     double V_PU;                /**<  mm3 */

     /* Temperature corrections */
     double T_correction_Ca;     /**<  dim */
     double T_correction_K;      /**<  dim */
     double T_correction_Na;     /**<  dim */
     double T_correction_BK;     /**<  uA/mm2 */

     /* Nernst potentials */
     double E_Na;                /**<  mV */
     double E_K;                 /**<  mV */
     double E_Cl;                /**<  mV */
     double E_NSCC;              /**<  mV */

     /* Activation gate time constants */
     double tau_d_ERG;           /**<  ms */
     double tau_d_Ltype;         /**<  ms */
     double tau_d_Na;            /**<  ms */
     double tau_d_VDDR;          /**<  ms */
     double tau_d_kv11;          /**<  ms */

     /* Inactivation gate time constants */
     double tau_f_Ltype;         /**<  ms */
     double tau_f_Na;            /**<  ms */
     double tau_f_VDDR;          /**<  ms */
     double tau_f_ca_Ltype;      /**<  ms */
     double tau_f_kv11;          /**<  ms */

     /* Speed ups */
     double e2FoRTdPsiMdPsiS;/**< speed-up constant*/
     double ebFoRTdPsiMdPsiS;/**< speed-up constant*/


public:
   /**
    * Constructor
    *
    * @param pSolver is a pointer to the ODE solver
    * @param pIntracellularStimulus is a pointer to the intracellular stimulus
    */
     CorriasBuistICCModified(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~CorriasBuistICCModified();

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
     * Sets the fraction of VDDR channels in the Pacemaker unit
     *
     * @param fraction the fraction of VDDr channels in the PU
     */
    void SetFractionOfVDDRInPU(double fraction);

    /**
     * Set the value of IP3 concentration in the cell
     *
     * @param concentration the concentration of IP3
     */
    void SetIP3Concentration(double concentration);

    /**
     * Set a multiplying factor for the influx of Ca2+ into the Er via the SERCA pump (set to zero will block the SW generation)
     *
     * @param scaleFactor the scale factor (=0 --> no Ca2+ uptake into the ER and, consequently, no SW)
     */
    void SetSercaPumpScaleFactor(double scaleFactor);

    /**
     * Set the carbon monoxide scale factor.
     * This will multiply the following currents: I_kv11, I_ERG, Ibk
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
CHASTE_CLASS_EXPORT(CorriasBuistICCModified)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const CorriasBuistICCModified * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, CorriasBuistICCModified * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)CorriasBuistICCModified(p_solver, p_stimulus);
        }
    }
}

#endif // CorriasBuistICCModified_HPP_
