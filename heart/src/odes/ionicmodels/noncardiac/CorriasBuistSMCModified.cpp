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


#include <cmath>
#include <cassert>
#include <memory>
#include "Exception.hpp"
#include "OdeSystemInformation.hpp"
#include "CorriasBuistSMCModified.hpp"
#include "HeartConfig.hpp"


    CorriasBuistSMCModified::CorriasBuistSMCModified(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(
                pSolver,
                14,
                0,
                pIntracellularStimulus)
    {
        mpSystemInfo = OdeSystemInformation<CorriasBuistSMCModified>::Instance();
        mScaleFactorCarbonMonoxide = 1.0; // initialise to 1 --> no effect
        mFakeIccStimulusPresent = true;//by default we use the fake ICC stimulus

        /* parameters */
        Cm = 77.0*1e-6;// 77 pF --> microF

        Asurf_in_cm_square = Cm / HeartConfig::Instance()->GetCapacitance();
        Asurf = Asurf_in_cm_square / 0.01;//cm2 --> mm2 /*cell surface area (mm2)*/

        VolCell   =      3.5e-6;      /*cell volume (mm3)*/
        hCa       =      2.014e-4;    /*conc for half inactivation of fCa */
        sCa       =      1.31e-05;    /*slope factor for inactivation of fCa */

        /* concentrations */
        Ki         =     164.0;       /*intra K conc (mM)*/
        Nai        =     10.0;        /*intra Na conc (mM)*/
        ACh        =     1e-05;       /*acetylcholine conc (mM)*/
        CaiRest    =     0.9e-04;     /*baseline Ca conc (mM)*/

        /* maximum conductances*/
        gLVA_max    =    2.33766E-05; /*max conductance of ILVA*/                 // (0.18 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gCaL_max    =    0.008441558; /*max conductance of ICaL*/                 // (65.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gBK_max     =    0.005935065; /*max conductance of IBK)*/                 // (45.7 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gKb_max     =    1.87013E-06; /*max conductance of IKb*/                  // (0.0144 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gKA_max     =    0.001168831; /*max conductance of IKA*/                  // (9.0  nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gKr_max     =    0.004545455; /*max conductance of IKr*/                  // (35.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gNa_max     =    0.00038961;  /*max conductance of INa*/                  // (3.0  nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gnsCC_max   =    0.006493506; /*max conductance of InsCC*/                // (50.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        gcouple     =     1.3*1e-6 / Asurf; /* coupling conductance bewteen fake ICC and SMC*/        // 1.3 nS * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2

        JCaExt_max   =   0.31705;     /*max flux of CaSR (mM/ms)*/

        /* Temperature corrections */
        Q10Ca        =   2.1;         /*(dim)*/
        Q10K         =   1.5;         /*(dim)*/ //1.365
        Q10Na        =   2.45;        /*(dim)*/
        Texp         =   297.0;       /*(degK)*/

         Ca_o   =         2.5;         // mM
         K_o     =        7.0;         // mM
         Na_o    =        137.0;       // mM

        /* Nernst parameters */
         R        =       8314.4;      // pJ/nmol/K
         T        =       310.0;       // degK
         F        =       96484.6;     // nC/nmol
         FoRT     =       0.03743;     // 1/mV
         RToF     =       26.7137;     // mV

        T_correct_Ca = pow(Q10Ca,(T-Texp)/10.0);/*temperature correction for Ca (dim)*/
        T_correct_K = pow(Q10K,(T-Texp)/10.0);  /*temperature correction for K (dim)*/
        T_correct_Na = pow(Q10Na,(T-Texp)/10.0);/*temperature correction for Na (dim)*/
        T_correct_gBK = gBK_max + 1.1*(T-Texp)*1e-6/Asurf; /*temperature correction for gBK*/  // (nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2

        /* Nernst potentials */
        EK = RToF*log(K_o/Ki);                  /*Nernst potential for K (mV)*/
        ENa = RToF*log(Na_o/Nai);               /*Nernst potential for Na (mV)*/
        EnsCC = -28.0;                          /*Nernst potential for nsCC (mV)*/

        Init();

    }

    CorriasBuistSMCModified::~CorriasBuistSMCModified()
    {
    }

    void CorriasBuistSMCModified::VerifyStateVariables()
    {}

    void CorriasBuistSMCModified::SetCarbonMonoxideScaleFactor(double scaleFactor)
    {
        mScaleFactorCarbonMonoxide = scaleFactor;
    }

    void CorriasBuistSMCModified::SetFakeIccStimulusPresent(bool present)
    {
        mFakeIccStimulusPresent = present;
    }

    bool CorriasBuistSMCModified::GetFakeIccStimulusPresent()
    {
        return mFakeIccStimulusPresent;
    }

    double  CorriasBuistSMCModified::GetCarbonMonoxideScaleFactor()
    {
        return mScaleFactorCarbonMonoxide;
    }

    double CorriasBuistSMCModified::GetIIonic(const std::vector<double>* pStateVariables)
    {
        if (!pStateVariables) pStateVariables = &rGetStateVariables();
        const std::vector<double>& rY = *pStateVariables;

        double ECa = 0.5*RToF*log(Ca_o/rY[13]);

        /* inward sodium current */
        double INa = gNa_max*rY[6]*rY[7]*(rY[0]-ENa);

        /* L-type calcium current */
        double ICaL = gCaL_max*rY[1]*rY[2]*rY[3]*(rY[0]-ECa);

        /* low voltage activated (T-type) calcium current */
        double ILVA = gLVA_max*rY[4]*rY[5]*(rY[0]-ECa);

        /* large conductance calcium activated potassium current */
        double Po_BK = 1.0/(1.0+exp(-(rY[0]/17.0)-2.0*log(rY[13]/0.001)));
        double IBK = T_correct_gBK*Po_BK*(rY[0]-EK);

        /* delayed rectifier potassium current */
        double IKr = mScaleFactorCarbonMonoxide*gKr_max*rY[9]*rY[10]*(rY[0]-EK);

        /* A-type potassium current */
        double IKA = mScaleFactorCarbonMonoxide*gKA_max*rY[11]*rY[12]*(rY[0]-EK);

        /* background (leakage) potassium current */
        double IKb = mScaleFactorCarbonMonoxide*gKb_max*(rY[0]-EK);

        /* non-specific cation current */
        double hCa_nsCC = 1.0/(1.0+pow((rY[13]/0.0002),-4.0));
        double rACh_nsCC = 1.0/(1.0+(0.01/ACh));
        double InsCC = gnsCC_max*rY[8]*rACh_nsCC*hCa_nsCC*(rY[0]-EnsCC);

        /* phenomenological calcium extrusion current */
        double JCaExt = JCaExt_max*pow(rY[13],1.34);

        //i_ionic_in microA/mm2
        double i_ionic = INa+ICaL+ILVA+IKr+IKA+IBK+IKb+InsCC+(JCaExt*2.0*F*VolCell/Asurf);

        assert(!std::isnan(i_ionic));
        /**
         * Now convert to microA over cm^2, the units that Chaste needs
         */
        return i_ionic / 0.01;
    }

    void CorriasBuistSMCModified::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {

        // index [0] = -69.8;           // Vm       (mV)
        // index [1] = 5.0e-6;          // d_CaL    (dim)
        // index [2] = 0.953;           // f_CaL    (dim)
        // index [3] = 1.0;             // fCa_CaL  (dim)
        // index [4] = 0.0202;          // d_LVA    (dim)
        // index [5] = 1.0;             // f_LVA    (dim)
        // index [6] = 0.0086;          // d_Na     (dim)
        // index [7] = 0.061;           // f_Na     (dim)
        // index [8] = 0.096;           // m_nsCC   (dim)
        // index [9] = 1.92e-4;         // xr1      (dim)
        // index [10] = 0.812;          // xr2      (dim)
        // index [11] = 0.00415;        // xa1      (dim)
        // index [12] = 0.716;          // xa2      (dim)
        // index [13] = 0.9e-04;        // Cai      (mM)

        double ECa = 0.5*RToF*log(Ca_o/rY[13]);

        /* inward sodium current */
        double inf_d_Na = 1.0/(1.0+exp(-(rY[0]+47.0)/4.8));
        double tau_d_Na = (0.44-0.017*rY[0])*T_correct_Na;
        double inf_f_Na = 1.0/(1.0+exp((rY[0]+78.0)/3.0));
        double tau_f_Na = (5.5-0.25*rY[0])*T_correct_Na;

        double INa = gNa_max*rY[6]*rY[7]*(rY[0]-ENa);

        /* L-type calcium current */
        double inf_d_CaL = 1.0/(1.0+exp(-(rY[0]+17.0)/4.3));
        double tau_d_CaL = 0.47*T_correct_Ca;

        double inf_f_CaL = 1.0/(1.0+exp((rY[0]+43.0)/8.9));
        double tau_f_CaL = 86.0*T_correct_Ca;

        double inf_fCa_CaL = 1.0-(1.0/(1.0+exp(-((rY[13]-CaiRest)-hCa)/sCa)));
        double tau_fCa_CaL = 2.0*T_correct_Ca;
        double ICaL = gCaL_max*rY[1]*rY[2]*rY[3]*(rY[0]-ECa);

        /* low voltage activated (T-type) calcium current */
        double inf_d_LVA = 1.0/(1.0+exp(-(rY[0]+27.5)/10.9));
        double tau_d_LVA = 3.0*T_correct_Ca;

        double inf_f_LVA = 1.0/(1.0+exp((rY[0]+15.8)/7.0));
        double tau_f_LVA = 7.58*exp(rY[0]*0.00817)*T_correct_Ca;

        double ILVA = gLVA_max*rY[4]*rY[5]*(rY[0]-ECa);

        /* large conductance calcium activated potassium current */
        double Po_BK = 1.0/(1.0+exp(-(rY[0]/17.0)-2.0*log(rY[13]/0.001)));
        double IBK = T_correct_gBK*Po_BK*(rY[0]-EK);

        /* delayed rectifier potassium current */
        double inf_xr1 = 1.0/(1.0+exp(-(rY[0]+27.0)/5.0));
        double tau_xr1 = 80.0*T_correct_K;

        double inf_xr2 = 0.2+0.8/(1.0+exp((rY[0]+58.0)/10.0));
        double tau_xr2 = (-707.0+1481.0*exp((rY[0]+36.0)/95.0))*T_correct_K;

        double IKr = mScaleFactorCarbonMonoxide*gKr_max*rY[9]*rY[10]*(rY[0]-EK);

        /* A-type potassium current */
        double inf_xa1 = 1.0/(1.0+exp(-(rY[0]+26.5)/7.9));
        double tau_xa1 = (31.8+175.0*exp(-0.5*pow(((rY[0]+44.4)/22.3),2.0)))*T_correct_K;

        double inf_xa2 = 0.1+0.9/(1.0+exp((rY[0]+65.0)/6.2));
        double tau_xa2 = 90.0*T_correct_K;

        double IKA = mScaleFactorCarbonMonoxide*gKA_max*rY[11]*rY[12]*(rY[0]-EK);

        /* background (leakage) potassium current */
        double IKb = mScaleFactorCarbonMonoxide*gKb_max*(rY[0]-EK);

        /* non-specific cation current */
        double inf_m_nsCC = 1.0/(1.0+exp(-(rY[0]+25.0)/20.0));
        double tau_m_nsCC = 150.0/(1.0+exp(-(rY[0]+66.0)/26.0));
        double hCa_nsCC = 1.0/(1.0+pow((rY[13]/0.0002),-4.0));
        double rACh_nsCC = 1.0/(1.0+(0.01/ACh));

        double InsCC = gnsCC_max*rY[8]*rACh_nsCC*hCa_nsCC*(rY[0]-EnsCC);

        /* phenomenological calcium extrusion current */
        double JCaExt = JCaExt_max*pow(rY[13],1.34);


         double t_ICCplateau = 7582.0; // time_units
         double V_decay = 37.25; // voltage_units
         double t_ICCpeak = 98.0; // time_units

         double period = 20000.0; // time_units
         double stim_start = ((time > (period * 1.0)) && (time <= (period * 2.0))) ? (period * 1.0) : ((time > (period * 2.0)) && (time <= (period * 3.0))) ? (period * 2.0) : ((time > (period * 3.0)) && (time <= (period * 4.0))) ? (period * 3.0) : ((time > (period * 4.0)) && (time <= (period * 5.0))) ? (period * 4.0) : 0.0; // time_units
         double local_time = time - (stim_start + t_ICCpeak); // time_units
         double t_ICC_stimulus = 10000.0; // time_units
         double delta_VICC = 59.0; // voltage_units

        double i_stim;
        //see whether we are running this in isolaation (and we need the fake ICC stimulus) or coupled to a real ICC model
        if (mFakeIccStimulusPresent)
        {
            //for single cell simulations where we want the fake ICC stimulus in
            i_stim = (local_time < t_ICCpeak) ? (gcouple * delta_VICC) : ((local_time >= t_ICCpeak) && (local_time <= t_ICCplateau)) ? (gcouple * delta_VICC * (1.0 / (1.0 + exp((local_time - 8000.0) / 1000.0)))) : ((local_time > t_ICCplateau) && (local_time < t_ICC_stimulus)) ? (gcouple * V_decay * (1.0 / (1.0 + exp((local_time - 8000.0) / 150.0)))) : 0.0; // current_units
        }
        else
        {
            i_stim = GetStimulus(time);//for tissue simulations with current injected into SMC
        }

        /* membrane potential */
        double Iion = INa+ICaL+ILVA+IKr+IKA+IBK+IKb+InsCC+(JCaExt*2.0*F*VolCell/Asurf);

        double voltage_derivative;
        if (mSetVoltageDerivativeToZero)
        {
            voltage_derivative = 0.0;
        }
        else
        {
            voltage_derivative = (-1.0 / 0.01) * (-i_stim + Iion);//microA/mm2---> microA/cm2
            //std::cout<<rY[0]<<std::endl;
            assert(!std::isnan(voltage_derivative));
        }

        rDY[0] =  voltage_derivative;/* Vm */
        rDY[1] = (inf_d_CaL - rY[1])/tau_d_CaL;
        rDY[2] = (inf_f_CaL - rY[2])/tau_f_CaL;
        rDY[3] = (inf_fCa_CaL - rY[3])/tau_fCa_CaL;
        rDY[4] = (inf_d_LVA - rY[4])/tau_d_LVA;
        rDY[5] = (inf_f_LVA - rY[5])/tau_f_LVA;
        rDY[6] = (inf_d_Na - rY[6])/tau_d_Na;
        rDY[7] = (inf_f_Na - rY[7])/tau_f_Na;
        rDY[8] = (inf_m_nsCC - rY[8])/tau_m_nsCC;
        rDY[9] = (inf_xr1 - rY[9])/tau_xr1;
        rDY[10] = (inf_xr2 - rY[10])/tau_xr2;
        rDY[11] = (inf_xa1 - rY[11])/tau_xa1;
        rDY[12] = (inf_xa2 - rY[12])/tau_xa2;
        rDY[13] = (-(ICaL+ILVA)*Asurf/(2.0*F*VolCell)-JCaExt); /* intracellular calcium *1000 M-> mM; /1000 F units*/
    }

template<>
void OdeSystemInformation<CorriasBuistSMCModified>::Initialise(void)
{
    // Time units: time_units
    //
    this->mSystemName = "SMC_model_Martincode";

    this->mVariableNames.push_back("Vm_SM");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(-69.8);           // Vm       (mV)

    this->mVariableNames.push_back("d_CaL");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(5.0e-6);          // d_CaL    (dim));

    this->mVariableNames.push_back("f_CaL");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.953);           // f_CaL    (dim)

    this->mVariableNames.push_back("fCa_CaL");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(1.0);             // fCa_CaL  (dim)

    this->mVariableNames.push_back("d_LVA");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.0202);          // d_LVA    (dim)

    this->mVariableNames.push_back("f_LVA");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(1.0);             // f_LVA    (dim)

    this->mVariableNames.push_back("d_Na");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.0086);          // d_Na     (dim)

    this->mVariableNames.push_back("f_Na");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.061);           // f_Na     (dim)

    this->mVariableNames.push_back("m_nsCC");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.096);           // m_nsCC   (dim)

    this->mVariableNames.push_back("xr1");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(1.92e-4);         // xr1      (dim)

    this->mVariableNames.push_back("xr2");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.812);          // xr2      (dim)

    this->mVariableNames.push_back("xa1");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.00415);        // xa1      (dim)

    this->mVariableNames.push_back("xa2");
    this->mVariableUnits.push_back("dim");
    this->mInitialConditions.push_back(0.716);          // xa2      (dim)

    this->mVariableNames.push_back("Cai");
    this->mVariableUnits.push_back("mM");
    this->mInitialConditions.push_back(0.9e-04);        // Cai      (mM)

    this->mInitialised = true;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CorriasBuistSMCModified)
