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

#include "CorriasBuistICCModified.hpp"
#include <cmath>
#include <cassert>
#include <memory>
#include "Exception.hpp"
#include "OdeSystemInformation.hpp"
#include "HeartConfig.hpp"

    CorriasBuistICCModified::CorriasBuistICCModified(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(
                pSolver,
                18,//reduced by 3 from original
                0,
                pIntracellularStimulus)
    {

        mpSystemInfo = OdeSystemInformation<CorriasBuistICCModified>::Instance();
        mFractionOfVDDRInPU = 0.0;
        mIP3Concentration = 0.0006;
        mScaleFactorSerca = 1.0;
        mScaleFactorCarbonMonoxide = 1.0; //initialise to 1 --> no effect
        //IP3  =  0.00065;// mM  *** no longer used ***

        /////////////////////
        //Constants
        ////////////////////

        /* Concentrations */
        Ca_o = 2.5    ;// mM
        Cl_o  =134.0  ;// mM
        K_o   =7.0    ;// mM
        Na_o  =137.0  ;// mM

        /* Nernst parameters */
        R =    8314.4 ;// pJ/nmol/K
        T =    310.0  ;// degK
        F =    96484.6;// nC/nmol
        FoRT   =    0.03743;// 1/mV
        RToF   =    26.7137;// mV

        Cm = 25.0*1e-6;// 25 pF --> microF

        Asurf_in_cm_square = Cm / HeartConfig::Instance()->GetCapacitance();
        Asurf = Asurf_in_cm_square / 0.01;//cm2 --> mm2

        Cl_i = 88.0   ;// mM
        K_i  = 120.0  ;// mM
        Na_i =  30.0   ;// mM
        P_cyto = 0.7;// dim
        Vol  =  1.0e-6 ;// mm3
        fc = 0.01   ;// dim
        fe = 0.01   ;// dim
        fm = 0.0003 ;// dim
        Q10Ca = 2.1;// dim
        Q10K  = 1.5;// dim
        Q10Na = 2.45   ;// dim
        T_exp = 297.0  ;// degK

        G_max_BK     =   23.0 * 1e-6 / Asurf;//  9.2e-3  mS/mm2           (23.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_CaCl   =   10.1 * 1e-6 / Asurf; //4.04e-3  mS/mm2           (10.1 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_ERG    =   2.5 * 1e-6 / Asurf; //1.0e-3   mS/mm2           ( 2.5 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_Ltype  =   2.0 * 1e-6 / Asurf;//0.8e-3    mS/mm2           ( 2.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_NSCC   =   12.15 * 1e-6 / Asurf;//4.86e-3 mS/mm2           (12.15nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_Na     =   20.0 * 1e-6 / Asurf;//8.0e-3  mS/mm2           (20.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_VDDR   =   3.0 * 1e-6 / Asurf;//1.2e-3   mS/mm2           ( 3.0 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_bk     =   0.15 * 1e-6 / Asurf;//0.06e-3 mS/mm2           (0.15 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2
        G_max_kv11   =   6.3 * 1e-6 / Asurf;//2.52e-3  mS/mm2           ( 6.3 nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2

        J_max_PMCA = 0.088464e-3 ;// mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms
        J_max_PMCA_PU =  0.33e-3;// mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms
        J_ERleak  = 1.666667e-3 ;// 1/ms    (1/s) * 1/1000 (ms/s) = 1/ms
        J_max_leak = 0.0;// 1/ms    (1/s) * 1/1000 (ms/s) = 1/ms
        Jmax_IP3  =  50000.0e-3  ;// 1/ms    (1/s) * 1/1000 (ms/s) = 1/ms
        Jmax_NaCa = 0.05e-3;// mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms
        Jmax_serca = 1.8333e-3   ;// mM/ms   (mM/s) * 1/1000 (s/ms) = mM/ms
        Jmax_uni  = 5000.0e-3   ;// 1/ms    (1/s) * 1/1000 (ms/s) = 1/ms

        NaPerm_o_Kperm = 1.056075    ;// dim
        L = 50.0   ;// dim
        P_ER  = 0.1;// dim
        P_PU =  0.001  ;// dim
        P_mito = 0.12871;// dim
        b = 0.5;// dim
        na = 2.8;// dim

        K_Ca  = 0.003  ;// mM
        K_Na  = 9.4;// mM
        K_act = 0.00038;// mM
        K_trans = 0.006  ;// mM
        k_serca = 0.00042;// mM
        conc  = 0.001  ;// mM
        d_ACT = 0.001  ;// mM
        d_IP3 = 0.00025;// mM
        d_INH = 0.0014 ;// mM

        tau_d_CaCl = 0.03e3 ;// ms(s) * 1000 (ms/s) = ms
        tau_d_NSCC = 0.35e3 ;// ms(s) * 1000 (ms/s) = ms
        tauh  = 4.0e3  ;// ms(s ) * 1000 (ms/s) = ms

        deltaPsi_B = 50.0   ;// mV
        deltaPsi_star =  91.0   ;// mV
        deltaPsi  = 164.000044  ;// mV

         /////////////////////
         //Calculated constants
         ////////////////////

        /* Volumes */
        V_cyto = Vol*P_cyto;
        V_MITO = Vol*P_mito;
        V_PU = Vol*P_PU;
        V_ER = Vol*P_ER;

        /* Temperature corrections */
        T_correction_Ca = pow(Q10Ca, (T-T_exp)/10.0);
        T_correction_K = pow(Q10K, (T-T_exp)/10.0);
        T_correction_Na = pow(Q10Na, (T-T_exp)/10.0);
        T_correction_BK = 1.1*(T-T_exp)*1e-6/Asurf;  //(nS) * 1e-6 (mS/nS) / Asurf (mm2) = mS/mm2

        /* Nernst potentials */
        E_Na = RToF*log(Na_o/Na_i);
        E_K = RToF*log(K_o/K_i);
        E_Cl = RToF*log(Cl_i/Cl_o);
        E_NSCC = RToF*log((K_o+Na_o*NaPerm_o_Kperm)/(K_i+Na_i*NaPerm_o_Kperm));

        /* Activation gate time constants s->ms */
        tau_d_ERG = T_correction_K*0.003*1000.0;
        tau_d_Ltype = T_correction_Ca*0.001*1000.0;
        tau_d_Na = T_correction_Na*0.003*1000.0;
        tau_d_VDDR = T_correction_Ca*0.006*1000.0;
        tau_d_kv11 = T_correction_K*0.005*1000.0;

        /* Inactivation gate time constants s->ms */
        tau_f_Ltype = T_correction_Ca*0.086*1000.0;
        tau_f_Na = T_correction_Na*0.0016*1000.0;
        tau_f_VDDR = T_correction_Ca*0.04*1000.0;
        tau_f_ca_Ltype = T_correction_Ca*0.002*1000.0;
        tau_f_kv11 = T_correction_K*0.005*1000.0;

        /* Speed ups */
        e2FoRTdPsiMdPsiS = exp(-2.0*FoRT*(deltaPsi-deltaPsi_star));
        ebFoRTdPsiMdPsiS = exp(b*FoRT*(deltaPsi-deltaPsi_star));


        Init();

    }

    CorriasBuistICCModified::~CorriasBuistICCModified()
    {
    }

    void CorriasBuistICCModified::VerifyStateVariables()
    {}

    void CorriasBuistICCModified::SetSercaPumpScaleFactor(double scaleFactor)
    {
        mScaleFactorSerca = scaleFactor;
    }

    void CorriasBuistICCModified::SetFractionOfVDDRInPU(double fraction)
    {
        mFractionOfVDDRInPU = fraction;
    }

    void CorriasBuistICCModified::SetIP3Concentration(double concentration)
    {
        mIP3Concentration = concentration;
    }

    void CorriasBuistICCModified::SetCarbonMonoxideScaleFactor(double scaleFactor)
    {
        mScaleFactorCarbonMonoxide = scaleFactor;
    }

    double  CorriasBuistICCModified::GetCarbonMonoxideScaleFactor()
    {
        return mScaleFactorCarbonMonoxide;
    }

    double CorriasBuistICCModified::GetIIonic(const std::vector<double>* pStateVariables)
    {
        if (!pStateVariables) pStateVariables = &rGetStateVariables();
        const std::vector<double>& rY = *pStateVariables;

        // index 0:  Vm         (mV)
        // index 1:  Ca_i       (mM)
        // index 2:  Ca_ER      (mM)
        // index 3:  Ca_PU      (mM)
        // index 4:  Ca_m       (mM)
        // index 5:  h          (dim)
        // index 6:  d_CaCl     (dim)
        // index 7:  d_ERG      (dim)
        // index 8:  d_Ltype    (dim)
        // index 9:  d_NSCC     (dim)
        // index 10: d_Na       (dim)
        // index 11: d_VDDR     (dim)
        // index 12: d_kv11     (dim)
        // index 13: f_Ltype    (dim)
        // index 14: f_Na       (dim)
        // index 15: f_VDDR     (dim)
        // index 16: f_ca_Ltype (dim)
        // index 17: f_kv11     (dim)

        double E_Ca = 0.5*RToF*log(Ca_o/rY[1]);
        /* --- INa --- */
        double I_Na = G_max_Na*rY[14]*rY[10]*(rY[0]-E_Na);
        /* --- ILtype --- */
        double I_Ltype = G_max_Ltype*rY[13]*rY[8]*rY[16]*(rY[0]-E_Ca);
        /* --- IVDDR --- */
        double I_VDDR = G_max_VDDR*rY[15]*rY[11]*(rY[0]-E_Ca);
        /* --- IKv1.1 --- */
        double I_kv11 = mScaleFactorCarbonMonoxide*G_max_kv11*rY[17]*rY[12]*(rY[0]-E_K);
        /* --- IERG --- */
        double I_ERG = mScaleFactorCarbonMonoxide*G_max_ERG*rY[7]*(rY[0]-E_K);
        /* --- IBK --- */
        double d_BK = 1.0/(1.0+((exp(rY[0]/-17.0))/((rY[1]/0.001)*(rY[1]/0.001))));
        double I_BK = (G_max_BK+T_correction_BK)*d_BK*(rY[0]-E_K);
        /* --- IKb --- */
        double I_bk = mScaleFactorCarbonMonoxide*G_max_bk*(rY[0]-E_K);
        /* --- ICaCL --- */
        double I_CaCl = G_max_CaCl*rY[6]*(rY[0]-E_Cl);
        /* --- INSCC --- */
        double I_NSCC = G_max_NSCC*rY[9]*(rY[0]-E_NSCC);
        /* --- JpmCa --- */
        double J_PMCA = J_max_PMCA*1.0/(1.0+(0.000298/rY[1]));

        //i_ionic_in microA/mm2
        double i_ionic = (I_Na+I_Ltype+I_VDDR+I_kv11+I_ERG+I_BK+I_CaCl+I_NSCC+I_bk+(J_PMCA*2.0*F*V_cyto/Asurf));
        assert(!std::isnan(i_ionic));
        /**
         * Now convert to microA over cm^2, the units that Chaste needs
         */
        return i_ionic / 0.01;
    }

    void CorriasBuistICCModified::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        // index 0:  Vm         (mV)
        // index 1:  Ca_i       (mM)
        // index 2:  Ca_ER      (mM)
        // index 3:  Ca_PU      (mM)
        // index 4:  Ca_m       (mM)
        // index 5:  h          (dim)
        // index 6:  d_CaCl     (dim)
        // index 7:  d_ERG      (dim)
        // index 8:  d_Ltype    (dim)
        // index 9:  d_NSCC     (dim)
        // index 10: d_Na       (dim)
        // index 11: d_VDDR     (dim)
        // index 12: d_kv11     (dim)
        // index 13: f_Ltype    (dim)
        // index 14: f_Na       (dim)
        // index 15: f_VDDR     (dim)
        // index 16: f_ca_Ltype (dim)
        // index 17: f_kv11     (dim)

        /* ----------------- */
        /* Membrane currents */
        /* ----------------- */

        double E_Ca = 0.5*RToF*log(Ca_o/rY[1]);

        /* --- INa --- */
        double d_inf_Na = 1.0/(1.0+exp((rY[0]+47.0)/-4.8));
        double f_inf_Na = 1.0/(1.0+exp((rY[0]+78.0)/7.0));
        double I_Na = G_max_Na*rY[14]*rY[10]*(rY[0]-E_Na);

        /* --- ILtype --- */
        double d_inf_Ltype = 1.0/(1.0+exp((rY[0]+17.0)/-4.3));
        double f_inf_Ltype = 1.0/(1.0+exp((rY[0]+43.0)/8.9));
        double f_ca_inf_Ltype = 1.0-1.0/(1.0+exp((rY[1]-0.0001-0.000214)/-0.0000131));
        double I_Ltype = G_max_Ltype*rY[13]*rY[8]*rY[16]*(rY[0]-E_Ca);

        /* --- IVDDR --- */
        double d_inf_VDDR = 1.0/(1.0+exp((rY[0]+26.0)/-6.0));
        double f_inf_VDDR = 1.0/(1.0+exp((rY[0]+66.0)/6.0));
        double I_VDDR = G_max_VDDR*rY[15]*rY[11]*(rY[0]-E_Ca);

        /* --- IKv1.1 --- */
        double d_inf_kv11 = 1.0/(1.0+exp((rY[0]+25.0)/-7.7));
        double f_inf_kv11 = 0.5+0.5/(1.0+exp((rY[0]+44.8)/4.4));
        double I_kv11 = mScaleFactorCarbonMonoxide*G_max_kv11*rY[17]*rY[12]*(rY[0]-E_K);

        /* --- IERG --- */
        double d_inf_ERG = 0.2+0.8/(1.0+exp((rY[0]+20.0)/-1.8));
        double I_ERG = mScaleFactorCarbonMonoxide*G_max_ERG*rY[7]*(rY[0]-E_K);

        /* --- IBK --- */
        //LUT d_BK = 1.0/(1.0+exp((rY[0]/-17.0)-2.0*log(rY[1]/0.001)));
        double d_BK = 1.0/(1.0+((exp(rY[0]/-17.0))/((rY[1]/0.001)*(rY[1]/0.001))));
        double I_BK = (G_max_BK+T_correction_BK)*d_BK*(rY[0]-E_K);

        /* --- IKb --- */
        double I_bk = mScaleFactorCarbonMonoxide*G_max_bk*(rY[0]-E_K);

        /* --- ICaCL --- */
        double tmp1 = 0.00014/rY[1];
        double d_inf_CaCl = 1.0/(1.0+(tmp1*tmp1*tmp1));
        double I_CaCl = G_max_CaCl*rY[6]*(rY[0]-E_Cl);

        /* --- INSCC --- */
        double d_inf_NSCC = 1.0/(1.0+pow(0.0000745/rY[3], -85.0));
        double I_NSCC = G_max_NSCC*rY[9]*(rY[0]-E_NSCC);

        /* --- JpmCa --- */
        double J_PMCA = J_max_PMCA*1.0/(1.0+(0.000298/rY[1]));

        /* ----------------- */
        /* ER fluxes         */
        /* ----------------- */

        //tmp1 = IP3/(IP3+d_IP3);
        tmp1 = mIP3Concentration/(mIP3Concentration+d_IP3);
        double tmp2 = rY[3]/(rY[3]+d_ACT);
        double J_ERout = (Jmax_IP3*tmp1*tmp1*tmp1*tmp2*tmp2*tmp2*rY[5]*rY[5]*rY[5]+J_ERleak)*(rY[2]-rY[3]);
        double J_SERCA = mScaleFactorSerca*Jmax_serca*rY[3]*rY[3]/(k_serca*k_serca+rY[3]*rY[3]);

        /* ----------------- */
        /* Mito fluxes       */
        /* ----------------- */

        /* Uniporter */
        tmp1 = 1.0+rY[3]/K_trans;
        double MWC = conc*(rY[3]/K_trans)*tmp1*tmp1*tmp1/(tmp1*tmp1*tmp1*tmp1+L/pow(1.0+rY[3]/K_act, na));
        double J_uni = Jmax_uni*(MWC-rY[4]*e2FoRTdPsiMdPsiS)*2.0*FoRT*(deltaPsi-deltaPsi_star)/(1.0-e2FoRTdPsiMdPsiS);

        /* NaCa Exchanger */
        double J_NaCa = Jmax_NaCa*ebFoRTdPsiMdPsiS/((1.0+K_Na*K_Na/(Na_i*Na_i))*(1.0+K_Ca/rY[4]));

        /* ----------------- */
        /* Cyto fluxes       */
        /* ----------------- */

        double J_leak = J_max_leak*(rY[3]-rY[1]); /* P.U.->Cai */

        /* ----------------- */
        /* Entrainment       */
        /* ----------------- */

        double E_Ca_PU = 0.5*RToF*log(Ca_o/rY[3]);
        double I_VDDR_PU = G_max_VDDR*rY[11]*rY[15]*(rY[0]-E_Ca_PU);
        //J_PMCA_PU = J_max_PMCA_PU*1.0/(1.0+exp(-(rY[3]-0.0001)/0.000015));
        double J_PMCA_PU = J_max_PMCA_PU*1.0/(1.0+exp(-(rY[3]-0.0001)/0.000015));

        double i_stim = GetStimulus(time);

        /* -------------------- */
        /* Resting Membrane, CO */
        /* -------------------- */

        //tmp1 = 2.8*spatVar[1]-0.1;
        //I_kv11 = I_kv11*tmp1;
        //I_ERG = I_ERG*tmp1;
        //I_bk = I_bk*tmp1;
        double voltage_derivative;
        if (mSetVoltageDerivativeToZero)
        {
            voltage_derivative = 0.0;
        }
        else
        {
            voltage_derivative = (-1.0 / 0.01) * (i_stim + I_Na+I_Ltype+I_VDDR+I_kv11+I_ERG+I_BK+I_CaCl+I_NSCC+I_bk+(J_PMCA*2.0*F*V_cyto/Asurf));
            assert(!std::isnan(voltage_derivative));
        }

        rDY[0] =  voltage_derivative;/* Vm */
        rDY[1] = fc*((-I_Ltype-I_VDDR)*Asurf/(2.0*F*V_cyto)+J_leak-J_PMCA);
        rDY[2] = fe*(J_SERCA-J_ERout);
        rDY[3] = fc*((J_NaCa-J_uni)*V_MITO/V_PU+(J_ERout-J_SERCA)*V_ER/V_PU-J_leak*V_cyto/V_PU);
        rDY[3]-= fc*(((mFractionOfVDDRInPU*I_VDDR_PU*Asurf)/(2.0*F*V_PU))+J_PMCA_PU); // *** new, 4% IVDDR ***
        rDY[4] = fm*(J_uni-J_NaCa);
        rDY[5] = 1.0*(d_INH-rY[5]*(rY[3]+d_INH))/tauh;
        rDY[6] = (d_inf_CaCl-rY[6])/tau_d_CaCl;
        rDY[7] = (d_inf_ERG-rY[7])/tau_d_ERG;
        rDY[8] = (d_inf_Ltype-rY[8])/tau_d_Ltype;
        rDY[9] = (d_inf_NSCC-rY[9])/tau_d_NSCC;
        rDY[10] = (d_inf_Na-rY[10])/tau_d_Na;
        rDY[11] = (d_inf_VDDR-rY[11])/tau_d_VDDR;
        rDY[12] = (d_inf_kv11-rY[12])/tau_d_kv11;
        rDY[13] = (f_inf_Ltype-rY[13])/tau_f_Ltype;
        rDY[14] = (f_inf_Na-rY[14])/tau_f_Na;
        rDY[15] = (f_inf_VDDR-rY[15])/tau_f_VDDR;
        rDY[16] = (f_ca_inf_Ltype-rY[16])/tau_f_ca_Ltype;
        rDY[17] = (f_inf_kv11-rY[17])/tau_f_kv11;
    }

template<>
void OdeSystemInformation<CorriasBuistICCModified>::Initialise(void)
{

    this->mSystemName = "ICC_model_Martincode";

    this->mVariableNames.push_back("Vm");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(-67.53988);         // Vm         (mV);

    this->mVariableNames.push_back("Ca_i");
    this->mVariableUnits.push_back("mM");
    this->mInitialConditions.push_back(0.00001);           // Ca_i       (mM));

    this->mVariableNames.push_back("Ca_ER");
    this->mVariableUnits.push_back("mM");
    this->mInitialConditions.push_back(0.00695);           // Ca_ER      (mM));

    this->mVariableNames.push_back("Ca_PU");
    this->mVariableUnits.push_back("mM");
    this->mInitialConditions.push_back(0.000095);          // Ca_PU      (mM));

    this->mVariableNames.push_back("Ca_m");
    this->mVariableUnits.push_back("mM");
    this->mInitialConditions.push_back(0.000138);          // Ca_m       (mM));

    this->mVariableNames.push_back("h");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.939443);          // h          (dim));

    this->mVariableNames.push_back("d_CaCl");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00038);           // d_CaCl     (dim));

    this->mVariableNames.push_back("d_ERG");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.2);               // d_ERG      (dim));

    this->mVariableNames.push_back("d_Ltype");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.000008);          // d_Ltype    (dim));

    this->mVariableNames.push_back("d_NSCC");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);               // d_NSCC     (dim));

    this->mVariableNames.push_back("d_Na");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.013778);         // d_Na       (dim));

    this->mVariableNames.push_back("d_VDDR");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00099);          // d_VDDR     (dim));

    this->mVariableNames.push_back("d_kv11");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.003992);         // d_kv11     (dim));

    this->mVariableNames.push_back("f_Ltype");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back( 0.940072);         // f_Ltype    (dim));

    this->mVariableNames.push_back("f_Na");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.182426);         // f_Na       (dim));

    this->mVariableNames.push_back("f_VDDR");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.562177);         // f_VDDR     (dim));

    this->mVariableNames.push_back("f_ca_Ltype");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);              // f_ca_Ltype (dim));

    this->mVariableNames.push_back("f_kv11");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.997143);         // f_kv11     (dim));

    this->mInitialised = true;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CorriasBuistICCModified)
