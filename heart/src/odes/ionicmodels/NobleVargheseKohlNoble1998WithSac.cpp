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

#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "HeartConfig.hpp"

CML_noble_varghese_kohl_noble_1998_basic_with_sac::CML_noble_varghese_kohl_noble_1998_basic_with_sac(
        boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCell(pSolver, 22, 0, pIntracellularStimulus)
{
    mStretch = 1.0;

    mpSystemInfo = OdeSystemInformation<CML_noble_varghese_kohl_noble_1998_basic_with_sac>::Instance();

    Init();

}

CML_noble_varghese_kohl_noble_1998_basic_with_sac::~CML_noble_varghese_kohl_noble_1998_basic_with_sac(void)
{
}

double CML_noble_varghese_kohl_noble_1998_basic_with_sac::GetIIonic(const std::vector<double>* pStateVariables)
{
    if (!pStateVariables) pStateVariables = &rGetStateVariables();
    const std::vector<double>& rY = *pStateVariables;
    double var_membrane__V = rY[0];
    // Units: millivolt; Initial value: -92.849333
    double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
    // Units: dimensionless; Initial value: 1.03e-5
    double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
    // Units: dimensionless; Initial value: 2e-7
    double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
    // Units: dimensionless; Initial value: 0.001302
    double var_fast_sodium_current_m_gate__m = rY[4];
    // Units: dimensionless; Initial value: 0.0016203
    double var_fast_sodium_current_h_gate__h = rY[5];
    // Units: dimensionless; Initial value: 0.9944036
    double var_L_type_Ca_channel_d_gate__d = rY[6];
    // Units: dimensionless; Initial value: 0
    double var_L_type_Ca_channel_f_gate__f = rY[7];
    // Units: dimensionless; Initial value: 1
    double var_L_type_Ca_channel_f2_gate__f2 = rY[8];
    // Units: dimensionless; Initial value: 0.9349197
    double var_L_type_Ca_channel_f2ds_gate__f2ds = rY[9];
    // Units: dimensionless; Initial value: 0.9651958
    double var_transient_outward_current_s_gate__s = rY[10];
    // Units: dimensionless; Initial value: 0.9948645
    double var_transient_outward_current_r_gate__r = rY[11];
    // Units: dimensionless; Initial value: 0
    double var_intracellular_sodium_concentration__Na_i = rY[14];
    // Units: millimolar; Initial value: 7.3321223
    double var_intracellular_potassium_concentration__K_i = rY[15];
    // Units: millimolar; Initial value: 136.5644281
    double var_intracellular_calcium_concentration__Ca_i = rY[16];
    // Units: millimolar; Initial value: 1.4e-5
    double var_intracellular_calcium_concentration__Ca_ds = rY[17];
    // Units: millimolar; Initial value: 1.88e-5

    const double var_membrane__R = 8314.472;
    const double var_membrane__T = 310.0;
    const double var_membrane__F = 96485.3415;
    double var_reversal_potentials__K_i = var_intracellular_potassium_concentration__K_i;
    double var_reversal_potentials__R = var_membrane__R;
    double var_reversal_potentials__T = var_membrane__T;
    double var_reversal_potentials__F = var_membrane__F;
    const double var_extracellular_potassium_concentration__K_o = 4.0;
    double var_reversal_potentials__K_o = var_extracellular_potassium_concentration__K_o;
    double var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i);
    double var_time_independent_potassium_current__E_K = var_reversal_potentials__E_K;
    double var_time_independent_potassium_current__K_o = var_extracellular_potassium_concentration__K_o;
    double var_time_independent_potassium_current__R = var_membrane__R;
    double var_time_independent_potassium_current__V = var_membrane__V;
    double var_time_independent_potassium_current__T = var_membrane__T;
    const double var_time_independent_potassium_current__K_mk1 = 10.0;
    const double var_time_independent_potassium_current__g_K1 = 0.5;
    double var_time_independent_potassium_current__F = var_membrane__F;
    double var_time_independent_potassium_current__i_K1 = (((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mk1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp((((var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K) - 10.0) * var_time_independent_potassium_current__F * 1.25) / (var_time_independent_potassium_current__R * var_time_independent_potassium_current__T)));
    double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
    double var_transient_outward_current__s = var_transient_outward_current_s_gate__s;
    double var_transient_outward_current__r = var_transient_outward_current_r_gate__r;
    const double var_transient_outward_current__g_to = 0.005;
    double var_transient_outward_current__V = var_membrane__V;
    double var_transient_outward_current__E_K = var_reversal_potentials__E_K;
    const double var_transient_outward_current__g_tos = 0.0;
    double var_transient_outward_current__i_to = var_transient_outward_current__g_to * (var_transient_outward_current__g_tos + (var_transient_outward_current__s * (1.0 - var_transient_outward_current__g_tos))) * var_transient_outward_current__r * (var_transient_outward_current__V - var_transient_outward_current__E_K);
    double var_membrane__i_to = var_transient_outward_current__i_to;
    const double var_rapid_delayed_rectifier_potassium_current__g_Kr2 = 0.0013;
    const double var_rapid_delayed_rectifier_potassium_current__g_Kr1 = 0.0021;
    double var_rapid_delayed_rectifier_potassium_current__xr1 = var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1;
    double var_rapid_delayed_rectifier_potassium_current__xr2 = var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2;
    double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
    double var_rapid_delayed_rectifier_potassium_current__E_K = var_reversal_potentials__E_K;
    double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((var_rapid_delayed_rectifier_potassium_current__g_Kr1 * var_rapid_delayed_rectifier_potassium_current__xr1) + (var_rapid_delayed_rectifier_potassium_current__g_Kr2 * var_rapid_delayed_rectifier_potassium_current__xr2)) * 1.0) / (1.0 + exp((var_rapid_delayed_rectifier_potassium_current__V + 9.0) / 22.4))) * (var_rapid_delayed_rectifier_potassium_current__V - var_rapid_delayed_rectifier_potassium_current__E_K);
    double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
    double var_slow_delayed_rectifier_potassium_current__xs = var_slow_delayed_rectifier_potassium_current_xs_gate__xs;
    const double var_extracellular_sodium_concentration__Na_o = 140.0;
    double var_reversal_potentials__Na_o = var_extracellular_sodium_concentration__Na_o;
    double var_reversal_potentials__Na_i = var_intracellular_sodium_concentration__Na_i;
    const double var_reversal_potentials__P_kna = 0.03;
    double var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i)));
    double var_slow_delayed_rectifier_potassium_current__E_Ks = var_reversal_potentials__E_Ks;
    const double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.0026;
    double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
    double var_slow_delayed_rectifier_potassium_current__i_Ks = var_slow_delayed_rectifier_potassium_current__g_Ks * pow(var_slow_delayed_rectifier_potassium_current__xs, 2.0) * (var_slow_delayed_rectifier_potassium_current__V - var_slow_delayed_rectifier_potassium_current__E_Ks);
    double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
    double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
    const double var_L_type_Ca_channel__FrICa = 1.0;
    double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
    double var_L_type_Ca_channel__K_o = var_extracellular_potassium_concentration__K_o;
    double var_L_type_Ca_channel__K_i = var_intracellular_potassium_concentration__K_i;
    double var_L_type_Ca_channel__F = var_membrane__F;
    const double var_L_type_Ca_channel__P_Ca_L = 0.1;
    double var_L_type_Ca_channel__T = var_membrane__T;
    const double var_L_type_Ca_channel__P_CaK = 0.002;
    double var_L_type_Ca_channel__V = var_membrane__V;
    double var_L_type_Ca_channel__f2 = var_L_type_Ca_channel_f2_gate__f2;
    double var_L_type_Ca_channel__R = var_membrane__R;
    double var_L_type_Ca_channel__i_Ca_L_K_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
    double var_L_type_Ca_channel__f2ds = var_L_type_Ca_channel_f2ds_gate__f2ds;
    double var_L_type_Ca_channel__i_Ca_L_K_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
    const double var_sodium_potassium_pump__i_NaK_max = 0.7;
    double var_sodium_potassium_pump__Na_i = var_intracellular_sodium_concentration__Na_i;
    double var_sodium_potassium_pump__K_o = var_extracellular_potassium_concentration__K_o;
    const double var_sodium_potassium_pump__K_mNa = 40.0;
    const double var_sodium_potassium_pump__K_mK = 1.0;
    double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_mK + var_sodium_potassium_pump__K_o)) * var_sodium_potassium_pump__Na_i) / (var_sodium_potassium_pump__K_mNa + var_sodium_potassium_pump__Na_i);
    double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
    const double var_fast_sodium_current__g_Na = 2.5;
    double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
    double var_fast_sodium_current__V = var_membrane__V;
    double var_reversal_potentials__E_mh = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__Na_o + (0.12 * var_reversal_potentials__K_o)) / (var_reversal_potentials__Na_i + (0.12 * var_reversal_potentials__K_i)));
    double var_fast_sodium_current__E_mh = var_reversal_potentials__E_mh;
    double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
    double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * (var_fast_sodium_current__V - var_fast_sodium_current__E_mh);
    double var_membrane__i_Na = var_fast_sodium_current__i_Na;
    double var_sodium_background_current__V = var_membrane__V;
    double var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i);
    double var_sodium_background_current__E_Na = var_reversal_potentials__E_Na;
    const double var_sodium_background_current__g_bna = 0.0006;
    double var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
    double var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
    const double var_persistent_sodium_current__g_pna = 0.004;
    double var_persistent_sodium_current__V = var_membrane__V;
    double var_persistent_sodium_current__E_Na = var_reversal_potentials__E_Na;
    double var_persistent_sodium_current__i_p_Na = ((var_persistent_sodium_current__g_pna * 1.0) / (1.0 + exp((-(var_persistent_sodium_current__V + 52.0)) / 8.0))) * (var_persistent_sodium_current__V - var_persistent_sodium_current__E_Na);
    double var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
    const double var_L_type_Ca_channel__P_CaNa = 0.01;
    double var_L_type_Ca_channel__Na_o = var_extracellular_sodium_concentration__Na_o;
    double var_L_type_Ca_channel__Na_i = var_intracellular_sodium_concentration__Na_i;
    double var_L_type_Ca_channel__i_Ca_L_Na_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
    double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
    double var_sodium_calcium_exchanger__Na_i = var_intracellular_sodium_concentration__Na_i;
    const double var_sodium_calcium_exchanger__n_NaCa = 3.0;
    const double var_sodium_calcium_exchanger__gamma = 0.5;
    double var_sodium_calcium_exchanger__F = var_membrane__F;
    double var_sodium_calcium_exchanger__Na_o = var_extracellular_sodium_concentration__Na_o;
    const double var_sodium_calcium_exchanger__FRiNaCa = 0.001;
    double var_sodium_calcium_exchanger__R = var_membrane__R;
    double var_sodium_calcium_exchanger__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_sodium_calcium_exchanger__T = var_membrane__T;
    double var_sodium_calcium_exchanger__V = var_membrane__V;
    const double var_sodium_calcium_exchanger__d_NaCa = 0.0;
    const double var_extracellular_calcium_concentration__Ca_o = 2.0;
    double var_sodium_calcium_exchanger__Ca_o = var_extracellular_calcium_concentration__Ca_o;
    const double var_sodium_calcium_exchanger__k_NaCa = 0.0005;
    double var_sodium_calcium_exchanger__i_NaCa_cyt = ((1.0 - var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_i))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_i * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_i / 0.0069)));
    double var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
    double var_sodium_calcium_exchanger__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
    double var_sodium_calcium_exchanger__i_NaCa_ds = (var_sodium_calcium_exchanger__FRiNaCa * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_ds))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_ds * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_ds / 0.0069)));
    double var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
    double var_L_type_Ca_channel__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_L_type_Ca_channel__Ca_o = var_extracellular_calcium_concentration__Ca_o;
    double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
    double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((var_L_type_Ca_channel__FrICa * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
    double var_reversal_potentials__Ca_o = var_extracellular_calcium_concentration__Ca_o;
    double var_reversal_potentials__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i);
    double var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca;
    const double var_calcium_background_current__g_bca = 0.00025;
    double var_calcium_background_current__V = var_membrane__V;
    double var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
    double var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;

    //////////////////////////////////////////////////////////////////////
    // new part of the model - addition of a Stretch-activated channel
    //////////////////////////////////////////////////////////////////////
    const double g_sac = 0.035; // uS
    const double E_sac = -10; // mV
    double f = (mStretch > 1.0) ? (mStretch-1.0)/0.15 : 0.0; // f = 0 if stretch < 1, scales linearly to f=1 at 15% stretch
    double sac_ionic_current = g_sac * f * (var_membrane__V - E_sac); // if g is uS, this is nA

    /*
     * The return value has to be scaled to match the units required by the mono/bidomain equations.
     * The cell model ionic current is in nano Amps, we require micro Amps/cm^2.
     * The estimate of the cell area is obtained by observing that Cm in the cell model and Cm in the bidomain equation are conceptually the same thing.
     * The Cm in the bidomain equation is expressed in capacitance units per area.
     * An estimate of the cell area is then the ratio of the two values of Cm.
     *
     */
    double value_in_nA = var_membrane__i_K1+var_membrane__i_to+var_membrane__i_Kr+var_membrane__i_Ks+var_membrane__i_Ca_L_K_cyt+var_membrane__i_Ca_L_K_ds+var_membrane__i_NaK+var_membrane__i_Na+var_membrane__i_b_Na+var_membrane__i_p_Na+var_membrane__i_Ca_L_Na_cyt+var_membrane__i_Ca_L_Na_ds+var_membrane__i_NaCa_cyt+var_membrane__i_NaCa_ds+var_membrane__i_Ca_L_Ca_cyt+var_membrane__i_Ca_L_Ca_ds+var_membrane__i_b_Ca;

//    std::cout << value_in_nA << " " << sac_ionic_current << "\n";
    value_in_nA += sac_ionic_current;

    double value_in_microA = 0.001*value_in_nA;
    double estimated_cell_surface_in_cm_square = 9.5e-05 / HeartConfig::Instance()->GetCapacitance();
    double value_in_microA_per_cm_square = value_in_microA/estimated_cell_surface_in_cm_square;
    return value_in_microA_per_cm_square;
}

void CML_noble_varghese_kohl_noble_1998_basic_with_sac::EvaluateYDerivatives (
        double var_environment__time,
        const std::vector<double> &rY,
        std::vector<double> &rDY)
{
    // Inputs:
    // Time units: second
    var_environment__time *= 0.001;
    double var_membrane__V = rY[0];
    // Units: millivolt; Initial value: -92.849333
    double var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = rY[1];
    // Units: dimensionless; Initial value: 1.03e-5
    double var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = rY[2];
    // Units: dimensionless; Initial value: 2e-7
    double var_slow_delayed_rectifier_potassium_current_xs_gate__xs = rY[3];
    // Units: dimensionless; Initial value: 0.001302
    double var_fast_sodium_current_m_gate__m = rY[4];
    // Units: dimensionless; Initial value: 0.0016203
    double var_fast_sodium_current_h_gate__h = rY[5];
    // Units: dimensionless; Initial value: 0.9944036
    double var_L_type_Ca_channel_d_gate__d = rY[6];
    // Units: dimensionless; Initial value: 0
    double var_L_type_Ca_channel_f_gate__f = rY[7];
    // Units: dimensionless; Initial value: 1
    double var_L_type_Ca_channel_f2_gate__f2 = rY[8];
    // Units: dimensionless; Initial value: 0.9349197
    double var_L_type_Ca_channel_f2ds_gate__f2ds = rY[9];
    // Units: dimensionless; Initial value: 0.9651958
    double var_transient_outward_current_s_gate__s = rY[10];
    // Units: dimensionless; Initial value: 0.9948645
    double var_transient_outward_current_r_gate__r = rY[11];
    // Units: dimensionless; Initial value: 0
    double var_calcium_release__ActFrac = rY[12];
    // Units: dimensionless; Initial value: 0.0042614
    double var_calcium_release__ProdFrac = rY[13];
    // Units: dimensionless; Initial value: 0.4068154
    double var_intracellular_sodium_concentration__Na_i = rY[14];
    // Units: millimolar; Initial value: 7.3321223
    double var_intracellular_potassium_concentration__K_i = rY[15];
    // Units: millimolar; Initial value: 136.5644281
    double var_intracellular_calcium_concentration__Ca_i = rY[16];
    // Units: millimolar; Initial value: 1.4e-5
    double var_intracellular_calcium_concentration__Ca_ds = rY[17];
    // Units: millimolar; Initial value: 1.88e-5
    double var_intracellular_calcium_concentration__Ca_up = rY[18];
    // Units: millimolar; Initial value: 0.4531889
    double var_intracellular_calcium_concentration__Ca_rel = rY[19];
    // Units: millimolar; Initial value: 0.4481927
    double var_intracellular_calcium_concentration__Ca_Calmod = rY[20];
    // Units: millimolar; Initial value: 0.0005555
    double var_intracellular_calcium_concentration__Ca_Trop = rY[21];
    // Units: millimolar; Initial value: 0.0003542


    // Mathematics
    const double var_membrane__R = 8314.472;
    const double var_membrane__T = 310.0;
    const double var_membrane__F = 96485.3415;
    const double var_membrane__Cm = 9.5e-05;
    double var_reversal_potentials__K_i = var_intracellular_potassium_concentration__K_i;
    double var_reversal_potentials__R = var_membrane__R;
    double var_reversal_potentials__T = var_membrane__T;
    double var_reversal_potentials__F = var_membrane__F;
    const double var_extracellular_potassium_concentration__K_o = 4.0;
    double var_reversal_potentials__K_o = var_extracellular_potassium_concentration__K_o;
    double var_reversal_potentials__E_K = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__K_o / var_reversal_potentials__K_i);
    double var_time_independent_potassium_current__E_K = var_reversal_potentials__E_K;
    double var_time_independent_potassium_current__K_o = var_extracellular_potassium_concentration__K_o;
    double var_time_independent_potassium_current__R = var_membrane__R;
    double var_time_independent_potassium_current__V = var_membrane__V;
    double var_time_independent_potassium_current__T = var_membrane__T;
    const double var_time_independent_potassium_current__K_mk1 = 10.0;
    const double var_time_independent_potassium_current__g_K1 = 0.5;
    double var_time_independent_potassium_current__F = var_membrane__F;
    double var_time_independent_potassium_current__i_K1 = (((var_time_independent_potassium_current__g_K1 * var_time_independent_potassium_current__K_o) / (var_time_independent_potassium_current__K_o + var_time_independent_potassium_current__K_mk1)) * (var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K)) / (1.0 + exp((((var_time_independent_potassium_current__V - var_time_independent_potassium_current__E_K) - 10.0) * var_time_independent_potassium_current__F * 1.25) / (var_time_independent_potassium_current__R * var_time_independent_potassium_current__T)));
    double var_membrane__i_K1 = var_time_independent_potassium_current__i_K1;
    double var_transient_outward_current__s = var_transient_outward_current_s_gate__s;
    double var_transient_outward_current__r = var_transient_outward_current_r_gate__r;
    const double var_transient_outward_current__g_to = 0.005;
    double var_transient_outward_current__V = var_membrane__V;
    double var_transient_outward_current__E_K = var_reversal_potentials__E_K;
    const double var_transient_outward_current__g_tos = 0.0;
    double var_transient_outward_current__i_to = var_transient_outward_current__g_to * (var_transient_outward_current__g_tos + (var_transient_outward_current__s * (1.0 - var_transient_outward_current__g_tos))) * var_transient_outward_current__r * (var_transient_outward_current__V - var_transient_outward_current__E_K);
    double var_membrane__i_to = var_transient_outward_current__i_to;
    const double var_rapid_delayed_rectifier_potassium_current__g_Kr2 = 0.0013;
    const double var_rapid_delayed_rectifier_potassium_current__g_Kr1 = 0.0021;
    double var_rapid_delayed_rectifier_potassium_current__xr1 = var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1;
    double var_rapid_delayed_rectifier_potassium_current__xr2 = var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2;
    double var_rapid_delayed_rectifier_potassium_current__V = var_membrane__V;
    double var_rapid_delayed_rectifier_potassium_current__E_K = var_reversal_potentials__E_K;
    double var_rapid_delayed_rectifier_potassium_current__i_Kr = ((((var_rapid_delayed_rectifier_potassium_current__g_Kr1 * var_rapid_delayed_rectifier_potassium_current__xr1) + (var_rapid_delayed_rectifier_potassium_current__g_Kr2 * var_rapid_delayed_rectifier_potassium_current__xr2)) * 1.0) / (1.0 + exp((var_rapid_delayed_rectifier_potassium_current__V + 9.0) / 22.4))) * (var_rapid_delayed_rectifier_potassium_current__V - var_rapid_delayed_rectifier_potassium_current__E_K);
    double var_membrane__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
    double var_slow_delayed_rectifier_potassium_current__xs = var_slow_delayed_rectifier_potassium_current_xs_gate__xs;
    const double var_extracellular_sodium_concentration__Na_o = 140.0;
    double var_reversal_potentials__Na_o = var_extracellular_sodium_concentration__Na_o;
    double var_reversal_potentials__Na_i = var_intracellular_sodium_concentration__Na_i;
    const double var_reversal_potentials__P_kna = 0.03;
    double var_reversal_potentials__E_Ks = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__K_o + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_o)) / (var_reversal_potentials__K_i + (var_reversal_potentials__P_kna * var_reversal_potentials__Na_i)));
    double var_slow_delayed_rectifier_potassium_current__E_Ks = var_reversal_potentials__E_Ks;
    const double var_slow_delayed_rectifier_potassium_current__g_Ks = 0.0026;
    double var_slow_delayed_rectifier_potassium_current__V = var_membrane__V;
    double var_slow_delayed_rectifier_potassium_current__i_Ks = var_slow_delayed_rectifier_potassium_current__g_Ks * pow(var_slow_delayed_rectifier_potassium_current__xs, 2.0) * (var_slow_delayed_rectifier_potassium_current__V - var_slow_delayed_rectifier_potassium_current__E_Ks);
    double var_membrane__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
    double var_L_type_Ca_channel__d = var_L_type_Ca_channel_d_gate__d;
    const double var_L_type_Ca_channel__FrICa = 1.0;
    double var_L_type_Ca_channel__f = var_L_type_Ca_channel_f_gate__f;
    double var_L_type_Ca_channel__K_o = var_extracellular_potassium_concentration__K_o;
    double var_L_type_Ca_channel__K_i = var_intracellular_potassium_concentration__K_i;
    double var_L_type_Ca_channel__F = var_membrane__F;
    const double var_L_type_Ca_channel__P_Ca_L = 0.1;
    double var_L_type_Ca_channel__T = var_membrane__T;
    const double var_L_type_Ca_channel__P_CaK = 0.002;
    double var_L_type_Ca_channel__V = var_membrane__V;
    double var_L_type_Ca_channel__f2 = var_L_type_Ca_channel_f2_gate__f2;
    double var_L_type_Ca_channel__R = var_membrane__R;
    double var_L_type_Ca_channel__i_Ca_L_K_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
    double var_L_type_Ca_channel__f2ds = var_L_type_Ca_channel_f2ds_gate__f2ds;
    double var_L_type_Ca_channel__i_Ca_L_K_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaK * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__K_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__K_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
    const double var_sodium_potassium_pump__i_NaK_max = 0.7;
    double var_sodium_potassium_pump__Na_i = var_intracellular_sodium_concentration__Na_i;
    double var_sodium_potassium_pump__K_o = var_extracellular_potassium_concentration__K_o;
    const double var_sodium_potassium_pump__K_mNa = 40.0;
    const double var_sodium_potassium_pump__K_mK = 1.0;
    double var_sodium_potassium_pump__i_NaK = (((var_sodium_potassium_pump__i_NaK_max * var_sodium_potassium_pump__K_o) / (var_sodium_potassium_pump__K_mK + var_sodium_potassium_pump__K_o)) * var_sodium_potassium_pump__Na_i) / (var_sodium_potassium_pump__K_mNa + var_sodium_potassium_pump__Na_i);
    double var_membrane__i_NaK = var_sodium_potassium_pump__i_NaK;
    const double var_fast_sodium_current__g_Na = 2.5;
    double var_fast_sodium_current__h = var_fast_sodium_current_h_gate__h;
    double var_fast_sodium_current__V = var_membrane__V;
    double var_reversal_potentials__E_mh = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log((var_reversal_potentials__Na_o + (0.12 * var_reversal_potentials__K_o)) / (var_reversal_potentials__Na_i + (0.12 * var_reversal_potentials__K_i)));
    double var_fast_sodium_current__E_mh = var_reversal_potentials__E_mh;
    double var_fast_sodium_current__m = var_fast_sodium_current_m_gate__m;
    double var_fast_sodium_current__i_Na = var_fast_sodium_current__g_Na * pow(var_fast_sodium_current__m, 3.0) * var_fast_sodium_current__h * (var_fast_sodium_current__V - var_fast_sodium_current__E_mh);
    double var_membrane__i_Na = var_fast_sodium_current__i_Na;
    double var_sodium_background_current__V = var_membrane__V;
    double var_reversal_potentials__E_Na = ((var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Na_o / var_reversal_potentials__Na_i);
    double var_sodium_background_current__E_Na = var_reversal_potentials__E_Na;
    const double var_sodium_background_current__g_bna = 0.0006;
    double var_sodium_background_current__i_b_Na = var_sodium_background_current__g_bna * (var_sodium_background_current__V - var_sodium_background_current__E_Na);
    double var_membrane__i_b_Na = var_sodium_background_current__i_b_Na;
    const double var_persistent_sodium_current__g_pna = 0.004;
    double var_persistent_sodium_current__V = var_membrane__V;
    double var_persistent_sodium_current__E_Na = var_reversal_potentials__E_Na;
    double var_persistent_sodium_current__i_p_Na = ((var_persistent_sodium_current__g_pna * 1.0) / (1.0 + exp((-(var_persistent_sodium_current__V + 52.0)) / 8.0))) * (var_persistent_sodium_current__V - var_persistent_sodium_current__E_Na);
    double var_membrane__i_p_Na = var_persistent_sodium_current__i_p_Na;
    const double var_L_type_Ca_channel__P_CaNa = 0.01;
    double var_L_type_Ca_channel__Na_o = var_extracellular_sodium_concentration__Na_o;
    double var_L_type_Ca_channel__Na_i = var_intracellular_sodium_concentration__Na_i;
    double var_L_type_Ca_channel__i_Ca_L_Na_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
    double var_L_type_Ca_channel__i_Ca_L_Na_ds = (((var_L_type_Ca_channel__FrICa * var_L_type_Ca_channel__P_CaNa * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Na_i * exp((50.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Na_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
    double var_sodium_calcium_exchanger__Na_i = var_intracellular_sodium_concentration__Na_i;
    const double var_sodium_calcium_exchanger__n_NaCa = 3.0;
    const double var_sodium_calcium_exchanger__gamma = 0.5;
    double var_sodium_calcium_exchanger__F = var_membrane__F;
    double var_sodium_calcium_exchanger__Na_o = var_extracellular_sodium_concentration__Na_o;
    const double var_sodium_calcium_exchanger__FRiNaCa = 0.001;
    double var_sodium_calcium_exchanger__R = var_membrane__R;
    double var_sodium_calcium_exchanger__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_sodium_calcium_exchanger__T = var_membrane__T;
    double var_sodium_calcium_exchanger__V = var_membrane__V;
    const double var_sodium_calcium_exchanger__d_NaCa = 0.0;
    const double var_extracellular_calcium_concentration__Ca_o = 2.0;
    double var_sodium_calcium_exchanger__Ca_o = var_extracellular_calcium_concentration__Ca_o;
    const double var_sodium_calcium_exchanger__k_NaCa = 0.0005;
    double var_sodium_calcium_exchanger__i_NaCa_cyt = ((1.0 - var_sodium_calcium_exchanger__FRiNaCa) * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_i))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_i * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_i / 0.0069)));
    double var_membrane__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
    double var_sodium_calcium_exchanger__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
    double var_sodium_calcium_exchanger__i_NaCa_ds = (var_sodium_calcium_exchanger__FRiNaCa * var_sodium_calcium_exchanger__k_NaCa * ((exp((var_sodium_calcium_exchanger__gamma * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_o) - (exp(((var_sodium_calcium_exchanger__gamma - 1.0) * (var_sodium_calcium_exchanger__n_NaCa - 2.0) * var_sodium_calcium_exchanger__V * var_sodium_calcium_exchanger__F) / (var_sodium_calcium_exchanger__R * var_sodium_calcium_exchanger__T)) * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa) * var_sodium_calcium_exchanger__Ca_ds))) / ((1.0 + (var_sodium_calcium_exchanger__d_NaCa * ((var_sodium_calcium_exchanger__Ca_ds * pow(var_sodium_calcium_exchanger__Na_o, var_sodium_calcium_exchanger__n_NaCa)) + (var_sodium_calcium_exchanger__Ca_o * pow(var_sodium_calcium_exchanger__Na_i, var_sodium_calcium_exchanger__n_NaCa))))) * (1.0 + (var_sodium_calcium_exchanger__Ca_ds / 0.0069)));
    double var_membrane__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
    double var_L_type_Ca_channel__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_L_type_Ca_channel__Ca_o = var_extracellular_calcium_concentration__Ca_o;
    double var_L_type_Ca_channel__i_Ca_L_Ca_cyt = ((((1.0 - var_L_type_Ca_channel__FrICa) * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2 * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
    double var_L_type_Ca_channel__i_Ca_L_Ca_ds = (((var_L_type_Ca_channel__FrICa * 4.0 * var_L_type_Ca_channel__P_Ca_L * var_L_type_Ca_channel__d * var_L_type_Ca_channel__f * var_L_type_Ca_channel__f2ds * (var_L_type_Ca_channel__V - 50.0) * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)) / (1.0 - exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T)))) * ((var_L_type_Ca_channel__Ca_i * exp((100.0 * var_L_type_Ca_channel__F) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))) - (var_L_type_Ca_channel__Ca_o * exp(((-(var_L_type_Ca_channel__V - 50.0)) * var_L_type_Ca_channel__F * 2.0) / (var_L_type_Ca_channel__R * var_L_type_Ca_channel__T))));
    double var_membrane__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
    double var_reversal_potentials__Ca_o = var_extracellular_calcium_concentration__Ca_o;
    double var_reversal_potentials__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_reversal_potentials__E_Ca = ((0.5 * var_reversal_potentials__R * var_reversal_potentials__T) / var_reversal_potentials__F) * log(var_reversal_potentials__Ca_o / var_reversal_potentials__Ca_i);
    double var_calcium_background_current__E_Ca = var_reversal_potentials__E_Ca;
    const double var_calcium_background_current__g_bca = 0.00025;
    double var_calcium_background_current__V = var_membrane__V;
    double var_calcium_background_current__i_b_Ca = var_calcium_background_current__g_bca * (var_calcium_background_current__V - var_calcium_background_current__E_Ca);
    double var_membrane__i_b_Ca = var_calcium_background_current__i_b_Ca;
    double var_membrane__i_Stim = GetStimulus((1.0/0.001)*var_environment__time);
    double var_rapid_delayed_rectifier_potassium_current_xr1_gate__V = var_rapid_delayed_rectifier_potassium_current__V;
    double var_rapid_delayed_rectifier_potassium_current_xr1_gate__alpha_xr1 = 50.0 / (1.0 + exp((-(var_rapid_delayed_rectifier_potassium_current_xr1_gate__V - 5.0)) / 9.0));
    double var_rapid_delayed_rectifier_potassium_current_xr1_gate__beta_xr1 = 0.05 * exp((-(var_rapid_delayed_rectifier_potassium_current_xr1_gate__V - 20.0)) / 15.0);
    double var_rapid_delayed_rectifier_potassium_current_xr2_gate__V = var_rapid_delayed_rectifier_potassium_current__V;
    double var_rapid_delayed_rectifier_potassium_current_xr2_gate__alpha_xr2 = 50.0 / (1.0 + exp((-(var_rapid_delayed_rectifier_potassium_current_xr2_gate__V - 5.0)) / 9.0));
    double var_rapid_delayed_rectifier_potassium_current_xr2_gate__beta_xr2 = 0.4 * exp(-pow((var_rapid_delayed_rectifier_potassium_current_xr2_gate__V + 30.0) / 30.0, 3.0));
    double var_slow_delayed_rectifier_potassium_current_xs_gate__V = var_slow_delayed_rectifier_potassium_current__V;
    double var_slow_delayed_rectifier_potassium_current_xs_gate__alpha_xs = 14.0 / (1.0 + exp((-(var_slow_delayed_rectifier_potassium_current_xs_gate__V - 40.0)) / 9.0));
    double var_slow_delayed_rectifier_potassium_current_xs_gate__beta_xs = 1.0 * exp((-var_slow_delayed_rectifier_potassium_current_xs_gate__V) / 45.0);
    double var_fast_sodium_current_m_gate__V = var_fast_sodium_current__V;
    double var_fast_sodium_current_m_gate__E0_m = var_fast_sodium_current_m_gate__V + 41.0;
    const double var_fast_sodium_current_m_gate__delta_m = 1e-05;
    double var_fast_sodium_current_m_gate__alpha_m = (fabs(var_fast_sodium_current_m_gate__E0_m) < var_fast_sodium_current_m_gate__delta_m) ? 2000.0 : ((200.0 * var_fast_sodium_current_m_gate__E0_m) / (1.0 - exp((-0.1) * var_fast_sodium_current_m_gate__E0_m)));
    double var_fast_sodium_current_m_gate__beta_m = 8000.0 * exp((-0.056) * (var_fast_sodium_current_m_gate__V + 66.0));
    double var_fast_sodium_current_h_gate__V = var_fast_sodium_current__V;
    const double var_fast_sodium_current_h_gate__shift_h = 0.0;
    double var_fast_sodium_current_h_gate__alpha_h = 20.0 * exp((-0.125) * ((var_fast_sodium_current_h_gate__V + 75.0) - var_fast_sodium_current_h_gate__shift_h));
    double var_fast_sodium_current_h_gate__beta_h = 2000.0 / (1.0 + (320.0 * exp((-0.1) * ((var_fast_sodium_current_h_gate__V + 75.0) - var_fast_sodium_current_h_gate__shift_h))));
    double var_L_type_Ca_channel__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
    const double var_L_type_Ca_channel__Km_f2 = 100000.0;
    const double var_L_type_Ca_channel__Km_f2ds = 0.001;
    const double var_L_type_Ca_channel__R_decay = 20.0;
    double var_L_type_Ca_channel_d_gate__V = var_L_type_Ca_channel__V;
    double var_L_type_Ca_channel_d_gate__E0_d = (var_L_type_Ca_channel_d_gate__V + 24.0) - 5.0;
    double var_L_type_Ca_channel_d_gate__alpha_d = (fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((30.0 * var_L_type_Ca_channel_d_gate__E0_d) / (1.0 - exp((-var_L_type_Ca_channel_d_gate__E0_d) / 4.0)));
    double var_L_type_Ca_channel_d_gate__beta_d = (fabs(var_L_type_Ca_channel_d_gate__E0_d) < 0.0001) ? 120.0 : ((12.0 * var_L_type_Ca_channel_d_gate__E0_d) / (exp(var_L_type_Ca_channel_d_gate__E0_d / 10.0) - 1.0));
    const double var_L_type_Ca_channel_d_gate__speed_d = 3.0;
    double var_L_type_Ca_channel_f_gate__V = var_L_type_Ca_channel__V;
    double var_L_type_Ca_channel_f_gate__E0_f = var_L_type_Ca_channel_f_gate__V + 34.0;
    const double var_L_type_Ca_channel_f_gate__delta_f = 0.0001;
    double var_L_type_Ca_channel_f_gate__alpha_f = (fabs(var_L_type_Ca_channel_f_gate__E0_f) < var_L_type_Ca_channel_f_gate__delta_f) ? 25.0 : ((6.25 * var_L_type_Ca_channel_f_gate__E0_f) / (exp(var_L_type_Ca_channel_f_gate__E0_f / 4.0) - 1.0));
    double var_L_type_Ca_channel_f_gate__beta_f = 12.0 / (1.0 + exp(((-1.0) * (var_L_type_Ca_channel_f_gate__V + 34.0)) / 4.0));
    const double var_L_type_Ca_channel_f_gate__speed_f = 0.3;
    double var_L_type_Ca_channel_f2_gate__Km_f2 = var_L_type_Ca_channel__Km_f2;
    double var_L_type_Ca_channel_f2_gate__Ca_i = var_L_type_Ca_channel__Ca_i;
    double var_L_type_Ca_channel_f2ds_gate__Km_f2ds = var_L_type_Ca_channel__Km_f2ds;
    double var_L_type_Ca_channel_f2ds_gate__R_decay = var_L_type_Ca_channel__R_decay;
    double var_L_type_Ca_channel_f2ds_gate__Ca_ds = var_L_type_Ca_channel__Ca_ds;
    double var_transient_outward_current_s_gate__V = var_transient_outward_current__V;
    double var_transient_outward_current_s_gate__alpha_s = 0.033 * exp((-var_transient_outward_current_s_gate__V) / 17.0);
    double var_transient_outward_current_s_gate__beta_s = 33.0 / (1.0 + exp((-0.125) * (var_transient_outward_current_s_gate__V + 10.0)));
    double var_transient_outward_current_r_gate__V = var_transient_outward_current__V;
    double var_sarcoplasmic_reticulum_calcium_pump__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_sarcoplasmic_reticulum_calcium_pump__Ca_up = var_intracellular_calcium_concentration__Ca_up;
    const double var_sarcoplasmic_reticulum_calcium_pump__alpha_up = 0.4;
    const double var_sarcoplasmic_reticulum_calcium_pump__beta_up = 0.03;
    const double var_sarcoplasmic_reticulum_calcium_pump__K_srca = 0.5;
    const double var_sarcoplasmic_reticulum_calcium_pump__K_xcs = 0.4;
    const double var_sarcoplasmic_reticulum_calcium_pump__K_cyca = 0.0003;
    double var_sarcoplasmic_reticulum_calcium_pump__K_1 = (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) / var_sarcoplasmic_reticulum_calcium_pump__K_srca;
    double var_sarcoplasmic_reticulum_calcium_pump__K_2 = var_sarcoplasmic_reticulum_calcium_pump__Ca_i + (var_sarcoplasmic_reticulum_calcium_pump__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_1) + (var_sarcoplasmic_reticulum_calcium_pump__K_cyca * var_sarcoplasmic_reticulum_calcium_pump__K_xcs) + var_sarcoplasmic_reticulum_calcium_pump__K_cyca;
    double var_sarcoplasmic_reticulum_calcium_pump__i_up = ((var_sarcoplasmic_reticulum_calcium_pump__Ca_i / var_sarcoplasmic_reticulum_calcium_pump__K_2) * var_sarcoplasmic_reticulum_calcium_pump__alpha_up) - (((var_sarcoplasmic_reticulum_calcium_pump__Ca_up * var_sarcoplasmic_reticulum_calcium_pump__K_1) / var_sarcoplasmic_reticulum_calcium_pump__K_2) * var_sarcoplasmic_reticulum_calcium_pump__beta_up);
    double var_calcium_translocation__Ca_rel = var_intracellular_calcium_concentration__Ca_rel;
    double var_calcium_translocation__Ca_up = var_intracellular_calcium_concentration__Ca_up;
    double var_calcium_translocation__i_trans = 50.0 * (var_calcium_translocation__Ca_up - var_calcium_translocation__Ca_rel);
    const double var_calcium_release__K_m_rel = 250.0;
    const double var_calcium_release__K_leak_rate = 0.05;
    double var_calcium_release__Ca_rel = var_intracellular_calcium_concentration__Ca_rel;
    double var_calcium_release__i_rel = ((pow(var_calcium_release__ActFrac / (var_calcium_release__ActFrac + 0.25), 2.0) * var_calcium_release__K_m_rel) + var_calcium_release__K_leak_rate) * var_calcium_release__Ca_rel;
    double var_calcium_release__V = var_membrane__V;
    double var_calcium_release__VoltDep = exp(0.08 * (var_calcium_release__V - 40.0));
    const double var_calcium_release__K_m_Ca_cyt = 0.0005;
    double var_calcium_release__Ca_i = var_intracellular_calcium_concentration__Ca_i;
    double var_calcium_release__CaiReg = var_calcium_release__Ca_i / (var_calcium_release__Ca_i + var_calcium_release__K_m_Ca_cyt);
    double var_calcium_release__Ca_ds = var_intracellular_calcium_concentration__Ca_ds;
    const double var_calcium_release__K_m_Ca_ds = 0.01;
    double var_calcium_release__CadsReg = var_calcium_release__Ca_ds / (var_calcium_release__Ca_ds + var_calcium_release__K_m_Ca_ds);
    double var_calcium_release__RegBindSite = var_calcium_release__CaiReg + ((1.0 - var_calcium_release__CaiReg) * var_calcium_release__CadsReg);
    double var_calcium_release__ActRate = (0.0 * var_calcium_release__VoltDep) + (500.0 * pow(var_calcium_release__RegBindSite, 2.0));
    double var_calcium_release__InactRate = 60.0 + (500.0 * pow(var_calcium_release__RegBindSite, 2.0));
    double var_calcium_release__PrecFrac = (1.0 - var_calcium_release__ActFrac) - var_calcium_release__ProdFrac;
    double var_calcium_release__SpeedRel = (var_calcium_release__V < (-50.0)) ? 5.0 : 1.0;
    const double var_intracellular_calcium_concentration__V_up_ratio = 0.01;
    const double var_intracellular_calcium_concentration__V_rel_ratio = 0.1;
    const double var_intracellular_calcium_concentration__V_e_ratio = 0.4;
    double var_intracellular_calcium_concentration__V_i_ratio = ((1.0 - var_intracellular_calcium_concentration__V_e_ratio) - var_intracellular_calcium_concentration__V_up_ratio) - var_intracellular_calcium_concentration__V_rel_ratio;
    const double var_intracellular_calcium_concentration__radius = 0.012;
    const double var_intracellular_calcium_concentration__length = 0.074;
    double var_intracellular_calcium_concentration__V_Cell = 3.141592654 * pow(var_intracellular_calcium_concentration__radius, 2.0) * var_intracellular_calcium_concentration__length;
    double var_intracellular_calcium_concentration__V_i = var_intracellular_calcium_concentration__V_Cell * var_intracellular_calcium_concentration__V_i_ratio;
    double var_intracellular_sodium_concentration__V_i = var_intracellular_calcium_concentration__V_i;
    double var_intracellular_sodium_concentration__F = var_membrane__F;
    double var_intracellular_sodium_concentration__i_Na = var_fast_sodium_current__i_Na;
    double var_intracellular_sodium_concentration__i_b_Na = var_sodium_background_current__i_b_Na;
    double var_intracellular_sodium_concentration__i_p_Na = var_persistent_sodium_current__i_p_Na;
    double var_intracellular_sodium_concentration__i_Ca_L_Na_cyt = var_L_type_Ca_channel__i_Ca_L_Na_cyt;
    double var_intracellular_sodium_concentration__i_Ca_L_Na_ds = var_L_type_Ca_channel__i_Ca_L_Na_ds;
    double var_intracellular_sodium_concentration__i_NaK = var_sodium_potassium_pump__i_NaK;
    double var_intracellular_sodium_concentration__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
    double var_intracellular_potassium_concentration__V_i = var_intracellular_calcium_concentration__V_i;
    double var_intracellular_potassium_concentration__i_K1 = var_time_independent_potassium_current__i_K1;
    double var_intracellular_potassium_concentration__i_Kr = var_rapid_delayed_rectifier_potassium_current__i_Kr;
    double var_intracellular_potassium_concentration__i_Ks = var_slow_delayed_rectifier_potassium_current__i_Ks;
    double var_intracellular_potassium_concentration__i_Ca_L_K_cyt = var_L_type_Ca_channel__i_Ca_L_K_cyt;
    double var_intracellular_potassium_concentration__i_Ca_L_K_ds = var_L_type_Ca_channel__i_Ca_L_K_ds;
    double var_intracellular_potassium_concentration__i_to = var_transient_outward_current__i_to;
    double var_intracellular_potassium_concentration__i_NaK = var_sodium_potassium_pump__i_NaK;
    double var_intracellular_potassium_concentration__F = var_membrane__F;
    const double var_intracellular_calcium_concentration__Calmod = 0.02;
    const double var_intracellular_calcium_concentration__Trop = 0.05;
    const double var_intracellular_calcium_concentration__alpha_Calmod = 100000.0;
    const double var_intracellular_calcium_concentration__beta_Calmod = 50.0;
    const double var_intracellular_calcium_concentration__alpha_Trop = 100000.0;
    const double var_intracellular_calcium_concentration__beta_Trop = 200.0;
    const double var_intracellular_calcium_concentration__V_ds_ratio = 0.1;
    const double var_intracellular_calcium_concentration__Kdecay = 10.0;
    double var_intracellular_calcium_concentration__i_up = var_sarcoplasmic_reticulum_calcium_pump__i_up;
    double var_intracellular_calcium_concentration__i_trans = var_calcium_translocation__i_trans;
    double var_intracellular_calcium_concentration__i_rel = var_calcium_release__i_rel;
    double var_intracellular_calcium_concentration__i_NaCa_cyt = var_sodium_calcium_exchanger__i_NaCa_cyt;
    double var_intracellular_calcium_concentration__i_NaCa_ds = var_sodium_calcium_exchanger__i_NaCa_ds;
    double var_intracellular_calcium_concentration__i_Ca_L_Ca_cyt = var_L_type_Ca_channel__i_Ca_L_Ca_cyt;
    double var_intracellular_calcium_concentration__i_Ca_L_Ca_ds = var_L_type_Ca_channel__i_Ca_L_Ca_ds;
    double var_intracellular_calcium_concentration__i_b_Ca = var_calcium_background_current__i_b_Ca;
    double var_intracellular_calcium_concentration__F = var_membrane__F;

    //////////////////////////////////////////////////////////////////////
    // new part of the model
    //////////////////////////////////////////////////////////////////////
    const double g_sac = 0.035; // uS
    const double E_sac = -10; // mV
    double f = (mStretch > 0) ? (mStretch-1.0)/0.15 : 0.0; // f = 0 if stretch < 1, scales linearly to f=1 at 15% stretch
    double sac_ionic_current = g_sac * f * (var_membrane__V - E_sac); // if g is uS, this is nA

    double d_dt_membrane__V;
    if (mSetVoltageDerivativeToZero)
    {
        d_dt_membrane__V = 0.0;
    }
    else
    {
        d_dt_membrane__V = ((-1.0) / var_membrane__Cm) * (sac_ionic_current + var_membrane__i_Stim + var_membrane__i_K1 + var_membrane__i_to + var_membrane__i_Kr + var_membrane__i_Ks + var_membrane__i_NaK + var_membrane__i_Na + var_membrane__i_b_Na + var_membrane__i_p_Na + var_membrane__i_Ca_L_Na_cyt + var_membrane__i_Ca_L_Na_ds + var_membrane__i_NaCa_cyt + var_membrane__i_NaCa_ds + var_membrane__i_Ca_L_Ca_cyt + var_membrane__i_Ca_L_Ca_ds + var_membrane__i_Ca_L_K_cyt + var_membrane__i_Ca_L_K_ds + var_membrane__i_b_Ca);
    }
    double d_dt_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1 = (var_rapid_delayed_rectifier_potassium_current_xr1_gate__alpha_xr1 * (1.0 - var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1)) - (var_rapid_delayed_rectifier_potassium_current_xr1_gate__beta_xr1 * var_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1);
    double d_dt_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2 = (var_rapid_delayed_rectifier_potassium_current_xr2_gate__alpha_xr2 * (1.0 - var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2)) - (var_rapid_delayed_rectifier_potassium_current_xr2_gate__beta_xr2 * var_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2);
    double d_dt_slow_delayed_rectifier_potassium_current_xs_gate__xs = (var_slow_delayed_rectifier_potassium_current_xs_gate__alpha_xs * (1.0 - var_slow_delayed_rectifier_potassium_current_xs_gate__xs)) - (var_slow_delayed_rectifier_potassium_current_xs_gate__beta_xs * var_slow_delayed_rectifier_potassium_current_xs_gate__xs);
    double d_dt_fast_sodium_current_m_gate__m = (var_fast_sodium_current_m_gate__alpha_m * (1.0 - var_fast_sodium_current_m_gate__m)) - (var_fast_sodium_current_m_gate__beta_m * var_fast_sodium_current_m_gate__m);
    double d_dt_fast_sodium_current_h_gate__h = (var_fast_sodium_current_h_gate__alpha_h * (1.0 - var_fast_sodium_current_h_gate__h)) - (var_fast_sodium_current_h_gate__beta_h * var_fast_sodium_current_h_gate__h);
    double d_dt_L_type_Ca_channel_d_gate__d = var_L_type_Ca_channel_d_gate__speed_d * ((var_L_type_Ca_channel_d_gate__alpha_d * (1.0 - var_L_type_Ca_channel_d_gate__d)) - (var_L_type_Ca_channel_d_gate__beta_d * var_L_type_Ca_channel_d_gate__d));
    double d_dt_L_type_Ca_channel_f_gate__f = var_L_type_Ca_channel_f_gate__speed_f * ((var_L_type_Ca_channel_f_gate__alpha_f * (1.0 - var_L_type_Ca_channel_f_gate__f)) - (var_L_type_Ca_channel_f_gate__beta_f * var_L_type_Ca_channel_f_gate__f));
    double d_dt_L_type_Ca_channel_f2_gate__f2 = 1.0 - (1.0 * ((var_L_type_Ca_channel_f2_gate__Ca_i / (var_L_type_Ca_channel_f2_gate__Km_f2 + var_L_type_Ca_channel_f2_gate__Ca_i)) + var_L_type_Ca_channel_f2_gate__f2));
    double d_dt_L_type_Ca_channel_f2ds_gate__f2ds = var_L_type_Ca_channel_f2ds_gate__R_decay * (1.0 - ((var_L_type_Ca_channel_f2ds_gate__Ca_ds / (var_L_type_Ca_channel_f2ds_gate__Km_f2ds + var_L_type_Ca_channel_f2ds_gate__Ca_ds)) + var_L_type_Ca_channel_f2ds_gate__f2ds));
    double d_dt_transient_outward_current_s_gate__s = (var_transient_outward_current_s_gate__alpha_s * (1.0 - var_transient_outward_current_s_gate__s)) - (var_transient_outward_current_s_gate__beta_s * var_transient_outward_current_s_gate__s);
    double d_dt_transient_outward_current_r_gate__r = 333.0 * ((1.0 / (1.0 + exp((-(var_transient_outward_current_r_gate__V + 4.0)) / 5.0))) - var_transient_outward_current_r_gate__r);
    double d_dt_calcium_release__ActFrac = (var_calcium_release__PrecFrac * var_calcium_release__SpeedRel * var_calcium_release__ActRate) - (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate);
    double d_dt_calcium_release__ProdFrac = (var_calcium_release__ActFrac * var_calcium_release__SpeedRel * var_calcium_release__InactRate) - (var_calcium_release__SpeedRel * 1.0 * var_calcium_release__ProdFrac);
    double d_dt_intracellular_sodium_concentration__Na_i = ((-1.0) / (1.0 * var_intracellular_sodium_concentration__V_i * var_intracellular_sodium_concentration__F)) * (var_intracellular_sodium_concentration__i_Na + var_intracellular_sodium_concentration__i_p_Na + var_intracellular_sodium_concentration__i_b_Na + (3.0 * var_intracellular_sodium_concentration__i_NaK) + (3.0 * var_intracellular_sodium_concentration__i_NaCa_cyt) + var_intracellular_sodium_concentration__i_Ca_L_Na_cyt + var_intracellular_sodium_concentration__i_Ca_L_Na_ds);
    double d_dt_intracellular_potassium_concentration__K_i = ((-1.0) / (1.0 * var_intracellular_potassium_concentration__V_i * var_intracellular_potassium_concentration__F)) * ((var_intracellular_potassium_concentration__i_K1 + var_intracellular_potassium_concentration__i_Kr + var_intracellular_potassium_concentration__i_Ks + var_intracellular_potassium_concentration__i_Ca_L_K_cyt + var_intracellular_potassium_concentration__i_Ca_L_K_ds + var_intracellular_potassium_concentration__i_to) - (2.0 * var_intracellular_potassium_concentration__i_NaK));
    double d_dt_intracellular_calcium_concentration__Ca_Trop = (var_intracellular_calcium_concentration__alpha_Trop * var_intracellular_calcium_concentration__Ca_i * (var_intracellular_calcium_concentration__Trop - var_intracellular_calcium_concentration__Ca_Trop)) - (var_intracellular_calcium_concentration__beta_Trop * var_intracellular_calcium_concentration__Ca_Trop);
    double d_dt_intracellular_calcium_concentration__Ca_Calmod = (var_intracellular_calcium_concentration__alpha_Calmod * var_intracellular_calcium_concentration__Ca_i * (var_intracellular_calcium_concentration__Calmod - var_intracellular_calcium_concentration__Ca_Calmod)) - (var_intracellular_calcium_concentration__beta_Calmod * var_intracellular_calcium_concentration__Ca_Calmod);
    double d_dt_intracellular_calcium_concentration__Ca_i = ((((((-1.0) / (2.0 * 1.0 * var_intracellular_calcium_concentration__V_i * var_intracellular_calcium_concentration__F)) * (((var_intracellular_calcium_concentration__i_Ca_L_Ca_cyt + var_intracellular_calcium_concentration__i_b_Ca) - (2.0 * var_intracellular_calcium_concentration__i_NaCa_cyt)) - (2.0 * var_intracellular_calcium_concentration__i_NaCa_ds))) + (var_intracellular_calcium_concentration__Ca_ds * var_intracellular_calcium_concentration__V_ds_ratio * var_intracellular_calcium_concentration__Kdecay) + ((var_intracellular_calcium_concentration__i_rel * var_intracellular_calcium_concentration__V_rel_ratio) / var_intracellular_calcium_concentration__V_i_ratio)) - d_dt_intracellular_calcium_concentration__Ca_Calmod) - d_dt_intracellular_calcium_concentration__Ca_Trop) - var_intracellular_calcium_concentration__i_up;
    double d_dt_intracellular_calcium_concentration__Ca_ds = (((-1.0) * var_intracellular_calcium_concentration__i_Ca_L_Ca_ds) / (2.0 * 1.0 * var_intracellular_calcium_concentration__V_ds_ratio * var_intracellular_calcium_concentration__V_i * var_intracellular_calcium_concentration__F)) - (var_intracellular_calcium_concentration__Ca_ds * var_intracellular_calcium_concentration__Kdecay);
    double d_dt_intracellular_calcium_concentration__Ca_up = ((var_intracellular_calcium_concentration__V_i_ratio / var_intracellular_calcium_concentration__V_up_ratio) * var_intracellular_calcium_concentration__i_up) - var_intracellular_calcium_concentration__i_trans;
    double d_dt_intracellular_calcium_concentration__Ca_rel = ((var_intracellular_calcium_concentration__V_up_ratio / var_intracellular_calcium_concentration__V_rel_ratio) * var_intracellular_calcium_concentration__i_trans) - var_intracellular_calcium_concentration__i_rel;

    rDY[0] = 0.001*d_dt_membrane__V;
    rDY[1] = 0.001*d_dt_rapid_delayed_rectifier_potassium_current_xr1_gate__xr1;
    rDY[2] = 0.001*d_dt_rapid_delayed_rectifier_potassium_current_xr2_gate__xr2;
    rDY[3] = 0.001*d_dt_slow_delayed_rectifier_potassium_current_xs_gate__xs;
    rDY[4] = 0.001*d_dt_fast_sodium_current_m_gate__m;
    rDY[5] = 0.001*d_dt_fast_sodium_current_h_gate__h;
    rDY[6] = 0.001*d_dt_L_type_Ca_channel_d_gate__d;
    rDY[7] = 0.001*d_dt_L_type_Ca_channel_f_gate__f;
    rDY[8] = 0.001*d_dt_L_type_Ca_channel_f2_gate__f2;
    rDY[9] = 0.001*d_dt_L_type_Ca_channel_f2ds_gate__f2ds;
    rDY[10] = 0.001*d_dt_transient_outward_current_s_gate__s;
    rDY[11] = 0.001*d_dt_transient_outward_current_r_gate__r;
    rDY[12] = 0.001*d_dt_calcium_release__ActFrac;
    rDY[13] = 0.001*d_dt_calcium_release__ProdFrac;
    rDY[14] = 0.001*d_dt_intracellular_sodium_concentration__Na_i;
    rDY[15] = 0.001*d_dt_intracellular_potassium_concentration__K_i;
    rDY[16] = 0.001*d_dt_intracellular_calcium_concentration__Ca_i;
    rDY[17] = 0.001*d_dt_intracellular_calcium_concentration__Ca_ds;
    rDY[18] = 0.001*d_dt_intracellular_calcium_concentration__Ca_up;
    rDY[19] = 0.001*d_dt_intracellular_calcium_concentration__Ca_rel;
    rDY[20] = 0.001*d_dt_intracellular_calcium_concentration__Ca_Calmod;
    rDY[21] = 0.001*d_dt_intracellular_calcium_concentration__Ca_Trop;
}


template<>
void OdeSystemInformation<CML_noble_varghese_kohl_noble_1998_basic_with_sac>::Initialise(void)
{
    // Time units: second
    //
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("millivolt");
    this->mInitialConditions.push_back(-92.849333);

    this->mVariableNames.push_back("xr1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.03e-5);

    this->mVariableNames.push_back("xr2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(2e-7);

    this->mVariableNames.push_back("xs");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.001302);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0016203);

    this->mVariableNames.push_back("h");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.9944036);

    this->mVariableNames.push_back("d");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("f");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1);

    this->mVariableNames.push_back("f2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.9349197);

    this->mVariableNames.push_back("f2ds");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.9651958);

    this->mVariableNames.push_back("s");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.9948645);

    this->mVariableNames.push_back("r");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0);

    this->mVariableNames.push_back("ActFrac");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0042614);

    this->mVariableNames.push_back("ProdFrac");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.4068154);

    this->mVariableNames.push_back("Na_i");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(7.3321223);

    this->mVariableNames.push_back("K_i");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(136.5644281);

    this->mVariableNames.push_back("Ca_i");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(1.4e-5);

    this->mVariableNames.push_back("Ca_ds");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(1.88e-5);

    this->mVariableNames.push_back("Ca_up");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.4531889);

    this->mVariableNames.push_back("Ca_rel");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.4481927);

    this->mVariableNames.push_back("Ca_Calmod");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.0005555);

    this->mVariableNames.push_back("Ca_Trop");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(0.0003542);

    this->mInitialised = true;
}

// Serialization for Boost>=1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CML_noble_varghese_kohl_noble_1998_basic_with_sac)
