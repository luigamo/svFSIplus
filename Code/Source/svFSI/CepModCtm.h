/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CEP_MOD_CTM_H 
#define CEP_MOD_CTM_H 

#include "Array.h"
#include "Vector.h"
#include <array>

#include <optional>
#include <functional>

//template <class T>
//T& make_ref(T&& x) { return x; }

/// @brief This module defines data structures for ten Tusscher-Panfilov
/// epicardial cellular activation model for cardiac electrophysiology
///
/// The classes defined here duplicate the data structures in the Fortran TPPMOD module defined 
/// in CEPMOD_TTP.f and PARAMS_TPP.f files. 
class CepModCtm
{
  public:
    CepModCtm();
    ~CepModCtm();


//--------------------------------------------------------------------
//
//     Constants for Courtemanche Myocyte Model.
//
//--------------------------------------------------------------------

//     Default model parameters
      /// Gas constant [J/mol/K]
      double Rc = 8.314472;

      /// Temperature [K]
      double Tc = 310.0;

      /// Faraday constant [C/mmol]
      double Fc = 96.4853415;

      /// Cell capacitance per unit surface area [uF/cm^{2}]
      double Cm = 100.0;

      /// Surface to volume ratio [um^{-1}]
      double sV = 0.2;

      /// Cytoplasmic volume [um^{3}]
      double V_c = 16.404E-3;

      /// Cellular resistivity [\f$\Omega\f$-cm]
      double rho = 162.0;

      /// Cytoplasmic volume [um^{3}]
      double Vcell = 20100.0;

      /// Intracelullar volume [um^{3}]
      double Vi = 13668.0;

      /// Sacroplasmic reticulum volume [um^{3}]
      double V_sr = 1.094E-3;

      /// Subspace volume [um^{3}]
      double V_ss = 5.468E-5;
      
      /// Vup: SR uptake compartment volume [um^{3}]
      double Vup = 1109.52;

      /// Vrel: SR uptake compartment volume [um^{3}]
      double Vrel = 96.48;

      /// Extracellular K concentration [mM]
      double K_o = 5.4;

      /// Extracellular Na concentration [mM]
      double Na_o = 140.0;

      /// Extracellular Ca concentration [mM]
      double Ca_o = 1.8;

      /// Maximal I_Na conductance [nS/pF]
      double G_Na = 7.8;

      /// Maximal I_K1 conductance [nS/pF]
      double G_K1 = 0.09;

      /// Maximal epicardial I_to conductance [nS/pF]
      Vector<double> G_to = {0.1652, 0.1652, 0.1652};

      /// Maximal I_Kr conductance [nS/pF]
      double G_Kr = 2.0*0.0294;

//     G_Kr for spiral wave breakup
//      double G_Kr = 0.172;     // units: nS/pF

      /// Maximal epicardial I_Ks conductance [nS/pF]
      Vector<double> G_Ks = {0.129, 0.129, 0.129};

//     G_Ks for spiral wave breakup (epi)
//      double G_Ks(3) = (/0.441, 0.392_RKIND, 0.098_RKIND/)

      /// Relative I_Ks permeability to Na [-]
      double p_KNa = 3.E-2;

      /// Maximal I_CaL conductance [cm^{3}/uF/ms]
      double G_CaL = 0.11142;

      /// Maximal Ib_Ca conductance [nS/pF]
      double Gb_Ca = 0.001131;

      /// Maximal Ib_Na conductance [nS/pF]
      double Gb_Na = 0.0006744375;

      /// Maximal IKur conductance [nS/pF]
      double G_Kur;

      /// Maximal I_NaCa [pA/pF]
      //double K_NaCa = 1000.;

      /// Maximal INaK current [pA/pF]
      double I_NaK_max = 0.6;

      /// Maximal INaCa current [pA/pF]
      double I_NaCa_max = 1600.0;

      /// Maximal IpCa current [pA/pF]
      double Ip_Ca_max = 0.275;

      /// Maximal Iup current [pA/pF]
      double I_up_max = 0.005;

      /// Q10- based temperatura adjustment factor
      double KQ10 = 3.0;

      /// Voltage dependent parameter of I_NaCa [-]
      double gamma = 0.35;

      /// Ca_i half-saturation constant for I_NaCa [mM]
      double Km_Ca = 1.38;

      /// Na_i half-saturation constant for I_NaCa [mM]
      double Km_Na = 87.5;

      /// Maximal release rate for Irel [1/ms]
      double krel = 30.0;

      /// Half-saturation constant of I_up [mM]
      double K_up = 0.00092;

      /// Saturation factor for I_NaCa [-]
      double K_sat = 0.1;

      /// Factor enhancing outward nature of I_NaCa [-]
      double alpha = 2.5;

      /// Maximal I_NaK [pA/pF]
      double p_NaK = 2.724;

      /// K_o half-saturation constant of I_NaK [mM]
      double Km_Ko = 1.5;

      /// Na_i half-saturation constant of I_NaK [mM]
      double Km_Nai = 10.0;

      /// Maximal I_pK conductance [nS/pF]
      double G_pK = 1.46E-2;

//     G_pK for spiral wave breakup
//      double G_pK = 2.19E-3;    // units: nS/pF

      /// Ca_i half-saturation constant for I_NaCa [mM]
      double K_mCa = 1.38;

      /// Na_i half-saturation constant for I_NaCa [mM]
      double K_mNai = 87.5;

      /// Maximal I_NaCa [pA/pF]
      double K_NaCa = 1000.;

      /// K_o half-saturation constant of I_NaK [mM]
      double K_mK = 1.;

      /// Na_i half-saturation constant of I_NaK [mM]
      double K_mNa = 40.;

      /// Maximal I_pCa conductance [pA/pF]
      double G_pCa = 0.1238;

//     G_pCa for spiral wave breakup
//      double G_pCa = 0.8666;    // units: pA/pF

      /// Half-saturation constant of I_pCa [mM]
      double K_pCa = 5.E-4;

      /// Maximal I_bNa conductance [nS/pF]
      double G_bNa = 2.9E-4;

      /// Maximal I_bCa conductance [nS/pF]
      double G_bCa = 5.92E-4;

      /// Maximal I_up conductance [mM/ms]
      double Vmax_up = 6.375E-3;

      /// Half-saturation constant of I_up [mM]
      //double K_up = 2.5E-4;

      /// Maximal I_rel conductance [mM/ms]
      double V_rel = 0.102;

      /// R to O and RI to I, I_rel transition rate [mM^{-2}/ms]
      double k1p = 0.15;

      /// O to I and R to RI, I_rel transition rate [mM^{-1}/ms]
      double k2p = 4.5E-2;

      /// O to R and I to RI, I_rel transition rate [ms^{-1}]
      double k3 = 6.E-2;

      /// I to O and Ri to I, I_rel transition rate [ms^{-1}]
      double k4 = 5.E-3;

      /// Ca_sr half-saturation constant of k_casr [mM]
      double EC = 1.5;

      /// Maximum value of k_casr [-]
      double max_sr = 2.5;

      /// Minimum value of k_casr [-]
      double min_sr = 1.;

      /// Maximal I_leak conductance [mM/ms]
      double V_leak = 3.6E-4;

      /// Maximal I_xfer conductance [mM/ms]
      double V_xfer = 3.8E-3;

      /// Total cytoplasmic buffer concentration [mM]
      double Buf_c = 0.2;

      /// Ca_i half-saturation constant for cytplasmic buffer [mM]
      double K_bufc = 1.E-3;

      /// Total sacroplasmic buffer concentration [mM]
      double Buf_sr = 10.;

      /// Ca_sr half-saturation constant for subspace buffer [mM]
      double K_bufsr = 0.3;

      /// Total subspace buffer concentration [mM]
      double Buf_ss = 0.4;

      /// Ca_ss half-saturation constant for subspace buffer [mM]
      double K_bufss = 2.5E-4;

      /// Resting potential [mV]
      double Vrest = -85.23;

      /// Maximal Ca concentration in uptake compartment [mM]
      double Caup_max = 15.0;

      /// Total calmodulin concentration in myoplasm [mM]
      double CMDN_max = 0.05;

      /// Total troponin concentration in myoplasm [mM]
      double TRPN_max = 0.07;

      /// Total calsequestrin concentration in myoplasm [mM]
      double CSQN_max = 10.0;

      /// Ca half-saturation constant for calmodulin [mM]
      double Km_Cmdn = 0.00238;

      /// Ca half-saturation constant for troponin [mM]
      double Km_Trpn = 0.0005;

      /// Ca half-saturation constant for calsequestrin [mM]
      double Km_Csqn = 0.8;

      /// Constants of I_KACh current [mM]
      double ACh = 0.005;                 
      double I_KACh_RA = 1.0;

      /// SACs constants
      // conductance: lambda = 1, g_sac = 0.0000049505 || lambda = 1.1, g_sac = 0.0000066597; || lambda = 1.2, g_sac = 0.0000089476;
      double lambda;
      double SAC_factor;
      double p_pNa = 1.0;
      double p_pK = 1.0;
      double p_pCa = 1.0;
      double p_ZNa = 1.0;
      double p_ZK = 1.0;
      double p_ZCa = 2.0;

      double g_sac;
      double Gsac = 0.0005; 
      double ksac = 100; 
      double alpha_sac = 3.0; 

//     Electromechanics coupling parameters: active stress model
      /// Resting Ca concentration [mM]
      double Ca_rest = 5.E-5;

      /// Critical Ca concentration [mM]
      double Ca_crit = 8.E-4;

      /// Saturation of concentration [MPa/mM]
      double eta_T = 12.5;

      /// Minimum activation [ms^{-1}]
      double eps_0 = 0.1;

      /// Maximum activation [ms^{-1}]
      double eps_i = 1. ;

      /// Transition rate [mM^{-1}]
      double xi_T = 4.E3;

//     Electromechanics coupling parameters: active strain model
//

      /// Active force of sacromere [-mM^{-2}]
      double alFa = -4.E6;

      /// Resting Ca concentration [mM]
      double c_Ca0 = 2.155E-4;

      /// Viscous-type constant [ms-mM^{-2}]
      double mu_Ca = 5.E6;

//     Force-length relationship parameters
      /// Initial length of sacromeres [um]
      double SL0 = 1.95;

      /// Min. length of sacromeres [um]
      double SLmin = 1.7;

      /// Max. length of sacromeres [um]
      double SLmax = 2.6;

      /// Fourier coefficients
      double f0  = -4333.618335582119;
      double fc1 =  2570.395355352195;
      double fs1 = -2051.827278991976;
      double fc2 =  1329.53611689133;
      double fs2 =  302.216784558222;
      double fc3 =  104.943770305116;
      double fs3 =  218.375174229422;

//-----------------------------------------------------------------------
//     Scaling factors
      /// Voltage scaling
      double Vscale  = 1.;

      /// Time scaling
      double Tscale  = 1.;

      /// Voltage offset parameter
      double Voffset = 0.;

//-----------------------------------------------------------------------
//     Variables
      /// Reverse potentials for Na, K, Ca
      double E_Na;
      double E_K;
      double E_Ca;
      double E_Ks;
//     Cellular transmembrane currents
      /// Fast sodium current
      double I_Na;

      /// inward rectifier outward current
      double I_K1;

      /// transient outward current
      double I_to;

      /// rapid delayed rectifier current
      double I_Kr;

      /// slow delayed rectifier current
      double I_Ks;

      /// L-type Ca current
      double I_CaL;

      /// Na-Ca exchanger current
      double I_NaCa;

      /// Na-K pump current
      double I_NaK;

      /// plateau Ca current
      double I_pCa;

      /// plateau K current
      double I_pK;

      /// background Ca current
      double I_bCa;

      /// background Na current
      double I_bNa;

      /// sacroplasmic reticulum Ca leak current
      double I_leak;

      /// sacroplasmic reticulum Ca pump current
      double I_up;

      /// Ca induced Ca release current
      double I_rel;

      /// diffusive Ca current
      double I_xfer;

            /// ultrarapid delayed rectifier K current 
      double I_Kur;

      /// inward rectifier K current
      double I_KACh;

      /// 
      double I_tr;

//-----------------------------------------------------------------------
//     State variables
      double V;
      double K_i;
      double Na_i;
      double Ca_i;
      double Ca_ss;
      double Ca_sr;
      double R_bar;
      double Kii;
      double Naii;
      double Caii;
      double Caupi;
      double Careli;

//     Gating variables (runtime, steady state)
      double xr1, xr1i;
      double xr2, xr2i;
      double f2, f2i;
      double fcass, fcassi;
      double s, si;
      double r, ri;
      double ui, u_infin, u_inf;
      double vii, v_inf;
      double wi, w_inf;
      double di, d_inf;
      double fCai, fCa_inf;
      double fi, f_inf;
      double hi, h_inf;
      double ji, j_inf;
      double mi, m_inf;
      double xri, xr_inf;
      double xsi, xs_inf;
      double oai, oa_inf;
      double uai, ua_inf;
      double uii, ui_inf;
      double oii, oi_inf;

//     Other variables
      double k1;
      double k2;
      double k_casr;
      double O;

//     Jacobian variables
      double E_Na_Nai, E_K_Ki, E_Ca_Cai, E_Ks_Ki, E_Ks_Nai;
      double I_Na_V, I_Na_Nai;
      double I_to_V, I_to_Ki;
      double I_K1_V, I_K1_Ki;
      double I_Kr_V, I_Kr_Ki;
      double I_Ks_V, I_Ks_Ki, I_Ks_Nai;
      double I_CaL_V, I_CaL_Cass;
      double I_NaCa_V, I_NaCa_Nai, I_NaCa_Cai;
      double I_NaK_V, I_NaK_Nai;
      double I_pCa_Cai;
      double I_pK_V, I_pK_Ki;
      double I_bCa_V, I_bCa_Cai;
      double I_bNa_V, I_bNa_Nai;
      double I_leak_Cai, I_leak_Casr;
      double I_up_Cai;
      double I_rel_Cass, I_rel_Casr, I_rel_Rbar;
      double I_xfer_Cai, I_xfer_Cass;
      double k_casr_sr, k1_casr, O_Casr, O_Cass, O_Rbar;

//    Land variables      
      //double I_TRPN = 0.0;
      double dlambda_dt = 0.0;

      double TRPN50 = 0.35;       /// Value of CaTRPN where B=0.5 in steady state 
      double TRPN_n = 2;          /// Cooperativity of Ca_TRPN binding rate 
      double k_trpn = 0.1;        /// Unbinding rate [/ms]
      double rs = 0.25;           /// Steady-state duty ratio 
      double rw = 0.5;            /// Steady-state ratio between U and W 

      double TOT_A = 25;          /// Effective distortion 
      double ktm_unblock = 1;     /// Tropomyosin rate constant [/ms]

      double beta_1 = -2.4;       /// Lenght-dipendent parameter 
      double beta_0 = 2.3;        /// Lenght-dipendent parameter
      double gamma_S = 0.0085;      /// S state rate 
      double gamma_wu = 0.615;    /// Distortion-dipendent unbinding rate of W state [/ms]
      double phi = 2.23; 

      double n_tm = 5;            /// Hill coefficient of cooperative activation 
      double ca50_ref = 0.86;     /// Calcium sensitivity adjusted for atria [μM]

      double k_ws = 3 * 0.012;    /// Crossbridge cycling rate adjusted for atria [μM]
      double k_uw = 0.182;        /// Crossbridge cycling rate
      double Tref = 120;          /// Maximal tension [kPa]
      
      /// Passive tension model parameters 
      double a = 2.1;             /// Parallel elastic element constant [kPa]
      double b = 9.1;             /// Parallel elastic element constant 
      double par_k = 7;           /// Series elastic element constant 
      double eta_l = 200;         /// Viscous dashpot element costant [/ms]
      double eta_s = 20;          /// Viscous dashpot element costant [/ms]

      double T;  
      double land_factor;
      double I_TRPN;
      double plotCaii;


    void actv_strn(const double c_Ca, const double I4f, const double dt, double& gf);
    void actv_strs(const double c_Ca, const double dt, double& Tact, double& epsX);

    void getf(const int i, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
        Vector<double>& dX, const double I_stim, const double K_sac, Vector<double>& RPAR);

    void getj(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, 
        Array<double>& JAC, const double Ksac);

    void init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg);

    void init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg,
        Vector<double>& X0, Vector<double>& Xg0);

    void integ_cn2(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg,
        const double Ts, const double dt, const double Istim, const double Ksac, 
        Vector<int>& IPAR, Vector<double>& RPAR);

    void integ_fe(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
        const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR);

    void integ_rk(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
        const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR);

    void update_g(const int i, const double dt, const int n, const int nG, const Vector<double>& X, 
        Vector<double>& Xg, Vector<double>& RPAR);

    void Land(Vector<double>& Y_land, Vector<double>& dY_land, double Caii, 
        double lambda, double dlambda_dt, double TRPN_max, double& T, double& I_TRPN);

};

#endif

