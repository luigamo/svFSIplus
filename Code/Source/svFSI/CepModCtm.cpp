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

#include "CepModCtm.h"

#include "mat_fun.h"
#include <math.h>
#include <iostream>
#include <fstream>

CepModCtm::CepModCtm()
{
}

CepModCtm::~CepModCtm()
{
}

/// @brief Compute macroscopic fiber strain based on sacromere force-length relationship and calcium concentration
void CepModCtm::actv_strn(const double c_Ca, const double I4f, const double dt, double& gf)
{
  // fiber length
  double SL = I4f * SL0;

  //  Sacromere force-length relationship
  if (SL >= SLmin && SL <= SLmax) {
    SL = 0.5*f0 + fc1*cos(SL) + fs1*sin(SL) + fc2*cos(2.0*SL) + fs2*sin(2.0*SL)  + fc3*cos(3.0*SL) + fs3*sin(3.0*SL);
  } else { 
    SL = 0.0;
  } 

  // Active force
  double Fa = alFa * (c_Ca-c_Ca0)*(c_Ca-c_Ca0) * SL;
  double rtmp = 2.0*I4f*(1.0/ pow(1.0+gf,3.0) - 1.0);
  gf = gf + dt*(Fa + rtmp)/(mu_Ca * c_Ca * c_Ca);
}

void CepModCtm::actv_strs(const double c_Ca, const double dt, double& Tact, double& epsX)
{
  epsX = exp(-exp(-xi_T*(c_Ca - Ca_crit)));
  epsX = eps_0 + (eps_i - eps_0)*epsX;
  double nr  = Tact + epsX*dt*eta_T*(c_Ca - Ca_rest);
  Tact = nr / (1.0 + epsX*dt);
}

/// @brief Compute currents and time derivatives of state variables
///
/// Note that is 'i' the myocardium zone id: 1, 2 or 3.
///
/// Reproduces Fortran 'GETF()'.
void CepModCtm::getf(const int i, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, Vector<double>& dX, 
    const double I_stim, const double K_sac, Vector<double>& RPAR)
{
  // Local copies of state variables
  double V      = X(0);
  double Kii    = X(1);
  double Naii   = X(2);
  double Caii   = X(3);
  double Caiupi = X(4);
  double Careli = X(5);
  double R_bar  = X(6);
  double XS     = X(7);
  double XW     = X(8);
  double Ca_TRPN = X(9);
  double TmBlocked = X(10);
  double ZETAS  = X(11);
  double ZETAW  = X(12);
  double Cd     = X(13);

  // Local copies of gating variables
  double ui   = Xg(0);
  double vii  = Xg(1);
  double wi   = Xg(2);
  double di   = Xg(3);
  double hi   = Xg(4);
  double xri  = Xg(5);
  double oii  = Xg(6);
  double uii  = Xg(7);
  double mi   = Xg(8);
  double ji   = Xg(9);
  double fi   = Xg(10);
  double xsi  = Xg(11);
  double oai  = Xg(12);
  double uai  = Xg(13);
  double fCai = Xg(14);
  double& T   = Xg(15);
  double& I_TRPN = Xg(16);

  // Stretch-activated currents
  //double I_sac = K_sac * (Vrest - V);


  double RT   = Rc * Tc / Fc;
  double E_K  = RT * log(K_o/Kii);
  double E_Na = RT * log(Na_o/Naii);
  double E_Ca = 0.5 * RT * log(Ca_o/Caii);
  double E_Ks = RT * log( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) );

  // I_Na: Fast sodium current
  double I_Na = G_Na * pow(mi,3.0) * hi * ji * (V - E_Na);

  // I_to: transient outward current
  double I_to = 1.0*(G_to[i-1] * (pow(oai,3.0)) * oii * (V - E_K));

  double I_K1 = 1.0 * (G_K1 * (V - E_K)/(1.0 + exp(0.07 * 
  (V + 80.0))));

  // I_Kr: rapid delayed rectifier current
  double I_Kr = 2.0 * (G_Kr * xri * (V - E_K)/(1.0 + 
  exp((V + 15.0)/22.4)));

  // I_Ks: slow delayed rectifier current
  double I_Ks = 1.0 * (G_Ks[i-1] * (pow(xsi,2.0)) * (V - E_K));

  // I_CaL: L-type Ca current
  double I_CaL = 0.90 * (G_CaL * di * fi * fCai * (V - 65.0));

  // I_NaCa: Na-Ca exchanger current

  double den = (((pow(Km_Na,3.0)) + (pow(Na_o,3.0))) * 
  (Km_Ca + Ca_o) * (1.0 + K_sat * 
  exp(((gamma - 1.0) * Fc * V)/(Rc * Tc))));
  double num1 = (exp((gamma * Fc * V)/(Rc * Tc)) * 
  (pow(Naii,3.0)) * Ca_o);
  double num2 = num1 - exp(((gamma - 1.0) * 
  Fc * V)/(Rc * Tc)) * (pow(Na_o,3.0)) * Caii; 
  double I_NaCa = I_NaCa_max * (num2/den);

  // I_NaK: Na-K pump current
  double sigma = (1.0/7.0) * (exp(Na_o/67.3) - 1.0);
  double fNaK = 1.0/(1.0 + 0.1245 * exp((-0.1 * Fc * V)/(Rc*Tc)) + 
  0.0365 * sigma * exp((-Fc * V)/(Rc * Tc)));
  double I_NaK = (K_o * I_NaK_max * fNaK/(1.0 + 
  (pow((Km_Nai/Naii),1.5))))/(K_o + Km_Ko);

  // I_pCa: plateau Ca current
  double I_pCa = Ip_Ca_max * Caii / (0.0005 + Caii);

  // I_bCa: background Ca current
  double I_bCa = Gb_Ca * (V - E_Ca);

  // I_bNa: background Na current
  double I_bNa = Gb_Na * (V - E_Na);

  // I_tr: transfer current from NSR to JSR
  double tau_tr = 180.0;
  double I_tr = (Caupi - Careli)/tau_tr;

  // I_Kur: ultrarapid delayed rectifier K current 
  double G_Kur = 0.005 + 0.05/(1.0 + exp((V - 15.0)/(-13.0)));
  double I_Kur = 1.0 * (G_Kur * (pow(uai,3.0))* uii * (V - E_K));

  // I_KACh: inward rectifier current
  double I_KACh = 1.0*(I_KACh_RA * (((0.08 + (0.04/(1.0 + 
  exp((V + 91.0)/12.0))))* (V - E_K))/(1.0 + (pow((0.03/ACh),2.1)))));

  // I_leak: Sacroplasmic Reticulum Ca leak current
  double I_leak = I_up_max * Caupi/Caup_max;

  // I_up: Sacroplasmic Reticulum Ca pump current
  double I_up = I_up_max/(1.0 + K_up/Caii);

  // I_rel: Ca induced Ca current (CICR)
  // double I_rel;
  if ((krel * (pow(ui,2.0) * vii * wi * (Careli - Caii))) < 1.0E-25) { 
            I_rel = 0.0;
  } else {
            I_rel = krel * (pow(ui,2.0)) * vii * wi * (Careli - Caii);
  }


  double g_sac = Gsac/(1 + ksac * exp(-(alpha_sac * (lambda - 1.0)))); 
   
  
  // Isac_Na: 
  double A_Na = p_pNa * g_sac * (pow(p_ZNa,2.0)) * (pow(Fc,2.0)) * 
  V/(Rc * Tc);
  double B_Na = (Naii - Na_o * exp(- p_ZNa * Fc * V/(Rc * Tc)));
  double C_Na = (1.0 - exp(- p_ZNa * Fc * V/(Rc * Tc)));
  double Isac_Na = SAC_factor*(A_Na*B_Na/C_Na);

  // Isac_K:
  double A_K = p_pK * g_sac * (pow(p_ZK,2.0)) * (pow(Fc,2.0)) * 
  V/(Rc * Tc);
  double B_K = (Kii - K_o * exp(- p_ZK * Fc * V/(Rc * Tc)));
  double C_K = (1.0 - exp(- p_ZK * Fc * V/(Rc * Tc)));
  double Isac_K = SAC_factor*(A_K * B_K/C_K);

  // Isac_Ca
  double A_Ca = p_pCa * g_sac * (pow(p_ZCa,2.0)) * (pow(Fc,2.0)) * 
  V/(Rc * Tc);
  double B_Ca = (Caii - Ca_o * exp(- p_ZCa * Fc * V/(Rc * Tc)));
  double C_Ca = (1.0 - exp(- p_ZCa * Fc * V/(Rc * Tc)));
  double Isac_Ca = SAC_factor*(A_Ca * B_Ca/C_Ca);

  // ISAC
  double Isac = Isac_Na + Isac_K + Isac_Ca;


  // Next variables will be used to update the intracelullar ion concentrations later
  // Ca buffers: 
  double Ca_Cmdn = (CMDN_max * Caii)/(Caii + Km_Cmdn);
  double Ca_Trpn = (TRPN_max * Caii)/(Caii + Km_Trpn);
  double Ca_Csqn = (CSQN_max * Careli)/(Careli + Km_Csqn);


  // Now compute time derivatives
  // dV/dt: rate of change of transmembrane voltage
  //
  dX(0)  = -(I_Na + I_to + I_K1 + I_Kr  + I_Ks + I_KACh + I_Kur + 
  I_CaL + I_NaK + I_pCa + I_bNa + I_bCa + I_NaCa + I_stim + Isac);

  // dKii/dt
  dX(1)  = (Cm/(Vi*Fc)) * (2.0 * I_NaK - (I_stim + I_K1 + I_to +
  I_Kr + I_Ks + I_KACh + I_Kur + Isac_K));

  //  dNaii/dt
  dX(2)  = (Cm/(Vi*Fc)) * (-3.0 * I_NaK - (3.0 * I_NaCa + I_Na + I_bNa + 
  Isac_Na));

  Land(X, dX, Caii,  lambda,  dlambda_dt, TRPN_max, T, I_TRPN);

  //double I_TRPN = dX(9) * TRPN_max;

  // dCaii/dt
  double B1 = Cm * (2.0 * I_NaCa - (I_CaL + I_pCa + I_bCa + 
  Isac_Ca))/(2.0 * Vi * Fc);
  double B11 = B1 + (Vup * (I_leak - I_up) + (I_rel*Vrel))/Vi - land_factor*I_TRPN;
  double B2 = 1.0 + (TRPN_max * Km_Trpn)/(pow((Caii + Km_Trpn),2.0)) + 
  (CMDN_max * Km_Cmdn)/(pow((Caii + Km_Cmdn),2.0));
  dX(3)  = B11/B2;

  // Plot Caii
  //plotCaii = (Caii + dX(3))*1000;

  // dCaupi/dt
  dX(4)  = (I_up - (I_leak + I_tr * Vrel/Vup));

  // dCareli/dt
  dX(5) = (I_tr - I_rel)/(1.0 + CSQN_max*Km_Csqn/(pow((Careli + Km_Csqn),2.0)));

  // Quantities to be written to file
  RPAR(2)  = I_Na;
  RPAR(3)  = I_K1;
  RPAR(4)  = I_to;
  RPAR(5)  = I_Kr;
  RPAR(6)  = I_Ks;
  RPAR(7)  = I_CaL;
  RPAR(8)  = I_NaCa;
  RPAR(9)  = I_NaK;
  RPAR(10) = I_pCa;
  RPAR(11) = I_bCa;
  RPAR(12) = I_bNa;
  RPAR(13) = I_tr;
  RPAR(14) = I_Kur;
  RPAR(15) = I_KACh;
  RPAR(16) = I_leak;
  RPAR(17) = I_up;
  RPAR(18) = I_rel;
}

void CepModCtm::getj(const int i, const int nX, const int nG, const Vector<double>& X, const Vector<double>& Xg, 
    Array<double>& JAC, const double Ksac)
{

  // Needs to be modify yet; it still has TTP equations 
  double RT, a, b, c, tau, sq5, e1, e2, e3, e4, n1, n2, d1, d2, d3;

  // Local copies of state variables
  double V     = X(0);
  double K_i   = X(1);
  double Na_i  = X(2);
  double Ca_i  = X(3);
  double Ca_ss = X(4);
  double Ca_sr = X(5);
  double R_bar = X(6);

  // Local copies of gating variables
  double xr1   = Xg(0);
  double xr2   = Xg(1);
  double xs    = Xg(2);
  double m     = Xg(3);
  double h     = Xg(4);
  double j     = Xg(5);
  double d     = Xg(6);
  double f     = Xg(7);
  double f2    = Xg(8);
  double fcass = Xg(9);
  double s     = Xg(10);
  double r     = Xg(11);

  RT = Rc * Tc / Fc;
  double E_K  = RT * log(K_o/K_i);
  double E_Na = RT * log(Na_o/Na_i);
  double E_Ca = 0.5 * RT * log(Ca_o/Ca_i);
  double E_Ks = RT * log( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) );

  E_K_Ki   = -RT / K_i;
  E_Na_Nai = -RT / Na_i;
  E_Ca_Cai = -RT / Ca_i / 2.0;
  E_Ks_Ki  = -RT / (K_i + p_KNa*Na_i);
  E_Ks_Nai = p_KNa * E_Ks_Ki;

  // I_Na: Fast sodium current
  I_Na = G_Na * pow(m,3.0) * h * j * (V - E_Na);
  I_Na_V   = G_Na * pow(m,3.0) * h * j;
  I_Na_Nai = I_Na_V * (-E_Na_Nai);

  // I_to: transient outward current
  I_to = G_to[i-1] * r * s * (V - E_K);
  I_to_V  = G_to[i-1] * r * s;
  I_to_Ki = I_to_V * (-E_K_Ki);

  // I_K1: inward rectifier outward current
  e1   = exp(0.060*(V - E_K - 200.0));
  e2   = exp(2.E-40*(V - E_K + 100.0));
  e3   = exp(0.10*(V - E_K - 10.0));
  e4   = exp(-0.50*(V - E_K));
  a    = 0.10/(1.0 + e1);
  b    = (3.0*e2 + e3) / (1.0 + e4);
  tau  = a / (a + b);
  sq5  = sqrt(K_o/5.40);
  n1   = -6.E-30*e1 / pow(1.0 + e1,2.0);
  n2   = (6.E-40*e2 + 0.10*e3 + 0.50*b*e4) / (1.0 + e4);
  n1   = (a + b)*n1 - a*(n1 + n2);
  d1   = pow(a + b,2.0);
  I_K1 = G_K1 * sq5 * tau * (V - E_K);
  I_K1_V  = G_K1 * sq5 * (tau + (V - E_K)*n1/d1);
  I_K1_Ki = I_K1_V * (-E_K_Ki);

  // I_Kr: rapid delayed rectifier current
  I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K);
  I_Kr_V   = G_Kr * sq5 * xr1 * xr2;
  I_Kr_Ki  = I_Kr_V * (-E_K_Ki);

  // I_Ks: slow delayed rectifier current
  I_Ks = G_Ks[i-1] * pow(xs,2.0) * (V - E_Ks);
  I_Ks_V   = G_Ks[i-1] * pow(xs,2.0);
  I_Ks_Ki  = I_Ks_V * (-E_Ks_Ki);
  I_Ks_Nai = I_Ks_V * (-E_Ks_Nai);

  //  I_CaL: L-type Ca current
  a = 2.0*(V-15.0)/RT;
  b = (0.250*Ca_ss*exp(a) - Ca_o) / (exp(a)-1.0);
  c = G_CaL * d * f * f2 * fcass * (2.0*a*Fc);
  n1 = (exp(a)/RT) / (exp(a)-1.0);
  I_CaL = c * b;
  I_CaL_V = I_CaL/(V-15.0) + n1*(c*0.50*Ca_ss - 2.0*I_CaL);
  I_CaL_Cass  = c * 0.250 * n1 * RT;

  //  I_NaCa: Na-Ca exchanger current
  e1 = exp(gamma*V/RT);
  e2 = exp((gamma-1.0)*V/RT);
  n1 = e1*pow(Na_i,3.0)*Ca_o - e2*pow(Na_o,3.0)*Ca_i*alpha;
  d1 = pow(K_mNai,3.0) + pow(Na_o,3.0);
  d2 = K_mCa + Ca_o;
  d3 = 1.0 + K_sat*e2;
  c  = 1.0/(d1*d2*d3);
  I_NaCa = K_NaCa * n1 * c;

  n1 = K_NaCa * c * ( e1*pow(Na_i,3.0)*Ca_o*(gamma/RT) - e2*pow(Na_o,3.0)*Ca_i*alpha*((gamma-1.0)/RT) );
  n2 = I_NaCa*K_sat*((gamma-1.0)/RT)*e2/d3;
  I_NaCa_V   = n1 - n2;
  I_NaCa_Nai =  K_NaCa * e1 * (3.0*pow(Na_i,2.0)) * Ca_o * c;
  I_NaCa_Cai = -K_NaCa * e2 * pow(Na_o,3.0) * alpha * c;

  // I_NaK: Na-K pump current
  e1 = exp(-0.10*V/RT);
  e2 = exp(-V/RT);
  n1 = p_NaK * K_o * Na_i;
  d1 = K_o + K_mK;
  d2 = Na_i + K_mNa;
  d3 = 1.0 + 0.12450*e1 + 0.03530*e2;
  I_NaK = n1 / (d1*d2*d3);
  n1 = (0.012450*e1 + 0.03530*e2)/RT;
  I_NaK_V = I_NaK * n1 / d3;
  I_NaK_Nai = I_NaK * K_mNa/(Na_i*d2);

  // I_pCa: plateau Ca current
  d1 = (K_pCa + Ca_i);
  I_pCa = G_pCa * Ca_i / d1;
  I_pCa_Cai = G_pCa * K_pCa/(d1*d1);

  //  I_pK: plateau K current
  e1 = exp((25.0-V)/5.980);
  I_pK  = G_pK * (V-E_K) / (1.0 + e1);
  I_pK_V  = (G_pK + I_pK*e1/5.980) / (1.0+e1);
  I_pK_Ki = G_pK * (-E_K_Ki) / (1.0 + e1);

  // I_bCa: background Ca current
  I_bCa = G_bCa * (V - E_Ca);
  I_bCa_V = G_bCa;
  I_bCa_Cai = G_bCa * (-E_Ca_Cai);

  //  I_bNa: background Na current
  I_bNa = G_bNa * (X(0) - E_Na);
  I_bNa_V = G_bNa;
  I_bNa_Nai = G_bNa * (-E_Na_Nai);

  // I_leak: Sacroplasmic Reticulum Ca leak current
  I_leak = V_leak * (Ca_sr - Ca_i);
  I_leak_Cai  = -V_leak;
  I_leak_Casr =  V_leak;

  // I_up: Sacroplasmic Reticulum Ca pump current
  d1 = 1.0 + pow(K_up/Ca_i,2.0);
  I_up  = Vmax_up / d1;
  I_up_Cai = (I_up / d1) * (2.0*pow(K_up,2.0) / pow(Ca_i,3.0));

  // I_rel: Ca induced Ca current (CICR)
  n1 = max_sr - min_sr;
  d1 = 1.0 + pow(EC/Ca_sr,2.0);
  k_casr = max_sr - (n1/d1);
  k1 = k1p / k_casr;
  n2 = Ca_ss*2.0;
  d2 = k3 + k1*n2;
  O = k1 * R_bar * n2 / d2;
  I_rel  = V_rel * O * (Ca_sr - Ca_ss);

  k_casr_sr = (n1 / pow(d1,2.0) ) * (2.0*pow(EC,2.0) / pow(Ca_sr,3.0));
  //k_casr_sr = (n1 / (d1**2.0)) * (2.0*EC**2.0 / Ca_sr**3.0);
  k1_casr   = -k1p * k_casr_sr / pow(k_casr,2.0);
  O_Cass = 2.0 * k3 * O / (Ca_ss * d2);
  O_Casr = k1_casr * n2 * (R_bar - O) / d2;
  O_Rbar = k1 * n2 / d2;

  I_rel_Cass = V_rel * (O_Cass*(Ca_sr - Ca_ss) - O);
  I_rel_Casr = V_rel * (O_Casr*(Ca_sr - Ca_ss) + O);
  I_rel_Rbar = V_rel * O_Rbar *(Ca_sr - Ca_ss);

  //  I_xfer: diffusive Ca current between Ca subspae and cytoplasm
  I_xfer = V_xfer * (Ca_ss - Ca_i);
  I_xfer_Cai  = -V_xfer;
  I_xfer_Cass =  V_xfer;

  // Compute Jacobian matrix
  //
  JAC = 0.0;
  c = Cm/(V_c*Fc);

  //  V
  JAC(0,0)  = -(I_Na_V + I_to_V + I_K1_V + I_Kr_V + I_Ks_V + I_CaL_V + I_NaCa_V + I_NaK_V + I_pK_V + I_bCa_V + I_bNa_V + Ksac);
  JAC(0,1)  = -(I_to_Ki + I_K1_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki);
  JAC(0,2)  = -(I_Na_Nai + I_Ks_Nai + I_NaCa_Nai + I_NaK_Nai + I_bNa_Nai);
  JAC(0,3)  = -(I_NaCa_Cai + I_pCa_Cai + I_bCa_Cai);
  JAC(0,4) = -I_CaL_Cass;

  // K_i
  JAC(1,0)  = -c * (I_K1_V + I_to_V + I_Kr_V + I_Ks_V + I_pK_V - 2.0*I_NaK_V );
  JAC(1,1)  = -c * (I_K1_Ki + I_to_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki);
  JAC(1,2)  = -c * (I_Ks_Nai - 2.0*I_NaK_Nai);

  // Na_i
  JAC(2,0)  = -c * (I_Na_V + I_bNa_V + 3.0*(I_NaK_V + I_NaCa_V));
  JAC(2,2)  = -c * (I_Na_Nai + I_bNa_Nai + 3.0*(I_NaK_Nai + I_NaCa_Nai));
  JAC(2,3)  = -c * (3.0*I_NaCa_Cai);

  //     Ca_i
  n1 = (I_leak - I_up)*V_sr/V_c + I_xfer - 0.50*c*(I_bCa + I_pCa - 2.0*I_NaCa);
  n2 = (I_leak_Cai - I_up_Cai)*V_sr/V_c + I_xfer_Cai - 0.50*c*(I_bCa_Cai + I_pCa_Cai - 2.0*I_NaCa_Cai);
  d1 = 1.0 + K_bufc*Buf_c / pow(Ca_i + K_bufc,2.0);
  d2 = 2.0*K_bufc*Buf_c / pow(Ca_i + K_bufc,3.0);
  JAC(3,0)  = -c * (I_bCa_V - 2.0*I_NaCa_V) / 2.0 / d1;
  JAC(3,2)  = c * I_NaCa_Nai / d1;
  JAC(3,3)  = (n2 + n1*d2/d1) / d1;
  JAC(3,4) = I_xfer_Cass / d1;
  JAC(3,5) = (I_leak_Casr*V_sr/V_c) / d1;

  //     Ca_ss
  a  = Cm/(2.0*Fc*V_ss);
  b  = V_sr / V_ss;
  c  = V_c / V_ss;
  n1 = -a*I_CaL + b*I_rel - c*I_xfer;
  n2 = -a*I_CaL_Cass + b*I_rel_Cass - c*I_xfer_Cass;
  d1 = 1.0 + K_bufss*Buf_ss/ pow(Ca_ss + K_bufss,2.0);
  d2 = 2.0*K_bufss*Buf_ss / pow(Ca_ss + K_bufss,3.0);
  JAC(4,0)  = -a * I_CaL_V / d1;
  JAC(4,3)  = -c * I_xfer_Cai / d1;
  JAC(4,4) = (n2 + n1*d2/d1) / d1;
  JAC(4,5) = b * I_rel_Casr / d1;
  JAC(4,6) = b * I_rel_Rbar / d1;

  // Ca_sr
  n1 = I_up - I_leak - I_rel;
  n2 = -(I_leak_Casr + I_rel_Casr);
  d1 = 1.0 + K_bufsr*Buf_sr / pow(Ca_sr + K_bufsr,2.0);
  d2 = 2.0*K_bufsr*Buf_sr / pow(Ca_sr + K_bufsr,3.0);
  JAC(5,3) = (I_up_Cai - I_leak_Cai) / d1;
  JAC(5,4) = -I_rel_Cass / d1;
  JAC(5,5) = (n2 + n1*d2/d1) / d1;
  JAC(5,6) = -I_rel_Rbar / d1;

  // Rbar: ryanodine receptor
  k2 = k2p * k_casr;
  JAC(6,4) = -k2 * R_bar;
  JAC(6,5) = -(k2p * k_casr_sr) * Ca_ss * R_bar;
  JAC(6,6) = -(k2*Ca_ss + k4);
}

void CepModCtm::init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg )
{
  switch (imyo) {

    // LA (left atrium)
    case 1:

      // // Initialize gating variables
      X(0)   = -81.2;           // V      (units: mV)
      X(1)   = 1.384571699E+02; // Kii    (units: mM)
      X(2)   = 11.81288054;     // Naii   (units: mM)
      X(3)   = 9.59608635E-05;  // Caii   (units: mM)
      X(4)   = 1.228880771;     // Caupi  (units: mM)
      X(5)   = 0.8208147476;    // Careli (units: mM)
      X(6)   =  0.9073;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  = 4.543092494E-53;   // ui    (dimensionless)
      Xg(1)  = 1.0;               // vii   (dimensionless)
      Xg(2)  = 0.9992293212;      // wi    (dimensionless)
      Xg(3)  = 1.257795372E-04;   // di    (dimensionless)
      Xg(4)  = 0.9696583637;      // hi    (dimensionless)
      Xg(5)  = 5.907522642E-04;   // xri   (dimensionless)
      Xg(6)  = 0.9993310315;      // oii   (dimensionless)
      Xg(7)  = 0.9910812856;      // uii   (dimensionless)
      Xg(8)  = 0.002606147731;    // mi    (dimensionless)
      Xg(9)  = 0.9807156912;      // ji    (dimensionless)
      Xg(10) = 0.9583644231;      // fi    (dimensionless)
      Xg(11) = 0.01851554132;     // xsi   (dimensionless)
      Xg(12) = 0.02933518307;     // oai   (dimensionless)
      Xg(13) = 0.004635522043;    // uai   (dimensionless)
      Xg(14) = 0.7847196193;      // fCai  (dimensionless) 

    break;

    // BBra (right Bachmann's bundle)
    case 2:

      // Initialize state variables
      X(0)   = -81.2;         // V      (units: mV)
      X(1)   =  137.42729671; // Kii    (units: mM)
      X(2)   =  12.48168081;  // Naii   (units: mM)
      X(3)   =  0.00013963;   // Caii   (units: mM)
      X(4)   =  1.88697642;   // Caupi  (units: mM)
      X(5)   =  1.42150433;   // Careli (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  =  1.58987004E-45;  // ui    (dimensionless)
      Xg(1)  =  0.99999999;      // vii   (dimensionless)
      Xg(2)  =  0.99919990;      // wi    (dimensionless)
      Xg(3)  =  0.00013621;      // di    (dimensionless)
      Xg(4)  =  0.96508770;      // hi    (dimensionless)
      Xg(5)  =  0.00316293;      // xri   (dimensionless)
      Xg(6)  =  0.99924435;      // oii   (dimensionless)
      Xg(7)  =  0.98688124;      // uii   (dimensionless)
      Xg(8)  =  0.00289452;      // mi    (dimensionless)
      Xg(9)  =  0.97746703;      // ji    (dimensionless)
      Xg(10) =  0.91112671;      // fi    (dimensionless)
      Xg(11) =  0.01959156;      // xsi   (dimensionless)
      Xg(12) =  0.03038954;      // oai   (dimensionless)
      Xg(13) =  0.00495269;      // uai   (dimensionless)
      Xg(14) =  0.71466706;      // fCai  (dimensionless)
    break;

    // TVR (right Tricuspid valve)
    case 3:

      // Initialize state variables
      X(0)   = -81.2;         // V      (units: mV)
      X(1)   = 138.36474905;  // Kii    (units: mM)
      X(2)   = 11.72135283;   // Naii   (units: mM)
      X(3)   = 8.95817415E-5; // Caii   (units: mM)
      X(4)   = 1.11033751;    // Caupi  (units: mM)
      X(5)   = 0.68712349;    // Careli (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  =  1.55290669E-50;   // ui    (dimensionless)
      Xg(1)  =  0.99999999;       // vii   (dimensionless)
      Xg(2)  =  0.99924076;      // wi    (dimensionless)
      Xg(3)  =  0.00012184;      // di    (dimensionless)
      Xg(4)  =  0.97130645;      // hi    (dimensionless)
      Xg(5)  =  0.00051798;      // xri   (dimensionless)
      Xg(6)  =  0.99936239;      // oii   (dimensionless)
      Xg(7)  =  0.99186450;      // uii   (dimensionless)
      Xg(8)  =  0.00249908;      // mi    (dimensionless)
      Xg(9)  =  0.98183379;      // ji    (dimensionless)
      Xg(10) =  0.96658221;      // fi    (dimensionless)
      Xg(11) =  0.01832886;      // xsi   (dimensionless)
      Xg(12) =  0.02892505;      // oai   (dimensionless)
      Xg(13) =  0.00451486;      // uai   (dimensionless)
      Xg(14) =  0.79610295;      // fCai  (dimensionless)
    break;

    // RAA (right atrial appendix)
    case 4:

      // Initialize state variables
      X(0)   = -81.2;           // V      (units: mV)
      X(1)   = 138.19164564;    // Kii    (units: mM)
      X(2)   = 11.85228357;     // Naii   (units: mM)
      X(3)   = 0.00010079;      // Caii   (units: mM)
      X(4)   = 1.29040614;      // Caupi  (units: mM)
      X(5)   = 0.84477390;      // Careli (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  =  3.75773862E-50;   // ui    (dimensionless)
      Xg(1)  =  0.99999999;       // vii   (dimensionless)
      Xg(2)  =  0.99923144;      // wi    (dimensionless)
      Xg(3)  =  0.00012504;      // di    (dimensionless)
      Xg(4)  =  0.96996417;      // hi    (dimensionless)
      Xg(5)  =  0.00106224;      // xri   (dimensionless)
      Xg(6)  =  0.99933671;      // oii   (dimensionless)
      Xg(7)  =  0.98970198;      // uii   (dimensionless)
      Xg(8)  =  0.00258605;      // mi    (dimensionless)
      Xg(9)  =  0.98091326;      // ji    (dimensionless)
      Xg(10) =  0.94911859;      // fi    (dimensionless)
      Xg(11) =  0.01860758;      // xsi   (dimensionless)
      Xg(12) =  0.02925946;      // oai   (dimensionless)
      Xg(13) =  0.00461313;      // uai   (dimensionless)
      Xg(14) =  0.77629516;      // fCai  (dimensionless) 
    break;

    // RA (right atrium)
    case 5:

      // Initialize state variables
      X(0)   = -81.2;      // V      (units: mV)                      
      X(1)   = 1.39E+02;  // Kii    (units: mM)                      
      X(2)   = 1.12E+01;  // Naii   (units: mM)                      
      X(3)   = 1.02E-04;  // Caii   (units: mM)                      
      X(4)   = 1.49;      // Caupi  (units: mM)                      
      X(5)   = 1.49;      // Careli (units: mM)  
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd                   
                           
      // Initialize gating variables
      Xg(0)  = 0.0;       // ui    (dimensionless)                       
      Xg(1)  = 1.0;       // vii   (dimensionless)                       
      Xg(2)  = 9.99E-01;  // wi    (dimensionless)                       
      Xg(3)  = 1.37E-04;  // di    (dimensionless)                       
      Xg(4)  = 9.65E-01;  // hi    (dimensionless)                       
      Xg(5)  = 3.29E-05;  // xri   (dimensionless)                       
      Xg(6)  = 9.99E-01;  // oii   (dimensionless)                       
      Xg(7)  = 9.99E-01;  // uii   (dimensionless)                       
      Xg(8)  = 2.91E-03;  // mi    (dimensionless)                       
      Xg(9)  = 9.78E-01;  // ji    (dimensionless)                       
      Xg(10) = 9.99E-01;  // fi    (dimensionless)                       
      Xg(11) = 1.87E-02;  // xsi   (dimensionless)                       
      Xg(12) = 3.04E-02;  // oai   (dimensionless)                       
      Xg(13) = 4.96E-03;  // uai   (dimensionless)                       
      Xg(14) = 7.75E-01;  // fCai  (dimensionless)  
    break;

    // PV (Pulmonary veins) 
    case 6:
      // Initialize state variables
      X(0)   = -81.2;         // V      (units: mV)
      X(1)   = 138.36252290;  // Kii    (units: mM)
      X(2)   = 11.71508629;   // Naii   (units: mM)
      X(3)   = 9.20684156E-5; // Caii   (units: mM)
      X(4)   = 1.14619006;    // Caupi  (units: mM)
      X(5)   = 0.71829680;    // Careli (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  =  1.83500055E-50;   // ui    (dimensionless)
      Xg(1)  =  0.99999999;       // vii   (dimensionless)
      Xg(2)  =  0.99920667;       // wi    (dimensionless)
      Xg(3)  =  0.00013377;       // di    (dimensionless)
      Xg(4)  =  0.96619626;       // hi    (dimensionless)
      Xg(5)  =  0.00064571;       // xri   (dimensionless)
      Xg(6)  =  0.99926566;       // oii   (dimensionless)
      Xg(7)  =  0.99134015;       // uii   (dimensionless)
      Xg(8)  =  0.00282656;       // mi    (dimensionless)
      Xg(9)  =  0.97830692;       // ji    (dimensionless)
      Xg(10) =  0.96323779;       // fi    (dimensionless)
      Xg(11) =  0.01891071;       // xsi   (dimensionless)
      Xg(12) =  0.03014627;       // oai   (dimensionless)
      Xg(13) =  0.00487862;       // uai   (dimensionless)
      Xg(14) =  0.79162114;       // fCai  (dimensionless) 
    break;

    // LAA (Left atrial appendix) 
    case 7:
     // Initialize state variables
      X(0)   = -81.2;         // V      (units: mV)
      X(1)   = 138.26628470;  // Kii    (units: mM)
      X(2)   = 11.79027177;   // Naii   (units: mM)
      X(3)   = 9.75452485E-5; // Caii   (units: mM)
      X(4)   = 1.23664245;    // Caupi  (units: mM)
      X(5)   = 0.79326261;    // Careli (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  =  2.08825591E-50;   // ui    (dimensionless)
      Xg(1)  =  0.99999999;       // vii   (dimensionless)
      Xg(2)  =  0.99922493;       // wi    (dimensionless)
      Xg(3)  =  0.00012731;       // di    (dimensionless)
      Xg(4)  =  0.96900345;       // hi    (dimensionless)
      Xg(5)  =  0.00080349;       // xri   (dimensionless)
      Xg(6)  =  0.99931851;       // oii   (dimensionless)
      Xg(7)  =  0.99098070;       // uii   (dimensionless)
      Xg(8)  =  0.00264801;       // mi    (dimensionless)
      Xg(9)  =  0.98025692;       // ji    (dimensionless)
      Xg(10) =  0.95331376;       // fi    (dimensionless)
      Xg(11) =  0.01868642;       // xsi   (dimensionless)
      Xg(12) =  0.02949291;       // oai   (dimensionless)
      Xg(13) =  0.00468233;       // uai   (dimensionless)
      Xg(14) =  0.78192307;       // fCai  (dimensionless) 
    break;

    // MV (Mitral valve) 
    case 8:
      // Initialize state variables
      X(0)   = -81.2;         // V      (units: mV)
      X(1)   = 138.35535119;  // Kii    (units: mM)
      X(2)   = 11.72897842;   // Naii   (units: mM)
      X(3)   = 9.75452485E-5; // Caii   (units: mM)
      X(4)   = 1.11778956;    // Caupi  (units: mM)
      X(5)   = 0.69356511;    // Careli (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  =  1.64253098E-50;   // ui    (dimensionless)
      Xg(1)  =  0.99999999;       // vii   (dimensionless)
      Xg(2)  =  0.99923193;       // wi    (dimensionless)
      Xg(3)  =  0.00012487;       // di    (dimensionless)
      Xg(4)  =  0.97004015;       // hi    (dimensionless)
      Xg(5)  =  0.00054710;       // xri   (dimensionless)
      Xg(6)  =  0.99933827;       // oii   (dimensionless)
      Xg(7)  =  0.99169150;       // uii   (dimensionless)
      Xg(8)  =  0.00258144;       // mi    (dimensionless)
      Xg(9)  =  0.98097530;       // ji    (dimensionless)
      Xg(10) =  0.96544845;       // fi    (dimensionless)
      Xg(11) =  0.01847208;       // xsi   (dimensionless)
      Xg(12) =  0.02924157;       // oai   (dimensionless)
      Xg(13) =  0.00460785;       // uai   (dimensionless)
      Xg(14) =  0.79510684;       // fCai  (dimensionless) 
    break;

    // BBla (left Bachmann's bundle) 
    case 9:
      // Initialize state variables
      X(0)   = -81.2;         // V      (units: mV)
      X(1)   = 137.39539547;  // Kii    (units: mM)
      X(2)   = 12.52015208;   // Naii   (units: mM)
      X(3)   = 9.75452485E-5; // Caii   (units: mM)
      X(4)   = 1.85733474;    // Caupi  (units: mM)
      X(5)   = 1.39071658;    // Careli (units: mM)
      X(6)   =  0.9068;     // R'     (dimensionless)
      X(7)   = 0.0;            // XS
      X(8)   = 0.0;            // XW
      X(9)   = 0.0;            // Ca_TRPN 
      X(10)   = 1.0;            // TmBlocked 
      X(11)   = 0.0;            // ZETAS 
      X(12)   = 0.0;            // ZETAW 
      X(13)   = 0.0;            // Cd

      // Initialize gating variables
      Xg(0)  =  1.25465339E-45;   // ui    (dimensionless)
      Xg(1)  =  0.99999999;       // vii   (dimensionless)
      Xg(2)  =  0.99919268;       // wi    (dimensionless)
      Xg(3)  =  0.00013884;       // di    (dimensionless)
      Xg(4)  =  0.96390313;       // hi    (dimensionless)
      Xg(5)  =  0.00307714;       // xri   (dimensionless)
      Xg(6)  =  0.99922226;       // oii   (dimensionless)
      Xg(7)  =  0.98628533;       // uii   (dimensionless)
      Xg(8)  =  0.00296827;       // mi    (dimensionless)
      Xg(9)  =  0.97662127;       // ji    (dimensionless)
      Xg(10) =  0.91592207;       // fi    (dimensionless)
      Xg(11) =  0.01961556;       // xsi   (dimensionless)
      Xg(12) =  0.03064754;       // oai   (dimensionless)
      Xg(13) =  0.00503184;       // uai   (dimensionless)
      Xg(14) =  0.71715014;       // fCai  (dimensionless) 
    break;
  }

}

void CepModCtm::init(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
    Vector<double>& X0, Vector<double>& Xg0)
{
  init(imyo, nX, nG, X, Xg);

  if (X0.size() != 0) {
    X = X0;
  }

  if (Xg0.size() != 0) {
    Xg = Xg0;
  }
}

/// @brief Time integration performed using Crank-Nicholson method
void CepModCtm::integ_cn2(const int imyo, const int nX, const int nG, Vector<double>& Xn, Vector<double>& Xg, 
    const double Ts, const double dt, const double Istim, const double Ksac, Vector<int>& IPAR, Vector<double>& RPAR)
{
  int itMax = IPAR(0);
  double atol  = RPAR(0);
  double rtol  = RPAR(1);
  auto Im = mat_fun::mat_id(nX);

  Vector<double> fn(nX);
  getf(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR);

  int k  = 0;
  auto Xk = Xn;
  bool l1 = false;
  bool l2 = false;
  bool l3 = false;
  double t = Ts + dt;
  double eps = std::numeric_limits<double>::epsilon();

  while (true) {
    k = k + 1;
    Vector<double> fk(nX);
    getf(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR);

    auto rK = Xk - Xn - 0.5*dt*(fk + fn);

    double rmsA = 0.0;
    double rmsR = 0.0;

    for (int i = 0; i < nX; i++) {
      rmsA = rmsA + pow(rK(i),2.0);
      rmsR = rmsR + pow(rK(i) / (Xk(i)+eps), 2.0);
    }

    rmsA = sqrt(rmsA / static_cast<double>(nX));
    rmsR = sqrt(rmsR / static_cast<double>(nX));

    l1 = (k > itMax);
    l2 = (rmsA <= atol);
    l3  = (rmsR <= rtol);
    if (l1 || l2 || l3) {
      break;
    }

    Array<double> JAC(nX,nX);
    getj(imyo, nX, nG, Xk, Xg, JAC, Ksac);

    JAC = Im - 0.5*dt*JAC;
    JAC = mat_fun::mat_inv(JAC, nX);
    rK  = mat_fun::mat_mul(JAC, rK);
    Xk  = Xk - rK;
  }

  Xn = Xk;

  update_g(imyo, dt, nX, nG, Xn, Xg, RPAR);
  getf(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR);

  if (!l2 && !l3) {
    IPAR(1) = IPAR(1) + 1;
  }
}

void CepModCtm::integ_fe(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
    const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR)
{
  Vector<double> f(nX);

  // Get time derivatives (RHS)
  getf(imyo, nX, nG, X, Xg, f, Istim, Ksac, RPAR);
  //CALL TTP_GETF(imyo, nX, nG, X, Xg, f, Istim, Ksac, RPAR)

  // Update gating variables
  update_g(imyo, dt, nX, nG, X, Xg, RPAR);
  //CALL TTP_UPDATEG(imyo, dt, nX, nG, X, Xg)

  //  Update state variables
  X = X + dt*f;
}

void CepModCtm::integ_rk(const int imyo, const int nX, const int nG, Vector<double>& X, Vector<double>& Xg, 
    const double Ts, const double dt, const double Istim, const double Ksac, Vector<double>& RPAR)
{
  double dt6 = dt / 6.0;

  Vector<double> frk1(nX), frk2(nX), frk3(nX), frk4(nX);

  // RK4: 1st pass
  auto Xrk = X;
  getf(imyo, nX, nG, Xrk, Xg, frk1, Istim, Ksac, RPAR);

  // Update gating variables by half-dt
  auto Xgr = Xg;
  update_g(imyo, 0.5*dt, nX, nG, X, Xgr, RPAR);

  // RK4: 2nd pass
  Xrk = X + 0.5*dt*frk1;
  getf(imyo, nX, nG, Xrk, Xgr, frk2, Istim, Ksac, RPAR);

  // RK4: 3rd pass
  Xrk = X + 0.5*dt*frk2;
  getf(imyo, nX, nG, Xrk, Xgr, frk3, Istim, Ksac, RPAR);

  // Update gating variables by full-dt
  Xgr = Xg;
  update_g(imyo, dt, nX, nG, X, Xgr, RPAR);

  // RK4: 4th pass
  Xrk = X + dt*frk3;
  getf(imyo, nX, nG, Xrk, Xgr, frk4, Istim, Ksac, RPAR);

  X = X + dt6 * (frk1 + 2.0*(frk2 + frk3) + frk4);
  Xg = Xgr;
}

/// @brief Update all the gating variables
void CepModCtm::update_g(const int i, const double dt, const int n, const int nG, const Vector<double>& X, Vector<double>& Xg, 
Vector<double>& RPAR)
{
  V  = X(0);
  Caii = X(3);

  I_NaCa = RPAR(8); 
  I_CaL  = RPAR(7);     
  I_rel  = RPAR(18);

  ui   =   Xg(0);   
  vii  =   Xg(1);       
  wi   =   Xg(2);       
  di   =   Xg(3);       
  hi   =   Xg(4);       
  xri  =   Xg(5);       
  oii  =   Xg(6);       
  uii  =   Xg(7);       
  mi   =   Xg(8);       
  ji   =   Xg(9);       
  fi   =   Xg(10);       
  xsi  =   Xg(11);       
  oai  =   Xg(12);       
  uai  =   Xg(13);       
  fCai =   Xg(14);  

  double Fn, tau_u, tau_v, tau_w, d_inf, tau_d, alfa_h, alfa_j, 
  beta_h, beta_j, tau_h, alfa_xr, beta_xr, tau_xr, alfa_oi, 
  beta_oi, tau_oi, alfa_ui, beta_ui, tau_ui, alfa_m, beta_m,
  tau_m, tau_j, tau_f, alfa_xs, beta_xs, tau_xs, alfa_oa,
  beta_oa, tau_oa, alfa_ua, beta_ua, tau_ua, tau_fCa;

  // Fn 
  Fn = 1.0E+03 * (1.0E-15 * Vrel * I_rel - Cm * 
  ((1.0E-15/(2.0 * Fc)) * (0.5 * I_CaL - 0.2 * I_NaCa)));
    
  // ui: Ca release current from JSR u gate
  tau_u = 8.0;
  u_infin = 1.0/(1.0 + exp(-(Fn - 3.4175E-13 - 13.0)/13.67E-16));
  if ((1.0/(1.0 + exp(-(Fn - 3.4175E-13)/13.67E-16)))< 1.0E-25) {
            u_inf = 0.0;
  } else { 
            u_inf = 1.0/(1.0 + exp(-(Fn - 3.4175E-13)/13.67E-16));
  }
  // Xg(0) = ui + dt * ((u_inf - ui)/tau_u);
  Xg(0) = u_inf - (u_inf - ui)*exp(-dt/tau_u);


  // vii: Ca release current from JSR v gate
  tau_v = 1.91 + 2.09/(1.0 + exp(-(Fn - 3.4175E-13)/13.67E-16));
  v_inf = 1.0 - 1.0/(1.0 + exp(-(Fn - 6.835E-14)/13.67E-16));
  // Xg(1) = vii + dt*((v_inf - vii)/tau_v);
  Xg(1) = v_inf - (v_inf - vii)*exp(-dt/tau_v); 


  // wi: Ca release current from JSR w gate
  tau_w;
  if (abs(V - 7.9)< 1.0E-10){
        tau_w = 6.0 * 0.2/1.3;
  } else {
        tau_w = 6.0 * (1.0 - exp(-(V - 7.9)/5.0))/((1.0 + 
        0.3 * exp(-(V - 7.9)/5.0)) * (V - 7.9));
  }
        w_inf = 1.0 - 1.0/(1.0 + exp(-(V - 40.0)/17.0));
 // Xg(2) = wi + dt*((w_inf - wi)/tau_w);
  Xg(2) = w_inf - (w_inf - wi)*exp(-dt/tau_w);


  // di: L type Ca channel d gate
  d_inf = 1.0/(1.0 + exp((V + 10.0)/(-8.0)));
  if (abs(V + 10.0)<1.0E-10){
        tau_d = 4.579/(1.0 + exp((V + 10.0)/(-6.24)));
  } else {
        tau_d = (1.0 - exp((V + 10.0)/(-6.24)))/
        (0.035 * (V + 10.0) * (1.0 + exp((V + 10.0)/(-6.24))));
  }
  // Xg(3) = di + dt * ((d_inf - di)/tau_d);
  Xg(3) = d_inf - (d_inf - di)*exp(-dt/tau_d);

  // hi: fast inactivation gate for I_Na
  if (V < -40.0){
      alfa_h = 0.135 * exp((V + 80.0)/(-6.8));
      beta_h = 3.56 * exp(0.079 * V) + 3.1E+05 * exp(0.35 * V);
      alfa_j = (-1.2714E+05 * exp(0.2444 * V) - 3.474E-05 * 
      exp(-0.04391 * V)) * (V + 37.78)/(1.0 + exp(0.311 * (V + 79.23)));
      beta_j = (0.1212 * exp(-0.01052 * V))/(1.0 + exp(-0.1378 * (V + 40.14)));
  } else {
      alfa_h = 0.0;
      beta_h = 1.0/(0.13 * (1.0 + exp((V + 10.66)/(-11.1))));
      alfa_j = 0.0;
      beta_j = (0.3 * exp(-2.535E-07 * V))/(1.0 + exp(-0.1 * (V + 32.0)));
  }
  h_inf = alfa_h/(alfa_h + beta_h);
  tau_h = 1.0/(alfa_h + beta_h);
  if ((hi + dt * ((h_inf - hi)/tau_h))< 1.0E-25){
      Xg(4)= 0.0;
  } else {
      // Xg(4) = hi + dt * ((h_inf - hi)/tau_h);
      Xg(4) = h_inf - (h_inf - hi)*exp(-dt/tau_h);
  }

  // xri: Rapid delayed rectifier K current xr gate
  if (abs(V + 14.1)< 1.0E-10){
      alfa_xr = 0.0015;
  } else {
      alfa_xr = 0.0003 * (V + 14.1)/(1.0 - exp((V + 14.1)/(-5.0)));
  }
  if (abs(V - 3.3328)< 1.0E-10){
      beta_xr = 3.7836118E-04;
  } else {
      beta_xr = 0.000073898 * (V - 3.3328)/(exp((V - 3.3328)/(5.1237))- 1.0);
  }
  tau_xr = 1.0/(alfa_xr + beta_xr);
  xr_inf = 1.0/(1.0 + exp((V + 14.1)/(-6.5)));
  //Xg(5) = xri + dt * ((xr_inf - xri)/tau_xr);
  Xg(5) = xr_inf - (xr_inf - xri)*exp(-dt/tau_xr);
  
  // oii: Transient outward K current oi gate
  alfa_oi = 1.0/(18.53 + exp((V + 113.7)/10.95));
  beta_oi = 1.0/(35.56 + exp((V + 1.26)/(-7.44)));
  tau_oi = 1.0/((alfa_oi + beta_oi) * KQ10);
  oi_inf = 1.0/(1.0 + exp((V + 43.1)/5.3));
  // Xg(6) = oii + dt * ((oi_inf - oii)/tau_oi);
  Xg(6) = oi_inf - (oi_inf - oii)*exp(-dt/tau_oi);

  // uii: Ultrarapid delayed rectifier K current ui gate 
  alfa_ui = 1.0/(21.0 + exp((V - 185.0)/(-28.0)));
  beta_ui = 1.0/(exp((V - 158.0)/(-16.0)));
  tau_ui = 1.0/((alfa_ui + beta_ui) * KQ10);
  ui_inf = 1.0/(1.0 + exp((V - 99.45)/(27.48)));
  // Xg(7) = uii + dt * ((ui_inf - uii)/tau_ui);
  Xg(7) = ui_inf - (ui_inf - uii)*exp(-dt/tau_ui);


  // mi: Fast sodium current m gate
  if (V == -47.13){
      alfa_m = 3.2;
  } else { 
      alfa_m = 0.32 * (V + 47.13)/(1.0 - exp(-0.1 * (V + 47.13)));
  }
  beta_m = 0.08 * exp(-V/11.0);
  m_inf = alfa_m/(alfa_m + beta_m);
  tau_m = 1.0/(alfa_m + beta_m);
  //Xg(8) = mi + dt * ((m_inf - mi)/tau_m);
  Xg(8) = m_inf - (m_inf - mi)*exp(-dt/tau_m); 


  // ji: Fast sodium current ji 
  j_inf = alfa_j/(alfa_j + beta_j);
  tau_j = 1.0/(alfa_j + beta_j);
  //Xg(9) = ji + dt * ((j_inf - ji)/tau_j);
  Xg(9) = j_inf - (j_inf - ji)*exp(-dt/tau_j);

  // fi: L type Ca channel f gate
  f_inf = exp(-(V + 28.0)/6.9)/(1.0 + exp(-(V + 28.0)/6.9));
  tau_f = 9.0/(0.0197 * exp(-(pow(0.0337,2.0)) * pow((V + 10.0),2.0)) + 0.02);
  // Xg(10) = fi + dt*((f_inf - fi)/tau_f); 
  Xg(10) = f_inf - (f_inf - fi)*exp(-dt/tau_f);

  // xsi: Slow delayed rectifier K current xs gate
  if (abs(V - 19.9)<1.0E-10){
      alfa_xs = 0.00068;
      beta_xs = 0.000315;
  } else {
      alfa_xs = 0.00004 * (V - 19.9)/(1.0 - exp((V - 19.9)/(-17.0)));
      beta_xs = 0.000035 * (V - 19.9)/(exp((V - 19.9)/9.0) - 1.0);
  }
  tau_xs = 0.5/((alfa_xs + beta_xs));
  xs_inf = 1.0/(sqrt(1.0 + exp((V - 19.9)/(-12.7))));
  // Xg(11) = xsi + dt * ((xs_inf - xsi)/tau_xs);
  Xg(11) = xs_inf - (xs_inf - xsi)*exp(-dt/tau_xs);

  // oai: Transient outward K current oa gate
  alfa_oa = 0.65/(exp((V + 10.0)/(-8.5)) + exp((V - 30.0)/(-59.0)));
  beta_oa = 0.65/(2.5 + exp((V + 82.0)/17.0));
  tau_oa = 1.0/((alfa_oa + beta_oa) * KQ10);
  oa_inf = 1.0/(1.0 + exp((V + 20.47)/(-17.54)));
  // Xg(12) = oai + dt * ((oa_inf - oai)/tau_oa);
  Xg(12) = oa_inf - (oa_inf - oai)*exp(-dt/tau_oa);

  // uai: Ultrarapid delayed rectifier K current ua gate
  alfa_ua = 0.65/(exp((V + 10.0)/(-8.5)) + exp((V - 30.0)/(-59.0)));
  beta_ua = 0.65/(2.5 + exp((V + 82.0)/17.0));
  tau_ua = 1.0/((alfa_ua + beta_ua) * KQ10);
  ua_inf = 1.0/(1.0 + exp((V + 30.3)/(-9.6)));
  // Xg(13) = uai + dt * ((ua_inf - uai)/tau_ua);
  Xg(13) = ua_inf - (ua_inf - uai)*exp(-dt/tau_ua);

  // fCai: L type Ca channel fCa gate
  fCa_inf = 1.0/(1.0 + Caii/0.00035);
  tau_fCa = 2.0;
  // Xg(14) = fCai + dt*((fCa_inf - fCai)/tau_fCa);
  Xg(14) = fCa_inf - (fCa_inf - fCai)*exp(-dt/tau_fCa); 
}
 
/// @brief
void CepModCtm::Land(Vector<double>& Y_land, Vector<double>& dY_land, double Caii, double lambda, double dlambda_dt, double TRPN_max, double& T, double& I_TRPN) {
    
    // State Variables
    double XS = std::max(0.0, Y_land(7));
    double XW = std::max(0.0, Y_land(8));
    double Ca_TRPN = std::max(0.0, Y_land(9));
    double TmBlocked = Y_land(10);
    double ZETAS = Y_land(11);
    double ZETAW = Y_land(12);

    double cdw = phi * k_uw * ((1 - rs) * (1 - rw)) / ((1 - rs) * rw);
    double cds = phi * k_ws * ((1 - rs) * rw) / rs;

    double k_wu = k_uw * (1 / rw - 1) - k_ws;
    double k_su = 3 * (k_ws * (1 / rs - 1) * rw);
    double Aw = (rs * TOT_A) / ((1 - rs) * rw + rs);
    double As = Aw;
    
    // XB model
    double lambda_value;
    lambda_value = std::min(1.2, lambda);
    double h_lambda_prima = 1 + beta_0 * (lambda_value + std::min(lambda_value, 0.87) - 1.87);
    double h_lambda = std::max(0.0, h_lambda_prima);

    double XU = (1 - TmBlocked) - XW - XS;

    // dXS and dXW operations
    double xb_ws = k_ws * XW;
    double xb_uw = k_uw * XU;
    double xb_wu = k_wu * XW;
    double xb_su = k_su * XS;

    double gamma_rate = gamma_S * std::max((ZETAS > 0) * ZETAS, (ZETAS < -1) * (-ZETAS - 1));
    double xb_su_gamma = gamma_rate * XS;

    double gamma_rate_w = gamma_wu * abs(ZETAW);
    double xb_wu_gamma = gamma_rate_w * XW;

    dY_land(7) = xb_ws - xb_su - xb_su_gamma;
    dY_land(8) = xb_uw - xb_wu - xb_ws - xb_wu_gamma;
    
    double min_val = std::min(lambda, 1.2);
    double ca50 = ca50_ref + beta_1 * (min_val - 1);
    dY_land(9) = k_trpn * (pow((Caii * 1000 / ca50), TRPN_n) * (1 - Ca_TRPN) - Ca_TRPN); 

    double ktm_block = ktm_unblock * pow(TRPN50, n_tm) / (1 - rs - (1 - rs) * rw);

    dY_land(10)= ktm_block * std::min(100.0, pow(Ca_TRPN, (-n_tm / 2))) * XU - ktm_unblock * pow(Ca_TRPN, (n_tm / 2)) * TmBlocked;

    // Velocity dependence
    dY_land(11) = As * dlambda_dt - cds * ZETAS;
    dY_land(12)= Aw * dlambda_dt - cdw * ZETAW;
    
    double Cd = Y_land(13);
    double C = lambda - 1;

    double eta;
    if (C - Cd < 0)
        eta = eta_s;
    else
        eta = eta_l;

    dY_land(13)= par_k * (C - Cd) / eta;

    double Fd = eta * dY_land(6);
    double F1 = (std::exp(b * C) - 1);
    double Tp = a * (F1 + Fd);

    // Active and Total Force
    double Ta = h_lambda * (Tref / rs) * ((ZETAS + 1) * XS + ZETAW * XW);

    T = Ta + Tp;

    // Compute Ca+ derivative value
    I_TRPN = dY_land(9) * TRPN_max;

}



