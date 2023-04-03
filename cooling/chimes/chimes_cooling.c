/****************************************************************************
 * This file is part of CHIMES.
 * Copyright (c) 2020 Alexander Richings (alexander.j.richings@durham.ac.uk)
 *
 * CHIMES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***************************************************************************/

#ifdef CHIMES_ENABLE_GNU_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 
#endif 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../allvars.h"
#include "chimes_vars.h"


#ifdef CHIMES // top-level flag

/** 
 * @brief Calculates the total number density. 
 * 
 * Sums the number densities of all ions and molecules 
 * in the network to get the total number density. 
 * 
 * @param my_abundances Abundance array. 
 * @param nH Total hydrogen number density. 
 * @param myGlobalVars The #globalVariables struct. 
 */ 
ChimesFloat calculate_total_number_density(ChimesFloat *my_abundances, ChimesFloat nH, struct globalVariables *myGlobalVars)
{
  int i;
  ChimesFloat result = 0.0f;
  
  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    result += my_abundances[i] * nH;
  
  return result;
}

/** 
 * @brief Calculates the mean molecular weight. 
 * 
 * Calculates the mean molecular weight from the abundance 
 * array. This is defined as: 
 * mu = (1.0 + sum(element_abundances[i] * Z[i])) / sum(abundances[i]). 
 * Where element_abundances[] are the abundances of He and the metals 
 * relative to H, Z[] are the atomic masses of each element, and 
 * abundances[] are the abundances of the individual ions and molecules 
 * relative to H. 
 * 
 * @param myGasVars The gasVariables struct. 
 * @param myGlobalVars The globalVariables struct. 
 */ 
ChimesFloat calculate_mean_molecular_weight(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  ChimesFloat denominator = 0.0;
  int i;
  
  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++)
    denominator += myGasVars->abundances[i];
  
  return (1.0 + myGasVars->element_abundances[0] * 4.0f + myGasVars->element_abundances[1] * 12.0f + myGasVars->element_abundances[2] * 14.0f + myGasVars->element_abundances[3] * 16.0f + myGasVars->element_abundances[4] * 20.0f + myGasVars->element_abundances[5] * 24.0f + myGasVars->element_abundances[6] * 28.0f + myGasVars->element_abundances[7] * 32.0f + myGasVars->element_abundances[8] * 40.0f + myGasVars->element_abundances[9] * 56.0f) / denominator;
}

/** 
 * @brief Compton cooling rate. 
 * 
 * Calculates the Compton cooling rate from the 
 * CMB radiation. It returns Lambda/nH^2, in 
 * units of erg cm^3 s^-1 
 * 
 * @param T Gas temperature. 
 * @param Tcmb CMB temperatures. 
 * @param xe Electron abundance. 
 * @param nH Total hydrogen number density. 
 */ 
ChimesFloat compton_cooling(ChimesFloat T, ChimesFloat Tcmb, ChimesFloat xe, ChimesFloat nH)
{
  return 1.017e-37 * Tcmb * Tcmb * Tcmb * Tcmb * (T - Tcmb) * xe / nH;	/* Lambda/nH^2 */
}

/** 
 * @brief OH rotational cooling. 
 * 
 * Calculates the cooling rate from rotational 
 * transitions of the OH molecule. This is calculated 
 * from Hollenbach and McKee 1979, ApJS, 41, 555 (see 
 * their equations 6.21 and 6.22, and the parameters 
 * in their Table 3). 
 * This function returns Lambda/nH^2, in units of 
 * erg cm^3 s^-1. 
 * 
 * @param myGasVars The #gasVariables struct. 
 * @param myGlobalVars The #globalVariables struct. 
 * @param data The #UserData struct. 
 */ 
ChimesFloat OH_rotational_cooling(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data)
{
  ChimesFloat dv, N_tau, tau_T, c_tau, n_cr, ym; 

  // Thermal velocity dispersion in cgs 
  dv = chimes_sqrt(3.0f * BOLTZMANNCGS * myGasVars->temperature / (PROTON_MASS * calculate_mean_molecular_weight(myGasVars, myGlobalVars))); 
  
  N_tau = 1.485e11f * dv;
  tau_T = 4.0f * data.OH_column / (10.0f * (myGasVars->temperature / 27.0f) * 6.8e-4f * N_tau);

  c_tau = tau_T * chimes_sqrt(2.0f * PI * chimes_log(2.13f + ((tau_T / EULERS_CONST) * (tau_T / EULERS_CONST)))) / ((chimes_exp(-data.extinction) / (1.0f + (data.extinction * data.extinction))) + 2.0f * data.extinction * chimes_sqrt(chimes_log(1.0f + (tau_T / EULERS_CONST))) * chimes_sqrt(chimes_log(tau_T / (data.extinction * EULERS_CONST))));

  if (isnan(c_tau) != 0)
    c_tau = 0.0f;

  n_cr = 1.5e10f * chimes_sqrt(myGasVars->temperature / 1.0e3f);
  ym = chimes_log(1.0f + (c_tau / (1.0f + 10.0f * (n_cr / myGasVars->nH_tot))));

  // Lambda / nH^2 (erg cm^3 s^-1) 
  return myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]] * (2.0f * BOLTZMANNCGS * myGasVars->temperature * myGasVars->temperature * 2.3e-2f / (myGasVars->nH_tot * 27.0f)) * ((2.0f + ym + 0.6f * ym * ym) / (1.0f + c_tau + (n_cr / myGasVars->nH_tot) + 1.5f * chimes_sqrt(n_cr / myGasVars->nH_tot))); 
}

/** 
 * @brief Update the cooling rates. 
 * 
 * Calculates the rates of the various cooling 
 * and heating channels and stores them in the 
 * #chimes_current_rates_struct struct within 
 * the #UserData struct. 
 *
 * @param myGasVars The #gasVariables struct. 
 * @param myGlobalVars The #globalVariables struct. 
 * @param data The #UserData struct. 
 */ 
void update_cooling_rates(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data) 
{ 
  int i, T_index, nHI_index, ne_index, nHII_index; 
  ChimesFloat dT, dnHI, dne, dnHII; 
  ChimesFloat log_T, log_nHI, log_ne, log_nHII; 

  const int N_T = chimes_table_bins.N_Temperatures; 
  log_T = (ChimesFloat) chimes_log10(myGasVars->temperature); 
  chimes_get_table_index(chimes_table_bins.Temperatures, N_T, log_T, &T_index, &dT);

  for (i = 0; i < chimes_table_cooling.N_coolants; i++) 
    data.chimes_current_rates->cooling_rate[i] = chimes_exp10(chimes_interpol_2d_fix_x(chimes_table_cooling.rates, i, T_index, dT, N_T)); 

  if (chimes_table_cooling.N_coolants_2d > 0) 
    {
      if (log_T < chimes_table_bins.cool_2d_Temperatures[chimes_table_bins.N_cool_2d_Temperatures - 1]) 
	{
	  const int N_Tcool_2d = chimes_table_bins.N_cool_2d_Temperatures;
	  const int N_ne_2d = chimes_table_bins.N_cool_2d_ElectronDensities; 
	  
	  log_ne = (ChimesFloat) chimes_log10(chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] * myGasVars->nH_tot, CHIMES_FLT_MIN)); 
	  chimes_get_table_index(chimes_table_bins.cool_2d_Temperatures, N_Tcool_2d, log_T, &T_index, &dT); 
	  chimes_get_table_index(chimes_table_bins.cool_2d_ElectronDensities, N_ne_2d, log_ne, &ne_index, &dne); 
      
	  for (i = 0; i < chimes_table_cooling.N_coolants_2d; i++) 
	    data.chimes_current_rates->cooling_rate_2d[i] = chimes_exp10(chimes_interpol_3d_fix_x(chimes_table_cooling.rates_2d, i, T_index, ne_index, dT, dne, N_Tcool_2d, N_ne_2d)); 
	}
      else 
	{
	  const int N_Tcool_hiT_2d = chimes_table_bins.N_cool_hiT_2d_Temperatures;
	  chimes_get_table_index(chimes_table_bins.cool_hiT_2d_Temperatures, N_Tcool_hiT_2d, log_T, &T_index, &dT); 
      
	  for (i = 0; i < chimes_table_cooling.N_coolants_2d; i++) 
	    data.chimes_current_rates->cooling_rate_2d[i] = chimes_exp10(chimes_interpol_2d_fix_x(chimes_table_cooling.rates_hiT_2d, i, T_index, dT, N_Tcool_hiT_2d)); 
	}
    }  
  
  if (chimes_table_cooling.N_coolants_4d > 0) 
    {
      if (log_T < chimes_table_bins.cool_4d_Temperatures[chimes_table_bins.N_cool_4d_Temperatures - 1]) 
	{
	  const int N_Tcool_4d = chimes_table_bins.N_cool_4d_Temperatures;
	  const int N_nHI_4d = chimes_table_bins.N_cool_4d_HIDensities;
	  const int N_ne_4d = chimes_table_bins.N_cool_4d_ElectronDensities;
	  const int N_nHII_4d = chimes_table_bins.N_cool_4d_HIIDensities; 
	  log_nHI = (ChimesFloat) chimes_log10(chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]] * myGasVars->nH_tot, CHIMES_FLT_MIN)); 
	  log_ne = (ChimesFloat) chimes_log10(chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]] * myGasVars->nH_tot, CHIMES_FLT_MIN)); 
	  log_nHII = (ChimesFloat) chimes_log10(chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]] * myGasVars->nH_tot, CHIMES_FLT_MIN)); 
	  
	  chimes_get_table_index(chimes_table_bins.cool_4d_Temperatures, N_Tcool_4d, log_T, &T_index, &dT); 
	  chimes_get_table_index(chimes_table_bins.cool_4d_HIDensities, N_nHI_4d, log_nHI, &nHI_index, &dnHI); 
	  chimes_get_table_index(chimes_table_bins.cool_4d_ElectronDensities, N_ne_4d, log_ne, &ne_index, &dne); 
	  chimes_get_table_index(chimes_table_bins.cool_4d_HIIDensities, N_nHII_4d, log_nHII, &nHII_index, &dnHII); 

	  for (i = 0; i < chimes_table_cooling.N_coolants_4d; i++) 
	    data.chimes_current_rates->cooling_rate_4d[i] = chimes_exp10(chimes_interpol_5d_fix_x(chimes_table_cooling.rates_4d, i, T_index, nHI_index, ne_index, nHII_index, dT, dnHI, dne, dnHII, N_Tcool_4d, N_nHI_4d, N_ne_4d, N_nHII_4d)); 
	}
      else 
	{
	  const int N_Tcool_hiT_4d = chimes_table_bins.N_cool_hiT_4d_Temperatures;
	  chimes_get_table_index(chimes_table_bins.cool_hiT_4d_Temperatures, N_Tcool_hiT_4d, log_T, &T_index, &dT); 
      
	  for (i = 0; i < chimes_table_cooling.N_coolants_4d; i++) 
	    data.chimes_current_rates->cooling_rate_4d[i] = chimes_exp10(chimes_interpol_2d_fix_x(chimes_table_cooling.rates_hiT_4d, i, T_index, dT, N_Tcool_hiT_4d)); 
	} 
    }
  
  return; 
}

/** 
 * @brief Calculates the total cooling rate. 
 * 
 * Sums all cooling and heating channels to calculate 
 * the total cooling rate. This function can return 
 * either the net cooling rate, cooling only, or heating 
 * only, depending on the mode argument. The rate is 
 * returned in units of erg cm^-3 s^-1. 
 * If the net cooling rate is returned, a positive value 
 * indicates net cooling, while a negative value indicates 
 * net heating. If only the cooling or heating rate is 
 * returned, the value will always be positive. 
 * 
 * @param myGasVars The #gasVariables struct. 
 * @param myGlobalVars The #globalVariables struct. 
 * @param data The #UserData struct. 
 * @param mode An integer flag. 0: return net cooling; 1: return cooling only; 2: return heating only. 
 */ 
ChimesFloat calculate_total_cooling_rate(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars, struct UserData data, int mode) 
{
  int i, xHII_index, Psi_index, T_index, T_mol_index, N_CO_rot_index, N_CO_vib_index; 
  int N_H2O_rot_index, N_H2O_vib_index, T_H2O_index; 
  ChimesFloat x_elec, log_xHII, d_xHII, log_Psi, dPsi, log_T, dT, G0; 
  ChimesFloat H2_lowDens, H2_LTE, dT_mol, xHI, xH2, H2_crit_density; 
  ChimesFloat dN_CO_rot, dN_CO_vib, log_N_eff; 
  ChimesFloat CO_rot_L0, CO_rot_Llte, CO_rot_nhalf, CO_rot_a, CO_rot_neff; 
  ChimesFloat CO_vib_L0, CO_vib_Llte, CO_vib_neff, H2O_rot_neff, H2O_vib_neff; 
  ChimesFloat dN_H2O_rot, dN_H2O_vib, dT_H2O; 
  ChimesFloat H2O_rot_L0, H2O_rot_Llte, H2O_rot_nhalf, H2O_rot_a, H2O_vib_L0, H2O_vib_Llte; 
  ChimesFloat cr_secondary; 
  ChimesFloat total_cooling = 0.0f; 
  ChimesFloat total_heating = 0.0f; 

  // Check that mode is valid.
  if (!((mode == 0) || (mode == 1) || (mode == 2)))
    {
      printf("CHIMES error: calculate_total_cooling_rate() called with mode == %d. Allowed modes: 0 (net cooling), 1 (cooling only), 2 (heating only).\n", mode); 
      chimes_exit();
    }
  
  update_cooling_rates(myGasVars, myGlobalVars, data); 
  
  x_elec = myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]]; 
  
  for (i = 0; i < chimes_table_cooling.N_coolants; i++) 
    total_cooling += data.chimes_current_rates->cooling_rate[i] * myGasVars->abundances[chimes_table_cooling.coolants[i]] * x_elec; 

  for (i = 0; i < chimes_table_cooling.N_coolants_2d; i++) 
    total_cooling += data.chimes_current_rates->cooling_rate_2d[i] * myGasVars->abundances[chimes_table_cooling.coolants_2d[i]] * x_elec; 
  
  for (i = 0; i < chimes_table_cooling.N_coolants_4d; i++) 
    total_cooling += data.chimes_current_rates->cooling_rate_4d[i] * myGasVars->abundances[chimes_table_cooling.coolants_4d[i]] / myGasVars->nH_tot; 

  // Photoheating 
  if (myGlobalVars->N_spectra > 0) 
    {
      for (i = 0; i < chimes_table_photoion_fuv.N_reactions[data.mol_flag_index]; i++) 
	total_heating += data.chimes_current_rates->photoion_fuv_heat_rate[i] * myGasVars->abundances[chimes_table_photoion_fuv.reactants[i]] / myGasVars->nH_tot; 

      for (i = 0; i < chimes_table_photoion_euv.N_reactions[data.mol_flag_index]; i++) 
	total_heating += data.chimes_current_rates->photoion_euv_heat_rate[i] * myGasVars->abundances[chimes_table_photoion_euv.reactants[i]] / myGasVars->nH_tot; 
    }

  // Cosmic ray heating 
  for (i = 0; i < chimes_table_cosmic_ray.N_reactions[data.mol_flag_index]; i++) 
    total_heating += 3.2e-11f * data.chimes_current_rates->cosmic_ray_rate[i] / myGasVars->nH_tot; 

  // Correct for secondary cosmic rays 
  log_xHII = chimes_log10(chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]], CHIMES_FLT_MIN));
  const int N_xHII = chimes_table_bins.N_secondary_cosmic_ray_xHII; 
  chimes_get_table_index(chimes_table_bins.secondary_cosmic_ray_xHII, N_xHII, log_xHII, &xHII_index, &d_xHII); 
  for (i = 0; i < 2; i++) 
    {
      cr_secondary = chimes_exp10(chimes_interpol_2d_fix_x(chimes_table_cosmic_ray.secondary_ratio, i, xHII_index, d_xHII, N_xHII)); 
      total_heating -= 3.2e-11f * data.chimes_current_rates->cosmic_ray_rate[chimes_table_cosmic_ray.secondary_base_reaction[i]] * cr_secondary / (myGasVars->nH_tot * (1.0f + cr_secondary)); 
    }

  // Compton cooling from the CMB 
  total_cooling += compton_cooling(myGasVars->temperature, myGlobalVars->cmb_temperature, myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]], myGasVars->nH_tot); 

  if (data.mol_flag_index == 1) 
    {
      // H2 rovibrational cooling 
      log_T = chimes_log10(myGasVars->temperature); 
      chimes_get_table_index(chimes_table_bins.mol_cool_Temperatures, chimes_table_bins.N_mol_cool_Temperatures, log_T, &T_mol_index, &dT_mol); 
      H2_lowDens = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_H2, T_mol_index, dT_mol)) * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]]; 
      H2_lowDens += chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_HI, T_mol_index, dT_mol)) * myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]]; 
      H2_lowDens += chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_HII, T_mol_index, dT_mol)) * myGasVars->abundances[myGlobalVars->speciesIndices[sp_HII]]; 
      H2_lowDens += chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_HeI, T_mol_index, dT_mol)) * myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI]]; 
      H2_lowDens += chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2_cool_lowDens_elec, T_mol_index, dT_mol)) * myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]]; 
      H2_lowDens *= myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]]; 

      if (H2_lowDens > 0.0f) 
	{
	  H2_LTE = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2_cool_LTE, T_mol_index, dT_mol)) * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]] / myGasVars->nH_tot; 
	  total_cooling += H2_LTE / (1.0f + (H2_LTE / H2_lowDens)); 
	}

      // H2 collis dissoc heating 
      xHI = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]]; 
      xH2 = myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]]; 
      
      total_cooling += 7.2e-12f * data.chimes_current_rates->H2_collis_dissoc_rate_coefficient[chimes_table_H2_collis_dissoc.Heating_reaction_index] * xHI * xH2;
      total_cooling += 7.2e-12f * data.chimes_current_rates->T_dependent_rate_coefficient[chimes_table_T_dependent.H2_collis_dissoc_heating_reaction_index] * xH2 * xH2; 

      // Gas-phase H2 formation 
      if (xHI + xH2 == 0) 
	H2_crit_density = 0.0f; 
      else 
	H2_crit_density = (xHI + xH2) / ((xHI / data.chimes_current_rates->H2_collis_dissoc_crit_H) + (xH2 / data.chimes_current_rates->H2_collis_dissoc_crit_H2)); 

      total_heating += (((2.93e-12f * data.chimes_current_rates->T_dependent_rate[chimes_table_T_dependent.H2_form_heating_reaction_index]) + (5.65e-12f * data.chimes_current_rates->constant_rate[chimes_table_constant.H2_form_heating_reaction_index])) * (1.0f / (myGasVars->nH_tot + H2_crit_density))); 

      // Dust-catalysed H2 formation 
      total_heating += 7.16e-12f * (data.chimes_current_rates->H2_dust_formation_rate / myGasVars->nH_tot) * (myGasVars->nH_tot / (myGasVars->nH_tot + H2_crit_density)); 

      // CO cooling 
      if ((myGlobalVars->element_included[0] == 1) && (myGlobalVars->element_included[2] == 1)) 
	{
	  // N_eff units: cm^-2 per km s^-1 
	  if (myGlobalVars->StaticMolCooling == 1)
	    log_N_eff = chimes_log10(chimes_max(1.0e5f * data.CO_column / (chimes_sqrt(3.0f * BOLTZMANNCGS * myGasVars->temperature / (PROTON_MASS * calculate_mean_molecular_weight(myGasVars, myGlobalVars)))), CHIMES_FLT_MIN));	
	  else
	    log_N_eff = chimes_log10(chimes_max(1.0e5f * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] * myGasVars->nH_tot / chimes_max(fabs(myGasVars->divVel), CHIMES_FLT_MIN), CHIMES_FLT_MIN)); 

	  const int N_COcool_rot_ColDens = chimes_table_bins.N_CO_cool_rot_ColumnDensities; 
	  const int N_COcool_vib_ColDens = chimes_table_bins.N_CO_cool_vib_ColumnDensities; 
	  chimes_get_table_index(chimes_table_bins.CO_cool_rot_ColumnDensities, N_COcool_rot_ColDens, log_N_eff, &N_CO_rot_index, &dN_CO_rot); 
	  chimes_get_table_index(chimes_table_bins.CO_cool_vib_ColumnDensities, N_COcool_vib_ColDens, log_N_eff, &N_CO_vib_index, &dN_CO_vib); 

	  CO_rot_L0 = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.CO_cool_rot_L0, T_mol_index, dT_mol)); 
	  CO_rot_Llte = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.CO_cool_rot_Llte, T_mol_index, N_CO_rot_index, dT_mol, dN_CO_rot, N_COcool_rot_ColDens)); 
	  CO_rot_nhalf = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.CO_cool_rot_nhalf, T_mol_index, N_CO_rot_index, dT_mol, dN_CO_rot, N_COcool_rot_ColDens)); 
	  CO_rot_a = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.CO_cool_rot_a, T_mol_index, N_CO_rot_index, dT_mol, dN_CO_rot, N_COcool_rot_ColDens)); 
      
	  CO_rot_neff = myGasVars->nH_tot * (xH2 + 9.857f * chimes_pow((myGasVars->temperature / 1.0e3f), 0.25f) * xHI + 680.13f * chimes_pow(myGasVars->temperature, -0.25f) * x_elec);
      
	  total_cooling += xH2 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] / ((1.0f / CO_rot_L0) + (CO_rot_neff / CO_rot_Llte) + (1.0f / CO_rot_L0) * chimes_pow((CO_rot_neff / CO_rot_nhalf), CO_rot_a) * (1.0f - (CO_rot_nhalf * CO_rot_L0 / CO_rot_Llte))); 

      
	  CO_vib_L0 = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.CO_cool_vib_L0, T_mol_index, dT_mol)); 
	  CO_vib_Llte = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.CO_cool_vib_Llte, T_mol_index, N_CO_vib_index, dT_mol, dN_CO_vib, N_COcool_vib_ColDens)); 
      
	  CO_vib_neff = myGasVars->nH_tot * (xH2 + 50.0f * xHI + 9035.09f * chimes_exp(68.0f / chimes_pow(myGasVars->temperature, 1.0f / 3.0f)) * chimes_pow(myGasVars->temperature / 300.0f, 0.938f) * x_elec); 

	  total_cooling += xH2 * myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]] / ((1.0f / CO_vib_L0) + (CO_vib_neff / CO_vib_Llte)); 
	}

      // H2O and OH cooling 
      if (myGlobalVars->element_included[2] == 1) 
	{
	  if (myGlobalVars->StaticMolCooling == 1)
	    log_N_eff = (ChimesFloat) chimes_log10(chimes_max(1.0e5f * data.H2O_column / (chimes_sqrt(3.0f * BOLTZMANNCGS * myGasVars->temperature / (PROTON_MASS * calculate_mean_molecular_weight(myGasVars, myGlobalVars)))), CHIMES_FLT_MIN)); 
	  else
	    log_N_eff = (ChimesFloat) chimes_log10(chimes_max(1.0e5f * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] * myGasVars->nH_tot / chimes_max(fabs(myGasVars->divVel), CHIMES_FLT_MIN), CHIMES_FLT_MIN)); 
	  
	  // rotational cooling
	  const int N_H2Ocool_rot_ColDens = chimes_table_bins.N_H2O_cool_rot_ColumnDensities; 
	  chimes_get_table_index(chimes_table_bins.H2O_cool_rot_ColumnDensities, N_H2Ocool_rot_ColDens, log_N_eff, &N_H2O_rot_index, &dN_H2O_rot);

	  H2O_rot_neff = myGasVars->nH_tot * (xH2 + 10.0f * xHI + chimes_exp10((-8.02f + (15.749f / chimes_pow(myGasVars->temperature, 1.0f / 6.0f)) - (47.137f / chimes_pow(myGasVars->temperature, 1.0f / 3.0f)) + (76.648f / chimes_sqrt(myGasVars->temperature)) - (60.191f / chimes_pow(myGasVars->temperature, 2.0f / 3.0f)))) * x_elec / (7.4e-12f * chimes_sqrt(myGasVars->temperature))); 
	  if (log_T >= 2.0f) 
	    {
	      chimes_get_table_index(chimes_table_bins.H2O_cool_hiT_Temperatures, chimes_table_bins.N_H2O_cool_hiT_Temperatures, log_T, &T_H2O_index, &dT_H2O); 
	      
	      H2O_rot_L0 = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2O_cool_rot_hiT_L0, T_H2O_index, dT_H2O)); 
	      H2O_rot_Llte = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2O_cool_rot_hiT_Llte, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 
	      H2O_rot_nhalf = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2O_cool_rot_hiT_nhalf, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 
	      H2O_rot_a = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2O_cool_rot_hiT_a, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 

	      total_cooling += myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] * xH2 / ((1.0f / H2O_rot_L0) + (H2O_rot_neff / H2O_rot_Llte) + (1.0f / H2O_rot_L0) * chimes_pow((H2O_rot_neff / H2O_rot_nhalf), H2O_rot_a) * (1.0f - (H2O_rot_nhalf * H2O_rot_L0 / H2O_rot_Llte))); 
	    }
	  else
	    {
	      chimes_get_table_index(chimes_table_bins.H2O_cool_lowT_Temperatures, chimes_table_bins.N_H2O_cool_lowT_Temperatures, log_T, &T_H2O_index, &dT_H2O); 
	      
	      H2O_rot_L0 = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2Oortho_cool_rot_lowT_L0, T_H2O_index, dT_H2O)); 
	      H2O_rot_Llte = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2Oortho_cool_rot_lowT_Llte, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 
	      H2O_rot_nhalf = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2Oortho_cool_rot_lowT_nhalf, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 
	      H2O_rot_a = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2Oortho_cool_rot_lowT_a, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 

	      total_cooling += 0.75f * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] * xH2 / ((1.0f / H2O_rot_L0) + (H2O_rot_neff / H2O_rot_Llte) + (1.0f / H2O_rot_L0) * chimes_pow((H2O_rot_neff / H2O_rot_nhalf), H2O_rot_a) * (1.0f - (H2O_rot_nhalf * H2O_rot_L0 / H2O_rot_Llte))); 

	      H2O_rot_L0 = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2Opara_cool_rot_lowT_L0, T_H2O_index, dT_H2O)); 
	      H2O_rot_Llte = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2Opara_cool_rot_lowT_Llte, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 
	      H2O_rot_nhalf = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2Opara_cool_rot_lowT_nhalf, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 
	      H2O_rot_a = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2Opara_cool_rot_lowT_a, T_H2O_index, N_H2O_rot_index, dT_H2O, dN_H2O_rot, N_H2Ocool_rot_ColDens)); 

	      total_cooling += 0.25f * myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] * xH2 / ((1.0f / H2O_rot_L0) + (H2O_rot_neff / H2O_rot_Llte) + (1.0f / H2O_rot_L0) * chimes_pow((H2O_rot_neff / H2O_rot_nhalf), H2O_rot_a) * (1.0f - (H2O_rot_nhalf * H2O_rot_L0 / H2O_rot_Llte))); 
	    }

	  // vibrational cooling 
	  const int N_H2Ocool_vib_ColDens = chimes_table_bins.N_H2O_cool_vib_ColumnDensities; 
	  chimes_get_table_index(chimes_table_bins.H2O_cool_vib_ColumnDensities, chimes_table_bins.N_H2O_cool_vib_ColumnDensities, log_N_eff, &N_H2O_vib_index, &dN_H2O_vib);

	  H2O_vib_L0 = chimes_exp10(chimes_interpol_1d(chimes_table_cooling.H2O_cool_vib_L0, T_H2O_index, dT_H2O)); 
	  H2O_vib_Llte = chimes_exp10(chimes_interpol_2d(chimes_table_cooling.H2O_cool_vib_Llte, T_H2O_index, N_H2O_vib_index, dT_H2O, dN_H2O_vib, N_H2Ocool_vib_ColDens)); 
	  
	  H2O_vib_neff = myGasVars->nH_tot * (xH2 + 10.0f * xHI + 4.0625e8f * chimes_exp(47.5f / chimes_pow(myGasVars->temperature, 1.0f / 3.0f)) * x_elec / chimes_sqrt(myGasVars->temperature));

	  total_cooling += myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]] * xH2 / ((1.0f / H2O_vib_L0) + (H2O_vib_neff / H2O_vib_Llte)); 

	  total_cooling += OH_rotational_cooling(myGasVars, myGlobalVars, data); 
	}


      chimes_get_table_index(chimes_table_bins.Temperatures, chimes_table_bins.N_Temperatures, log_T, &T_index, &dT);

      if (myGlobalVars->N_spectra > 0) 
	{
	  // H2 photodissoc heating 
	  total_heating += 6.4e-13f * data.chimes_current_rates->H2_photodissoc_rate[0] / myGasVars->nH_tot; 
	  
	  // H2 UV pumping 
	  total_heating += 2.7e-11f * data.chimes_current_rates->H2_photodissoc_rate[0] * (1.0f / (myGasVars->nH_tot + H2_crit_density)); 

	  // Photoelectric dust heating & grain recombination cooling 
	  // Note that dust processes are only included when 
	  // the molecular network is switched on. 
	  G0 = 0.0f; 
	  for (i = 0; i < myGlobalVars->N_spectra; i++) 
	    G0 += myGasVars->isotropic_photon_density[i] * LIGHTSPEED * myGasVars->G0_parameter[i]; 

	  G0 *= chimes_exp(-data.extinction * G0_GAMMA); 
      
	  if ((G0 > 0.0f) && (x_elec > 0.0f)) 
	    {
	      // Include a factor phi_pah = 0.5 in the denominator 
	      log_Psi = (ChimesFloat) chimes_log10(chimes_max(G0 * chimes_sqrt(myGasVars->temperature) / chimes_max(myGasVars->nH_tot * x_elec * 0.5f, CHIMES_FLT_MIN), CHIMES_FLT_MIN));

	      chimes_get_table_index(chimes_table_bins.Psi, chimes_table_bins.N_Psi, log_Psi, &Psi_index, &dPsi); 
	  
	      total_heating += chimes_exp10(chimes_interpol_2d(chimes_table_cooling.photoelectric_heating, T_index, Psi_index, dT, dPsi, chimes_table_bins.N_Psi)) * G0 * myGasVars->dust_ratio / myGasVars->nH_tot; 
	      total_cooling += chimes_exp10(chimes_interpol_2d(chimes_table_cooling.grain_recombination, T_index, Psi_index, dT, dPsi, chimes_table_bins.N_Psi)) * myGasVars->dust_ratio * myGasVars->abundances[myGlobalVars->speciesIndices[sp_elec]]; 
	    }
	}

      // Gas-grain energy transfer 
      total_cooling += chimes_exp10(chimes_interpol_1d(chimes_table_cooling.gas_grain_transfer, T_index, dT)) * myGasVars->dust_ratio * (myGasVars->temperature - myGlobalVars->grain_temperature); 
    }

  if (myGlobalVars->hybrid_cooling_mode == 1) 
    {
      if (myGlobalVars->hybrid_cooling_fn == NULL)
	{
	  printf("CHIMES error: hybrid_cooling_mode == %d, but the hybrid_cooling_fn has not been set. \n", myGlobalVars->hybrid_cooling_mode);
	  chimes_exit();
	}

      total_heating += (ChimesFloat) (*myGlobalVars->hybrid_cooling_fn)(myGasVars, myGlobalVars); 
    }

  // Convert units to erg/cm3/s 
  total_cooling *= myGasVars->nH_tot * myGasVars->nH_tot; 
  total_heating *= myGasVars->nH_tot * myGasVars->nH_tot; 
  total_heating += myGasVars->constant_heating_rate; 

  if (mode == 0) 
    return total_cooling - total_heating; 
  else if (mode == 1) 
    return total_cooling; 
  else  // mode == 2
    return total_heating; 
}

/** 
 * @brief Evolve the cooling for equilibrium abundances. 
 * 
 * If the abundances have been set to equilibrium, from the 
 * pre-computed equilibrium tables, and thermal evolution 
 * is switched on, then this routine is used to evolve 
 * the temperature, using a bisection integration method. 
 * As the temperature evolves, the equilibrium abundances 
 * are re-computed from the tables, and are then used 
 * to update the net cooling rate. 
 * 
 * @param data The #UserData struct. 
 */
void do_equilibrium_cooling(struct UserData data)
{
  int i, maxIter;
  ChimesFloat u, u_old, u_upper, u_lower, du, LambdaNet, dt;

  dt = data.myGasVars->hydro_timestep;

  data.myGasVars->temperature = chimes_max(data.myGasVars->temperature, data.myGasVars->TempFloor);

  set_equilibrium_abundances_from_tables(data);

  update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
  update_rates(data.myGasVars, data.myGlobalVars, data);

  LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 

  if (data.myGasVars->temperature <= data.myGasVars->TempFloor && LambdaNet <= 0.0f)
    return; 

  u = data.myGasVars->temperature * 1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS;
  u_old = u;
  u_upper = u;
  u_lower = u;

  /* If the cooling rate is small, take explicit solution. */
  if (fabs(LambdaNet * dt) < 0.10f * u_old)
    {
      u = u_old + LambdaNet * dt;
      data.myGasVars->temperature = chimes_max(u / (1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS), data.myGasVars->TempFloor);
      set_equilibrium_abundances_from_tables(data);

      // Check that explicit solution is valid
      update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
      update_rates(data.myGasVars, data.myGlobalVars, data);
      
      LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 

      if (fabs(LambdaNet * dt) < 0.10f * u_old && fabs(LambdaNet * dt) < 0.10f * u)
	return;
      else
	{
	  /* New cooling rate has increased, and explicit 
	   * solution is no longer valid. Reset and 
	   * continue with implicit solution. */
	  u = u_old; 
	  data.myGasVars->temperature = chimes_max(u / (1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS), data.myGasVars->TempFloor);
	  set_equilibrium_abundances_from_tables(data);
	}
    }

  i = 0;
  maxIter = 150;

  /* Bracketing */
  if (u - u_old - LambdaNet * dt < 0.0f) /* heating */
    {
      u_upper *= chimes_sqrt(1.2f);
      u_lower /= chimes_sqrt(1.2f);

      data.myGasVars->temperature = u_upper / (1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS);
      set_equilibrium_abundances_from_tables(data);

      update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
      update_rates(data.myGasVars, data.myGlobalVars, data);
      
      LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 

      while (u_upper - u_old - LambdaNet * dt < 0.0f && i < maxIter)
	{
	  u_upper *= 1.2f;
	  u_lower *= 1.2f;

	  data.myGasVars->temperature = u_upper / (1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS);
	  set_equilibrium_abundances_from_tables(data);

	  update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
	  update_rates(data.myGasVars, data.myGlobalVars, data);
	  
	  LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 

	  i++;
	}
	  
      if (i == maxIter)
	printf("WARNING: Problem with eqm cooling finding the upper bound.\n");
    }
  else /* cooling */ 
    {
      u_upper *= chimes_sqrt(1.2f);
      u_lower /= chimes_sqrt(1.2f);

      data.myGasVars->temperature = u_lower / (1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS);
      if (data.myGasVars->temperature <= data.myGasVars->TempFloor)
	{
	  data.myGasVars->temperature = data.myGasVars->TempFloor; 
	  u_lower = data.myGasVars->TempFloor * 1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS;
	  set_equilibrium_abundances_from_tables(data);

	  update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
	  update_rates(data.myGasVars, data.myGlobalVars, data);
	  
	  LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 
	}
      else
	{
	  set_equilibrium_abundances_from_tables(data);

	  update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
	  update_rates(data.myGasVars, data.myGlobalVars, data);
	  
	  LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 
	  
	  while (u_lower - u_old - LambdaNet * dt > 0.0f && i < maxIter)
	    {
	      u_upper /= 1.2f;
	      u_lower /= 1.2f;

	      data.myGasVars->temperature = u_lower / (1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS);
	      if (data.myGasVars->temperature <= data.myGasVars->TempFloor)
		{
		  data.myGasVars->temperature = data.myGasVars->TempFloor;
		  u_lower = data.myGasVars->TempFloor * 1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS;
		  set_equilibrium_abundances_from_tables(data);

		  update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
		  update_rates(data.myGasVars, data.myGlobalVars, data);
	  
		  LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 
		  break;
		}

	      set_equilibrium_abundances_from_tables(data);

	      update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
	      update_rates(data.myGasVars, data.myGlobalVars, data);
	      
	      LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 

	      i++;
	    }
	  if (i == maxIter)
	    printf("WARNING: Problem with eqm cooling finding the lower bound.\n");
	}
      if (u_lower - u_old - LambdaNet * dt > 0.0f && i < maxIter)
	return; /* u_lower reached TempFloor, but is still above converged solution. */
    }

  /* Iterate to convergence */
  i = 0;
  
  do
    {
      u = 0.5f * (u_lower + u_upper);
	  
      data.myGasVars->temperature = u / (1.5f * calculate_total_number_density(data.myGasVars->abundances, data.myGasVars->nH_tot, data.myGlobalVars) * BOLTZMANNCGS);
      set_equilibrium_abundances_from_tables(data);

      update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn);
      update_rates(data.myGasVars, data.myGlobalVars, data);
	      
      LambdaNet = -calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 

      if (u - u_old - LambdaNet * dt > 0.0f)
	u_upper = u;
      else
	u_lower = u;

      du = u_upper - u_lower;
      i++;
    }
  while (fabs(du / u) > 1.0e-6f && i < maxIter);

  if (i >= maxIter)
    printf("WARNING: eqm cooling failed to converge.\n");

  return;
}



#endif

