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
#include <sys/types.h>
#include <time.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h> 
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h> 
#include <sunmatrix/sunmatrix_dense.h> 
#include "../../allvars.h"
#include "chimes_vars.h"


#ifdef CHIMES


/** 
 * @brief Sets equilibrium abundances. 
 * 
 * Sets the abundances to their equilibrium values, taken 
 * from the pre-computed equilibrium abundance tables. 
 * 
 * @param data The #UserData struct containing the input data. 
 */ 
void set_equilibrium_abundances_from_tables(struct UserData data)
{
  // This is used when ForceEqOn == 1 
  int T_index, nH_index, Z_index, i;
  ChimesFloat dT, dnH, dZ;

  const int N_T = chimes_table_eqm_abundances.N_Temperatures;
  const int N_nH = chimes_table_eqm_abundances.N_Densities;
  const int N_Z = chimes_table_eqm_abundances.N_Metallicities; 
  chimes_get_table_index(chimes_table_eqm_abundances.Temperatures, N_T, chimes_log10(data.myGasVars->temperature), &T_index, &dT); 
  chimes_get_table_index(chimes_table_eqm_abundances.Densities, N_nH, chimes_log10(data.myGasVars->nH_tot), &nH_index, &dnH); 
  chimes_get_table_index(chimes_table_eqm_abundances.Metallicities, N_Z, chimes_log10(chimes_max(data.myGasVars->metallicity, CHIMES_FLT_MIN)), &Z_index, &dZ); 

  /* Note that the equilibrium tables tabulate
   * ionisation (or molecular) fraction, and 
   * NOT the abundance wrt H. Now we need to 
   * multiply by the appropriate element abundance. */
  for (i = 0; i < data.myGlobalVars->totalNumberOfSpecies; i++) 
    data.myGasVars->abundances[i] = chimes_exp10(chimes_interpol_4d_fix_x(chimes_table_eqm_abundances.Abundances, i, T_index, nH_index, Z_index, dT, dnH, dZ, N_T, N_nH, N_Z)) * data.species[i].element_abundance;

  /* Enforce constraint equations. */
  check_constraint_equations(data.myGasVars, data.myGlobalVars);

  return;
}

/** 
 * @brief Prints the gasVariables struct. 
 * 
 * Prints everything in the #gasVariables struct. 
 * 
 * @param log_file Output file to print to (typically, you would set this to stderr). 
 * @param myGasVars The #gasVariables struct. 
 * @param myGlobalVars The #globalVariables struct. 
 */ 
void chimes_print_gas_vars(FILE *log_file, struct gasVariables *myGasVars, struct globalVariables *myGlobalVars) 
{
  int i; 

  fprintf(log_file, "**************\n"); 
  fprintf(log_file, "ChimesGasVars:\n"); 
  fprintf(log_file, "**************\n"); 

  for (i = 0; i < 10; i++) 
    fprintf(log_file, "element_abundances[%d] = %.6e \n", i, myGasVars->element_abundances[i]); 

  fprintf(log_file, "nH_tot = %.6e \n", myGasVars->nH_tot); 
  fprintf(log_file, "temperature = %.6e \n", myGasVars->temperature); 
  fprintf(log_file, "TempFloor = %.6e \n", myGasVars->TempFloor); 
  fprintf(log_file, "divVel = %.6e \n", myGasVars->divVel); 
  fprintf(log_file, "doppler_broad = %.6e \n", myGasVars->doppler_broad); 

  for (i = 0; i < myGlobalVars->N_spectra; i++) 
    {
      fprintf(log_file, "isotropic_photon_density[%d] = %.6e \n", i, myGasVars->isotropic_photon_density[i]); 
      fprintf(log_file, "G0_parameter[%d] = %.6e \n", i, myGasVars->G0_parameter[i]); 
      fprintf(log_file, "H2_dissocJ[%d] = %.6e \n", i, myGasVars->H2_dissocJ[i]); 
    }

  fprintf(log_file, "cr_rate = %.6e \n", myGasVars->cr_rate); 
  fprintf(log_file, "metallicity = %.6e \n", myGasVars->metallicity); 
  fprintf(log_file, "dust_ratio = %.6e \n", myGasVars->dust_ratio); 
  fprintf(log_file, "cell_size = %.6e \n", myGasVars->cell_size); 
  fprintf(log_file, "hydro_timestep = %.6e \n", myGasVars->hydro_timestep); 
  fprintf(log_file, "ForceEqOn = %d \n", myGasVars->ForceEqOn); 
  fprintf(log_file, "ThermEvolOn = %d \n", myGasVars->ThermEvolOn); 
  fprintf(log_file, "InitIonState = %d \n", myGasVars->InitIonState); 
  fprintf(log_file, "constant_heating_rate = %.6e \n", myGasVars->constant_heating_rate); 

  for (i = 0; i < myGlobalVars->totalNumberOfSpecies; i++) {fprintf(log_file, "abundances[%d] = %.6e \n", i, myGasVars->abundances[i]);}
  fprintf(log_file, "++++++++++++++\n"); 
}

/** 
 * @brief Error handler function.
 * 
 * Error handler function for the 
 * CVODE error and warning messages. 
 * 
 * @param error_code The error code.
 * @param module CVODE module. 
 * @param function Function where the error occurred. 
 * @param msg The error message. 
 * @param eh_data Pointer to the user data. 
 */ 
void chimes_err_handler_fn(int error_code, const char *module, const char *function, char *msg, void *eh_data)
{
  struct UserData *user_data;

  user_data = (struct UserData *) eh_data;

  if (user_data->myGlobalVars->chimes_debug != 0)
    {
      if (!((user_data->myGasVars->temp_floor_mode == 1) && ((error_code == CV_RHSFUNC_FAIL) || (error_code == CV_LSETUP_FAIL))))
	{
	  cvErrHandler(error_code, module, function, msg, user_data->cvode_mem);

	  if (user_data->myGlobalVars->chimes_debug == 2)
	    {
	      fprintf(stderr, "CHIMES CVode error occurred for the following particle: \n"); 
	      chimes_print_gas_vars(stderr, user_data->myGasVars, user_data->myGlobalVars);
	    }
	}
    }

  return;
}

/** 
 * @brief Evolves the CHIMES network. 
 * 
 * This is the main CHIMES routine that actually integrates 
 * the chemical abundances and, if required, the temperature. 
 * 
 * @param myGasVars The #gasVariables struct. 
 * @param myGlobalVars The #globalVariables struct. 
 */ 
void chimes_network(struct gasVariables *myGasVars, struct globalVariables *myGlobalVars)
{
  realtype reltol, abstol_scalar, t;
  N_Vector abstol_vector, y;
  void *cvode_mem;

  ChimesFloat internal_energy;
  int total_network_size, nonmolecular_network_size, i, j;
  struct Species_Structure species[myGlobalVars->totalNumberOfSpecies];  
  struct UserData data;

  struct chimes_current_rates_struct chimes_current_rates; 
  allocate_current_rates_memory(&chimes_current_rates, myGlobalVars); 

  set_species_structures(species, myGasVars, &total_network_size, &nonmolecular_network_size, myGlobalVars);

  /* Set up structure to pass user
   * data to the solver. */
  data.myGasVars = myGasVars;
  data.myGlobalVars = myGlobalVars;
  data.species = species; 
  data.chimes_current_rates = &chimes_current_rates; 

  if (myGasVars->temperature <= myGlobalVars->T_mol)
    {
      data.mol_flag_index = 1; 
      data.network_size = total_network_size;
    }
  else
    {
      zero_molecular_abundances(species, myGasVars, myGlobalVars); 
      data.mol_flag_index = 0; 
      data.network_size = nonmolecular_network_size;
    }

  /* Enforce constraint equations. */
  check_constraint_equations(myGasVars, myGlobalVars);

  if (myGlobalVars->cellSelfShieldingOn > 0)
    {
      data.HI_column = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HI]] * myGasVars->cell_size * myGasVars->nH_tot;
      data.H2_column = myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2]] * myGasVars->cell_size * myGasVars->nH_tot;
      data.HeI_column = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeI]] * myGasVars->cell_size * myGasVars->nH_tot;
      data.HeII_column = myGasVars->abundances[myGlobalVars->speciesIndices[sp_HeII]] * myGasVars->cell_size * myGasVars->nH_tot;
      if (myGlobalVars->speciesIndices[sp_CO] > -1) 
	data.CO_column = chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_CO]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else 
	data.CO_column = 0.0f; 
      if (myGlobalVars->speciesIndices[sp_H2O] > -1) 
	data.H2O_column = chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_H2O]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else 
	data.H2O_column = 0.0f; 
      if (myGlobalVars->speciesIndices[sp_OH] > -1) 
	data.OH_column = chimes_max(myGasVars->abundances[myGlobalVars->speciesIndices[sp_OH]], 0.0) * myGasVars->cell_size * myGasVars->nH_tot;
      else 
	data.OH_column = 0.0f; 
      data.extinction = DUSTEFFSIZE * myGasVars->cell_size * myGasVars->nH_tot * myGasVars->dust_ratio;
    }
  else
    {
      data.HI_column = 0.0f;
      data.H2_column = 0.0f;
      data.HeI_column = 0.0f;
      data.HeII_column = 0.0f;
      data.CO_column = 0.0f;
      data.H2O_column = 0.0f;
      data.OH_column = 0.0f;
      data.extinction = 0.0f;
    }

  /* To determine whether to use case A or 
   * case B recombination, consider tau_HI 
   * and tau_HeI. Cross sections are taken 
   * from Verner et al. (1996). */ 
  if ((6.3463e-18f * data.HI_column) < 1.0f) 
    data.case_AB_index[0] = 0; 
  else 
    data.case_AB_index[0] = 1; 

  if ((7.4347e-18f * data.HeI_column) < 1.0f) 
    data.case_AB_index[1] = 0; 
  else 
    data.case_AB_index[1] = 1; 

  set_initial_rate_coefficients(myGasVars, myGlobalVars, data); 

  if (myGasVars->ForceEqOn == 1)
    {
      if (myGasVars->ThermEvolOn == 0)
	set_equilibrium_abundances_from_tables(data);
      else	  
	do_equilibrium_cooling(data); 

      free_current_rates_memory(&chimes_current_rates, myGlobalVars); 

      return;
    }


  /***************************** 
   * Try the explicit solution * 
   *****************************/ 

  // Update rates 
  int indices[CHIMES_TOTSIZE]; 
  update_rate_coefficients(data.myGasVars, data.myGlobalVars, data, data.myGasVars->ThermEvolOn); 
  update_rates(data.myGasVars, data.myGlobalVars, data); 
  
  // Zero all species rates 
  i = 0; 
  for (j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
    {
      if (data.species[j].include_species == 1)
	{
	  data.species[i].creation_rate = 0.0f;
	  data.species[i].destruction_rate = 0.0f;
	  indices[i] = j; 
	  i++; 
	}
    }

  // Compute creation and destruction rates 
  update_rate_vector(data.species, data.myGasVars, data.myGlobalVars, data); 
  
  ChimesFloat new_abundances[CHIMES_TOTSIZE]; 
  ChimesFloat old_energy, cool_rate, relative_change, this_absolute_tolerance; 
  ChimesFloat new_energy = 0.0f; 
  ChimesFloat max_relative_change = 0.0f; 
  
  for (i = 0; i < data.network_size; i++) 
    {
      new_abundances[indices[i]] = myGasVars->abundances[indices[i]] + ((data.species[indices[i]].creation_rate - data.species[indices[i]].destruction_rate) * myGasVars->hydro_timestep); 

      if (myGlobalVars->scale_metal_tolerances == 1) 
	this_absolute_tolerance = myGlobalVars->absoluteTolerance * data.species[indices[i]].element_abundance; 
      else 
	this_absolute_tolerance = myGlobalVars->absoluteTolerance; 
      
      if ((new_abundances[indices[i]] > this_absolute_tolerance) || (myGasVars->abundances[indices[i]] > this_absolute_tolerance)) 
	{
	  relative_change = fabs(new_abundances[indices[i]] - myGasVars->abundances[indices[i]]) / chimes_max(myGasVars->abundances[indices[i]], CHIMES_FLT_MIN); 
	  if (relative_change > max_relative_change) 
	    max_relative_change = relative_change; 
	}
    }

  if (data.myGasVars->ThermEvolOn == 1) 
    {
      if (data.myGasVars->temperature > data.myGasVars->TempFloor) 
	cool_rate = calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0); 
      else 
	cool_rate = chimes_min(calculate_total_cooling_rate(data.myGasVars, data.myGlobalVars, data, 0), 0.0f); 

      old_energy = myGasVars->temperature * 1.5f * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS;

      new_energy = old_energy - (cool_rate * myGasVars->hydro_timestep); 

      relative_change = fabs(new_energy - old_energy) / chimes_max(old_energy, CHIMES_FLT_MIN); 
      if (relative_change > max_relative_change) 
	max_relative_change = relative_change; 
    } 

  if (max_relative_change < myGlobalVars->explicitTolerance) 
    { 
      for (i = 0; i < data.network_size; i++) 
	myGasVars->abundances[indices[i]] = new_abundances[indices[i]]; 

      /* Enforce constraint equations. */
      check_constraint_equations(myGasVars, myGlobalVars);

      if (data.myGasVars->ThermEvolOn == 1) 
	myGasVars->temperature = chimes_max(new_energy / (1.5f * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS), myGasVars->TempFloor); 
      
      free_current_rates_memory(&chimes_current_rates, myGlobalVars); 

      return; 
    }
  else 
    {
      /************************************** 
       * Explicit solution is insufficient. * 
       * Use implicit solver.               *
       **************************************/ 
      
      /* Create a serial vector of length network_size
       * for the initial conditions. */
      if (myGasVars->ThermEvolOn == 0)
	{
	  y = N_VNew_Serial(data.network_size);
	  abstol_vector = N_VNew_Serial(data.network_size);
	}
      else
	{
	  y = N_VNew_Serial(data.network_size + 1);
	  internal_energy = myGasVars->temperature * 1.5f * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS;
	  NV_Ith_S(y, data.network_size) = (realtype) internal_energy;
	  abstol_vector = N_VNew_Serial(data.network_size + 1);

	  /* For the integration of the thermal energy,
	   * set the absolute tolerance to the minimum 
	   * float value. */ 
	  NV_Ith_S(abstol_vector, data.network_size) = (realtype) CHIMES_FLT_MIN;
	}
      
      i = 0;
      for (j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
	{
	  if (species[j].include_species == 1)
	    {  
	      NV_Ith_S(y, i) = (realtype) myGasVars->abundances[j];

	      if (myGlobalVars->scale_metal_tolerances == 1) 
		NV_Ith_S(abstol_vector, i) = (realtype) (myGlobalVars->absoluteTolerance * species[j].element_abundance); 
	      else if (myGasVars->ThermEvolOn == 1)
		NV_Ith_S(abstol_vector, i) = (realtype) myGlobalVars->absoluteTolerance;

	      i++;
	    }
	}

      /* Set up the solver */
      /* Set the tolerances*/
      reltol = (realtype) myGlobalVars->relativeTolerance;
      abstol_scalar = (realtype) myGlobalVars->absoluteTolerance;
    
      /* Use CVodeCreate to create the solver 
       * memory and specify the Backward Differentiation
       * Formula. Note that CVODE now uses Newton iteration
       * iteration by default, so no need to specify this. */
      cvode_mem = CVodeCreate(CV_BDF);
      data.cvode_mem = cvode_mem;

      /* Set the user data for CVode */
      CVodeSetUserData(cvode_mem, &data);
      
      /* Use CVodeSetMaxNumSteps to set the maximum number
       * of steps CVode takes. */      
      CVodeSetMaxNumSteps(cvode_mem, MAXSTEPS);

      /* Set the error handler function. */ 
      CVodeSetErrHandlerFn(cvode_mem, chimes_err_handler_fn, &data); 

      /* Use CVodeInit to initialise the integrator 
       * memory and specify the right hand side 
       * function in y' = f(t,y) (i.e. the rate
       * equations), the initial time 0.0 and the 
       * initial conditions, in y. */  
      CVodeInit(cvode_mem, f, 0.0f, y);

      /* Use CVodeSVtolerances to specify the scalar
       * relative and absolute tolerances. */
      if ((myGasVars->ThermEvolOn == 0) && (myGlobalVars->scale_metal_tolerances == 0))
	CVodeSStolerances(cvode_mem, reltol, abstol_scalar);
      else
	CVodeSVtolerances(cvode_mem, reltol, abstol_vector); 

      /* Create a dense SUNMatrix to use in the 
       * linear solver. */ 
      SUNMatrix A_sun; 

      if (myGasVars->ThermEvolOn == 0)
	A_sun = SUNDenseMatrix(data.network_size, data.network_size);
      else 
	A_sun = SUNDenseMatrix(data.network_size + 1, data.network_size + 1);

      /* Create a denst SUNLinearSolver object 
       * to use in CVode. */ 
      SUNLinearSolver LS_sun; 
      LS_sun = SUNLinSol_Dense(y, A_sun);

      /* Attach the matrix and linear 
       * solver to CVode. */ 
      CVodeSetLinearSolver(cvode_mem, LS_sun, A_sun);
      
      /* Specify the maximum number of convergence 
       * test failures. */
      CVodeSetMaxConvFails(cvode_mem, 5000);

      /* Call CVode() to integrate the chemistry. */ 
      CVode(cvode_mem, (realtype) myGasVars->hydro_timestep, y, &t, CV_NORMAL);

      /* Write the output abundances to the gas cell 
       * Note that species not included in the reduced
       * network are kept constant in the GasVars struct. */
      i = 0;
      for (j = 0; j < myGlobalVars->totalNumberOfSpecies; j++)
	{
	  if (species[j].include_species == 1)
	    {
	      myGasVars->abundances[j] = (ChimesFloat) NV_Ith_S(y, i);
	      i++;
	    }
	}
      
      /* Enforce constraint equations. */
      check_constraint_equations(myGasVars, myGlobalVars);

      if (myGasVars->ThermEvolOn == 1)
	myGasVars->temperature = chimes_max(((ChimesFloat) NV_Ith_S(y, data.network_size)) / (1.5f * calculate_total_number_density(myGasVars->abundances, myGasVars->nH_tot, myGlobalVars) * BOLTZMANNCGS), myGasVars->TempFloor); 

      SUNLinSolFree(LS_sun);
      SUNMatDestroy(A_sun);
      N_VDestroy_Serial(y);
      N_VDestroy_Serial(abstol_vector);
      CVodeFree(&cvode_mem);

      free_current_rates_memory(&chimes_current_rates, myGlobalVars); 

      return;
  }
}  


#endif
