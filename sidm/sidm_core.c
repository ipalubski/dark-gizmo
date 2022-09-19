#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include <unistd.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#define GSLWORKSIZE 100000

/*! \file sidm_routines.c
 *  \brief Fuctions and routines needed for the calculations of dark matter self interactions
 *
 *  This file contains the functions and routines necesary for the computation of
 *  the self-interaction probabilities and the velocity kicks due to the interactios.
 *  Originally written by Miguel Rocha, rocham@uci.edu. Oct 2010. Updated on 2014 & re-written by PFH March 2018
 */

/*! This function calculates the interaction probability between two particles.
 *  It checks if comoving integration is on and does the necesary change of
 *  variables and units.
 */

#ifdef DM_SIDM


double prob_of_interaction(double mass, double r, double h_si, double dV[3], double dt, double rholoc)
{
    double dVmag = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]) / All.cf_atime; // velocity in physical
    double rho_eff = mass / (h_si*h_si*h_si) * All.cf_a3inv; // density in physical
    //double rho_eff = rholoc;
    
    //    if(rholoc > 1.01e+8){
    //    printf("sidm density is %.15f\n",rholoc);
    double cx_eff = All.DM_InteractionCrossSection * g_geo(r/h_si); // effective cross section (physical) scaled to cgs
    double units = UNIT_SURFDEN_IN_CGS; // needed to convert everything to cgs
    if(All.DM_InteractionVelocityScale>0) {double x=dVmag/All.DM_InteractionVelocityScale; cx_eff/=1+x*x*x*x;} // take velocity dependence
    return rho_eff * cx_eff * dVmag * dt * units; // dimensionless probability
}

/*! This routine sets the kicks for each particle after it has been decided that they will
 *  interact. It uses an algorithm tha conserves energy and momentum but picks a random direction so it does not conserves angular momentum. */
#if !defined(GRAIN_COLLISIONS) /* if using the 'grain collisions' module, these functions will be defined elsewhere [in the grains subroutines] */
void calculate_interact_kick(double dV[3], double kick[3], double m)
{
    double dVmag = (1-All.DM_DissipationFactor)*sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]);
    if(dVmag<0) {dVmag=0;}
    if(All.DM_KickPerCollision>0) {double v0=All.DM_KickPerCollision; dVmag=sqrt(dVmag*dVmag+v0*v0);}
    double cos_theta = 2.0*gsl_rng_uniform(random_generator)-1.0, sin_theta = sqrt(1.-cos_theta*cos_theta), phi = gsl_rng_uniform(random_generator)*2.0*M_PI;
    kick[0] = 0.5*(dV[0] + dVmag*sin_theta*cos(phi));
    kick[1] = 0.5*(dV[1] + dVmag*sin_theta*sin(phi));
    kick[2] = 0.5*(dV[2] + dVmag*cos_theta);
}
#endif

#ifdef DM_SIDM_AREPO
// This function initializes a temporary structure to hold neighbor particle data for the SIDM routines
void sidm_start(void)
{
    int i, Ndm;
    Ndm_active = 0;

    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if((1 << P[i].Type) & (DM_SIDM)){
            P[i].IndexMapToTempsidmStruc = Ndm_active;
            Ndm_active++;
        }
    }
    //printf("Ndm_active = %i\n", Ndm_active);
    if(Ndm_active > 0)
    {
        SIDMtempInfo = (struct SIDM_temp_particle_data *) mymalloc("SIDMTempInfo", Ndm_active * sizeof(struct SIDM_temp_particle_data));
    } else {
        SIDMtempInfo = (struct SIDM_temp_particle_data *) mymalloc("SIDMTempInfo", 1 * sizeof(struct SIDM_temp_particle_data));
    }

    memset(&SIDMtempInfo[0], 0, Ndm_active * sizeof(struct SIDM_temp_particle_data));
}
/* Free temporary structure and other final SIDM operations*/
void sidm_end(void)
{
    myfree(SIDMtempInfo);
}

#endif

/*! This function returns the value of the geometrical factor needed for the calculation of the interaction probability. */
double g_geo(double r)
{
    double f, u; int i; u = r / 2.0 * GEOFACTOR_TABLE_LENGTH; i = (int) u;
    if(i >= GEOFACTOR_TABLE_LENGTH) {i = GEOFACTOR_TABLE_LENGTH - 1;}
    if(i <= 1) {f = 0.992318  + (GeoFactorTable[0] - 0.992318)*u;} else {f = GeoFactorTable[i - 1] + (GeoFactorTable[i] - GeoFactorTable[i - 1]) * (u - i);}
    return f;
}

/* This routine returns sigma(v_rel) form table */
double sigma(double dVmag)
{
  double f, vi; int i;
  vi = dVmag*SIGFACTOR_TABLE_LENGTH/500;
  i = (int) vi;
  if (vi > SIGFACTOR_TABLE_LENGTH)
    {f = 0;}
  else
    {f = SigFactorTable[i-1] + (SigFactorTable[i] - SigFactorTable[i - 1]) * (vi - i);}
  return f;
}

/* This routine integrates the differential cross-section for velocities ranging from 0 to 500 km/s */
void init_sigma_table(void){
  int i;
  for(i = 0; i < SIGFACTOR_TABLE_LENGTH; i++)
    {
      double params[2];
      params[0] = (double)i/SIGFACTOR_TABLE_LENGTH*500;
      params[1] = All.DM_InteractionVelocityScale;
      gsl_function F; gsl_integration_workspace *workspace; workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);
      double result, abserr;
      F.function = &prob_angle_integ;
      F.params = &params;
      gsl_integration_qag(&F, 0.0, M_PI, 0, 1.0e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41,workspace, &result, &abserr);
      SigFactorTable[i] = result;
      gsl_integration_workspace_free(workspace);
    }
}

/*Define the differential cross-section below*/
double prob_angle_integ (double x, void * params) {
  double v = *(double *) params;
  double w = *(double *) (params + sizeof(double));
  double v2 = v*v;
  double w2 = w*w;
  double f = sin(x)/2/((1+v2/w2*sin(x/2)*sin(x/2))*(1+v2/w2*sin(x/2)*sin(x/2)));
  return f;
}

/*! This routine initializes the table that will be used to get the geometrical factor
 *  as a function of the two particle separations. It populates a table with the results of the numerical integration */
void init_geofactor_table(void)
{
    int i; double result, abserr,r;
    gsl_function F; gsl_integration_workspace *workspace; workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);
    for(i = 0; i < GEOFACTOR_TABLE_LENGTH; i++)
    {
        r =  2.0/GEOFACTOR_TABLE_LENGTH * (i + 1);
        F.function = &geofactor_integ;
        F.params = &r;
        gsl_integration_qag(&F, 0.0, 1.0, 0, 1.0e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41,workspace, &result, &abserr);
        GeoFactorTable[i] = 2*M_PI*result;
    }
    gsl_integration_workspace_free(workspace);
}

/*! This function returns the integrand of the numerical integration done on init_geofactor_table(). */
double geofactor_integ(double x, void * params)
{
    double result, abserr, r, newparams[2];
    r = *(double *) params; newparams[0] = r; newparams[1] = x;
    gsl_function F; gsl_integration_workspace *workspace; workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);
    F.function = &geofactor_angle_integ; F.params = newparams;
    
    gsl_integration_qag(&F, -1.0, 1.0, 0, 1.0e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41,workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    /*! This function returns the value W(x). The values of the density kernel as a funtion of x=r/h */
    double wk=0; if(x<1) kernel_main(x, 1, 1, &wk, &wk, -1);
    return x*x*wk*result;
}

/*! This function returns the integrand of the angular part of the integral done on init_geofactor_table(). */
double geofactor_angle_integ(double u, void * params)
{
    double x,r,f;
    r = *(double *) params;
    x = *(double *) (params + sizeof(double));
    f = sqrt(x*x + r*r + 2*x*r*u);
    double wk=0; if(f<1) kernel_main(f, 1, 1, &wk, &wk, -1); /*! This function returns the value W(x). The values of the density kernel as a funtion of x=r/h */
    return wk;
}

/*! This function simply initializes some variables to prevent memory errors */
void init_self_interactions() {int i; for(i = 0; i < NumPart; i++) {P[i].dtime_sidm = 0; P[i].NInteractions = 0;}}

void log_self_interactions(double r)
{
  char name[50];
  //  pid_t pid = getpid();
  sprintf(name,"Ndmsi%d.txt",getpid());
  FILE *f = fopen(name, "a");
  if (f){
    fprintf(f, "%.10f \t %.10g \t %.10g\n", r, All.Time, All.Time + All.TimeStep);
    //printf("\nNumber of DM self-interactions at this time-step: %lu\n", All.Ndmsi);
  }
  fclose(f);
}

#endif
