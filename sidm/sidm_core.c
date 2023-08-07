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
//#include "sidm.h"

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

double prob_of_interaction(double mass, double r, double h_si, double dV[3], double dt)
{
    double dVmag = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]) / All.cf_atime; // velocity in physical
    double rho_eff = mass / (h_si*h_si*h_si) * All.cf_a3inv; // density in physical
    double cx_eff = All.DM_InteractionCrossSection * g_geo(r/h_si); // effective cross section (physical) scaled to cgs
    double units = UNIT_SURFDEN_IN_CGS; // needed to convert everything to cgs
    if(All.DM_InteractionVelocityScale>0) {double x=dVmag/All.DM_InteractionVelocityScale; double x2 = x*x; cx_eff*=4*M_PI*(1/(1+x2)-log(1+x2)/(x2*(2+x2)));} // take velocity dependence
    return rho_eff * cx_eff * dVmag * dt * units; // dimensionless probability
}

double prob_of_interaction_robertson(double mass, double r, double h_si, double dV[3], double dt)
{
    double dVmag = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]) / All.cf_atime; // velocity in physical
    double rho_eff = mass/(1.33*M_PI*h_si*h_si*h_si) * All.cf_a3inv; // density in physical
    double cx_eff = All.DM_InteractionCrossSection; // effective cross section (physical) scaled to cgs
    double units = UNIT_SURFDEN_IN_CGS; // needed to convert everything to cgs
    if(All.DM_InteractionVelocityScale>0) {double x=dVmag/All.DM_InteractionVelocityScale; double x2 = x*x; cx_eff*=4*M_PI*(1/(1+x2)-log(1+x2)/(x2*(2+x2)));} // take velocity dependence
    //if((rho_eff * cx_eff * dVmag * dt * units) > 1.0){printf("P = %f\n",rho_eff * cx_eff * dVmag * dt * units);}
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
    #if defined(SIDM_ISOTROPIC)
    double cos_theta = 2.0*gsl_rng_uniform(random_generator)-1.0, sin_theta = sqrt(1.-cos_theta*cos_theta), phi = gsl_rng_uniform(random_generator)*2.0*M_PI;
    #endif
    #ifdef SIDM_ANISOTROPIC
    double R = gsl_rng_uniform(random_generator);
    double vw2 = All.DM_InteractionVelocityScale*All.DM_InteractionVelocityScale/dVmag/dVmag, cos_theta = (R*(1.0+2.0*vw2)-1.0-vw2)/(R-1.0-vw2);
    //double vw2 = dVmag*dVmag/All.DM_InteractionVelocityScale/All.DM_InteractionVelocityScale, cos_theta = (1.0+vw2-R*(vw2+2))/(1.0+vw2-R*vw2);//cos_theta = (R*(1.0+2.0*vw2)-1.0-vw2)/(R-1.0-vw2);
    double sin_theta = sqrt(1.-cos_theta*cos_theta), phi = acos(cos_theta);
    //printf("phi = %f, vw2 = %f, R = %f, cos_theta = %f \n", phi,vw2,R,cos_theta);
    //printf("num = %f, den = %f \n",R*(1+2*vw2)-1.0-vw2,R-1.0-vw2);
    #endif
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

/*int cmp(const void *a, const void *b)
    {
        struct ngb_list_data *a1 = (struct ngb_list_data *)a;
        struct ngb_list_data *a2 = (struct ngb_list_data *)b;
        if ((*a1).r < (*a2).r)
        return -1;
        else if ((*a1).r > (*a2).r)
        return 1;
        else
        return 0;
    }
*/
#endif

/*! This function returns the value of the geometrical factor needed for the calculation of the interaction probability. */
double g_geo(double r)
{
    double f, u; int i; u = r / 2.0 * GEOFACTOR_TABLE_LENGTH; i = (int) u;
    if(i >= GEOFACTOR_TABLE_LENGTH) {i = GEOFACTOR_TABLE_LENGTH - 1;}
    if(i <= 1) {f = 0.992318  + (GeoFactorTable[0] - 0.992318)*u;} else {f = GeoFactorTable[i - 1] + (GeoFactorTable[i] - GeoFactorTable[i - 1]) * (u - i);}
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

#ifdef NSI
void log_self_interactions(double dx, double x, double y, double z)
{
  char name[50];
  sprintf(name,"Ndmsi%d.txt",getpid());
  FILE *f = fopen(name, "a");
  if (f){
    fprintf(f, "%.10f \t %.10g \t %.10f \t %.10f \t %.10f\n",dx , All.Time, x, y, z);
    //printf("\nNumber of DM self-interactions at this time-step: %lu\n", All.Ndmsi);
  }
  fclose(f);
}

#define ARRSZ   2   /* use 8 or more, set to 2 here to force realloc */
#define MAXC 1024

/* simple implementation of strdup - in the event you don't have it */
char *mystrdupnsi (const char *s)
{
    if (!s)     /* validate s not NULL */
        return NULL;

    size_t len = strlen (s);            /* get length */
    char *sdup = malloc (len + 1);      /* allocate length + 1 */

    if (!sdup)          /* validate */
        return NULL;

    return memcpy (sdup, s, len + 1);   /* pointer to copied string */
}

/* simple function to free all data when done */
void freemydatansi (mydata_t *data, size_t n)
{
    for (size_t i = 0; i < n; i++) {    /* free allocated strings */
        free (data[i].col1);
        free (data[i].col2);
        free (data[i].col3);
        free (data[i].col4);
        free (data[i].col5);
    }
    free (data);    /* free structs */
}

void catnsi(void){
    system("cat Nd* > NSI.txt");
    system("rm Ndmsi*");  
}
/*#ifdef FOF
void log_halo_cm(void)
{
  char name[50];
  //  pid_t pid = getpid();                                                                                                                                                                                    
  //sprintf(name,"Halo_Canter.txt");
  FILE *f = fopen("Halo_Center.txt", "a");
  if (f){
    fprintf(f, "%.3f \t %.3f \t %.3f\n", Halo_CM[0],Halo_CM[1],Halo_CM[2]);
  }
  fclose(f);
}
#endif*/

void write_NSI(void) {
    char buf[MAXC];
    size_t arrsz = ARRSZ, line = 0, row = 0;
    mydata_t *data = NULL;
    /* use filename provided as 1st argument (stdin by default) */
    FILE *fp = fopen ("NSI.txt", "r");

    if (!fp) {  /* validate file open for reading */
        printf("file open failed");
    }

    /* allocate an 'arrsz' initial number of struct */
    if (!(data = malloc (arrsz * sizeof *data))) {
        printf("malloc-data");
    }

    while (fgets (buf, MAXC, fp)) {         /* read each line from file */
        char c1[MAXC], c2[MAXC], c3[MAXC], c4[MAXC], c5[MAXC];  /* temp strings for c1,2,5 */
        //int c3, c4;                         /* temp ints for c3,4 */                                                                                                                                                         
        size_t len = strlen (buf);          /* length for validation */

        line++;     /* increment line count */

        /* validate line fit in buffer */
        if (len && buf[len-1] != '\n' && len == MAXC - 1) {
            printf("error: line %zu exceeds MAXC chars.\n", line);
        }

        if (row == arrsz) { /* check if all pointers used */
            void *tmp = realloc (data, arrsz * 2 * sizeof *data);
            if (!tmp) {     /* validate realloc succeeded */
                printf("realloc-data\n");
                break;      /* break, don't exit, data still valid */
            }
            data = tmp;     /* assign realloc'ed block to data */
            arrsz *= 2;     /* update arrsz to reflect new allocation */
        }

        /* parse buf into fields, handle error on invalid format of line */
        if (sscanf (buf, "%1023s %1023s %1023s %1023s %1023s",
                    &c1, &c2, &c3, &c4, &c5) != 5) {
            printf("error: invalid format line %zu\n", line);
            continue;   /* get next line */
        }

        /* allocate copy strings, assign allocated blocks to pointers */
        if (!(data[row].col1 = mystrdupnsi (c1))) { /* validate copy of c1 */
            fprintf (stderr, "error: malloc-c1 line %zu\n", line);
            break;      /* same reason to break not exit */
        }
        if (!(data[row].col2 = mystrdupnsi (c2))) { /* validate copy of c2 */
            fprintf (stderr, "error: malloc-c1 line %zu\n", line);
            break;      /* same reason to break not exit */
        }
        if (!(data[row].col3 = mystrdupnsi (c3))) { /* validate copy of c5 */
            fprintf (stderr, "error: malloc-c1 line %zu\n", line);
            break;      /* same reason to break not exit */
        }
        if (!(data[row].col4 = mystrdupnsi (c4))) { /* validate copy of c5 */
            fprintf (stderr, "error: malloc-c1 line %zu\n", line);
            break;      /* same reason to break not exit */
        }
        if (!(data[row].col5 = mystrdupnsi (c5))) { /* validate copy of c5 */
            fprintf (stderr, "error: malloc-c1 line %zu\n", line);
            break;      /* same reason to break not exit */
        }
        row++;      /* increment number of row pointers used */
    }
    if (fp != stdin)    /* close file if not stdin */
        fclose (fp);

    remove("NSI.txt");
    char *eptr;
    FILE *f = fopen("NSI.txt", "a");
    for (size_t i = 0; i < row; i++){
      double time;
      time = strtod(data[i].col2, &eptr);
      if(time < All.Time){
        fprintf(f, "%.10s \t %.10s \t %.10s \t %.10s \t %.10s\n", data[i].col1, data[i].col2, data[i].col3, data[i].col4, data[i].col5);
      }
    }
    freemydatansi (data, row);
    fclose(f);
}
#endif
#endif
