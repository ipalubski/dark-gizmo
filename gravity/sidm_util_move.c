#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../../allvars.h"
#include "../../proto.h"

/* This file containes some functions for SIDM setup */

/* This function initializes a temporary structure to hold neighbor particle data for the SIDM routines*/

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
    
    if(Ndm_active > 0)
    {
        SIDMtempInfo = (struct SIDM_temp_particle_data *) mymalloc("SIDMtempInfo", Ndm_active * sizeof(struct SIDM_temp_particle_data));
    } else {
        SIDMtempInfo = (struct SIDM_temp_particle_data *) mymalloc("SIDMtempInfo", 1 * sizeof(struct SIDM_temp_particle_data));
    }
    memset( &SIDMtempInfo[0], 0, Ndm_active * sizeof(struct SIDM_temp_particle_data));

}

void sidm_end(void)
{
    myfree(SIDMtempInfo);
}