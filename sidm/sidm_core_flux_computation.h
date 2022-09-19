/* here is where we call the core of the SIDM calculation for DM particle-particle interactions */
#ifdef DM_SIDM
{
  //MyIDType SItarget = local.SItarget;
    /* check if target+neighbor are an SIDM candidate, and against self-interaction */
  double Pj_dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j);
    //if( ((1 << local.Type) & (DM_SIDM)) && ((1 << P[j].Type) & (DM_SIDM)) && (local.ID != P[j].ID) && (local.dtime <= Pj_dtime))
    //{
    //if((local.dtime==Pj_dtime) && (local.ID > P[j].ID)) continue; // ensures interaction will only be calculated once for each pair //
    double m_si = 0.5*(local.Mass + P[j].Mass);
    //double rho = 2;//kernel.DMRho;
	/*#ifdef GRAIN_COLLISIONS
	  double prob = prob_of_grain_interaction(local.Grain_CrossSection_PerUnitMass , local.Mass, kernel.r, h_si, kernel.dv, local.dtime, j);
	  #else
	  double prob = prob_of_interaction(m_si, kernel.r, h_si, kernel.dv, local.dtime, rho);
	  #endif*/
        //if(prob > 0.2) {out.dtime_sidm = DMIN(out.dtime_sidm , local.dtime*(0.2/prob));} // timestep condition not being met as desired, warn code to lower timestep next turn //
    //if(local.SIprob > 0.02) {out.dtime_sidm = DMIN(out.dtime_sidm , local.dtime*(0.02/local.SIprob));}
    
    if (local.R < local.SIprob)
    {
	//printf("PROB = %f, ID = %i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n", local.SIprob,local.ID);
	//printf("R = %f\n",local.R);
	if(local.SItarget == local.ID){printf("target is self!!!!\n");}
	//if(sizeof(local.SItarget) != sizeof(MyIDType)){printf("wrong id size\n");}
#ifdef WAKEUP
	if(!(TimeBinActive[P[j].TimeBin])) {if(WAKEUP*local.dtime < Pj_dtime) {
            #pragma omp atomic write
	    PPPZ[j].wakeup=1;
            #pragma omp atomic write
	    NeedToWakeupParticles_local = 1;
	  }}
#endif
	//double kick[3]; calculate_interact_kick(kernel.dv, kick, m_si);
	double kick[3]; calculate_interact_kick(kernel.dv, kick, m_si);
	int k; for(k=0;k<3;k++) {
	  out.sidm_kick[k] -= (P[j].Mass/m_si)*kick[k];
          #pragma omp atomic
	  P[j].Vel[k] += (local.Mass/m_si)*kick[k]; // this variable is modified here so need to do this carefully here to ensure we don't multiply-write at the same time
	}
	out.si_count++;
	double ri = sqrt(local.Pos[0]*local.Pos[0] + local.Pos[1]*local.Pos[1] + local.Pos[2]*local.Pos[2]);
	//double ri = sqrt(P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2]);
	log_self_interactions(ri);
        #pragma omp atomic
	P[j].NInteractions++;
    }
}
#endif
