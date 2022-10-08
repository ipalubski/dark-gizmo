/* here is where we call the core of the SIDM calculation for DM particle-particle interactions */
#ifdef DM_SIDM
{
	double Pj_dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j);
	//if( ((1 << local.Type) & (DM_SIDM)) && ((1 << P[j].Type) & (DM_SIDM)) && (local.ID != P[j].ID) && (local.dtime <= Pj_dtime))
    //{
		double h_si = kernel.h_i, m_si = 0.5*(local.Mass + P[j].Mass); //0.5*(kernel.h_i + kernel.h_j)
#ifndef DM_SIDM_AREPO
	/* check if target+neighbor are an SIDM candidate, and against self-interaction */
		if((local.dtime==Pj_dtime) && (local.ID > P[j].ID)) continue; // ensures interaction will only be calculated once for each pair //
#ifdef GRAIN_COLLISIONS
	  	double prob = prob_of_grain_interaction(local.Grain_CrossSection_PerUnitMass , local.Mass, kernel.r, h_si, kernel.dv, local.dtime, j);
#endif
#ifdef DM_SIDM_GIZMO
	  	double prob = prob_of_interaction(m_si, kernel.r, h_si, kernel.dv, local.dtime);
#endif
#ifdef DM_SIDM_ROBERTSON
		double prob = prob_of_interaction_robertson(m_si, kernel.r, h_si, kernel.dv, local.dtime);
		//double prob1= prob_of_interaction(m_si, kernel.r, h_si, kernel.dv, local.dtime);
		//printf("P/P1 = %f\n",prob/prob1);
		//if(kernel.r > h_si){printf("h = %f, r/h = %f\n",h_si, kernel.r/h_si);}
#endif
    	if(prob > 0.002) {out.dtime_sidm = DMIN(out.dtime_sidm , local.dtime*(0.002/prob));} // timestep condition not being met as desired, warn code to lower timestep next turn //
#endif

		int scatter = 0;
#ifdef DM_SIDM_AREPO
		if(local.R < local.SIprob){
			scatter = 1;
		}
#else
		if( ((1 << local.Type) & (DM_SIDM)) && ((1 << P[j].Type) & (DM_SIDM)) && (local.ID != P[j].ID) && (local.dtime <= Pj_dtime))
    	{
			if(gsl_rng_uniform(random_generator) < prob){
				scatter = 1;
			}
		}
#endif

		if(scatter == 1) //If we scatter with this particles give them both a kick
		{
#ifdef WAKEUP
			if(!(TimeBinActive[P[j].TimeBin])) {if(WAKEUP*local.dtime < Pj_dtime) {
				#pragma omp atomic write
				PPPZ[j].wakeup=1;
				#pragma omp atomic write
				NeedToWakeupParticles_local = 1;
			}}
#endif
			double kick[3]; calculate_interact_kick(kernel.dv, kick, m_si);
			int k; for(k=0;k<3;k++) {
			out.sidm_kick[k] -= (P[j].Mass/m_si)*kick[k];
			#pragma omp atomic
			P[j].Vel[k] += (local.Mass/m_si)*kick[k]; // this variable is modified here so need to do this carefully here to ensure we don't multiply-write at the same time
			}
			out.si_count++;
			//double ri = sqrt(local.Pos[0]*local.Pos[0] + local.Pos[1]*local.Pos[1] + local.Pos[2]*local.Pos[2]);
			//log_self_interactions(ri);
			#pragma omp atomic
			P[j].NInteractions++;
		}
    //}
}
#endif