/* here is where we call the core of the SIDM calculation for DM particle-particle interactions */
#ifdef DM_SIDM
{
  double Pj_dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j); double prob = 0;
  double m_si = 0.5*(local.Mass + P[j].Mass);
  //if((local.dtime==Pj_dtime) && (local.ID > P[j].ID)) continue;
#ifndef DM_SIDM_AREPO
  /* check if target+neighbor are an SIDM candidate, and against self-interaction */
  if((local.dtime==Pj_dtime) && (local.ID > P[j].ID)) continue; // ensures interaction will only be calculated once for each pair //
#ifdef GRAIN_COLLISIONS
  double h_si = 0.5*(kernel.h_i + kernel.h_j);//, m_si = 0.5*(local.Mass + P[j].Mass);
  prob = prob_of_grain_interaction(local.Grain_CrossSection_PerUnitMass , local.Mass, kernel.r, h_si, kernel.dv, local.dtime, j);
#endif
#ifdef DM_SIDM_GIZMO
  double h_si = 0.5*(kernel.h_i + kernel.h_j);//, m_si = 0.5*(local.Mass + P[j].Mass);
  prob = prob_of_interaction(m_si, kernel.r, h_si, kernel.dv, local.dtime);
#endif
#ifdef DM_SIDM_ROBERTSON
  double h_si = kernel.h_i;//, m_si = 0.5*(local.Mass + P[j].Mass);
  prob = prob_of_interaction_robertson(m_si, kernel.r, h_si, kernel.dv, local.dtime);
  //if(prob > 0.02){printf("P = %f, too high, lower timestep!\n",prob);}
#endif
    	if(prob > All.Ptarget) {out.dtime_sidm = DMIN(out.dtime_sidm , local.dtime*(All.Ptarget/prob));} // timestep condition not being met as desired, warn code to lower timestep next turn //
#endif

		int scatter = 0;
#ifdef DM_SIDM_AREPO
		if(local.R < local.SIprob){
		  scatter = 1;
		}
#else
		if( ((1 << local.Type) & (DM_SIDM)) && ((1 << P[j].Type) & (DM_SIDM)) && (local.ID != P[j].ID) && (local.dtime <= Pj_dtime))
		{
		  //if((local.dtime==Pj_dtime) && (local.ID > P[j].ID)) continue; 
		    if(gsl_rng_uniform(random_generator) < prob){
		      //if(prob > 0.02){printf("P = %f, too high, lower timestep!\n",prob);}
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
	                #pragma omp atomic
			P[j].NInteractions++;

			#ifdef NSI
			//double ri = sqrt(pow(local.Pos[0],2) + pow(local.Pos[1],2) + pow(local.Pos[2],2));
			//double rj = sqrt(pow(P[j].Pos[0],2) + pow(P[j].Pos[1],2) + pow(P[j].Pos[2],2));
			double dx = sqrt(pow(local.Pos[0]-P[j].Pos[0],2) + pow(local.Pos[1]-P[j].Pos[1],2) + pow(local.Pos[2]-P[j].Pos[2],2));
			log_self_interactions(dx,local.Pos[0],local.Pos[1],local.Pos[2]);
                        #endif
		}
//}
}
//#endif
#endif
