#ifdef DM_SIDM_AREPO

extern int Ndm_active;

extern struct SIDM_temp_particle_data{
    int ngbcount;
    MyIDType SItarget;
    double SIprob;
    double R;
    MyIDType ngblist_sum[350];
    double ngbprob[350];
#ifdef DM_NGB_SORT
    double ngbr[350];
#endif

}
*SIDMtempInfo;

void sidm_start(void);
void sidm_end(void);
#endif