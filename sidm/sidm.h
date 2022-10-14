#ifdef DM_SIDM_AREPO

extern int Ndm_active;

extern struct SIDM_temp_particle_data{
    int ngbcount;
    MyIDType SItarget;
    double SIprob;
    double R;
    MyIDType ngblist_sum[350];
    double ngbprob[350];
    double ngbr[350];

}
*SIDMtempInfo;

void sidm_start(void);
void sidm_end(void);
#endif