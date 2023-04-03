#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        HYDRO_MESHLESS_FINITE_VOLUME\n"
"        ADAPTIVE_GRAVSOFT_FORALL=2\n"
"        DM_SIDM=2\n"
"        NSI\n"
"        DM_SIDM_AREPO\n"
"        PMAXHIGH\n"
"        HAVE_HDF5\n"
"        OUTPUT_POTENTIAL\n"
"        EVALPOTENTIAL\n"
"        USE_FFTW3\n"
"\n");
}
