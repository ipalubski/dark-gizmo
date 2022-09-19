#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        HYDRO_MESHLESS_FINITE_MASS\n"
"        ADAPTIVE_GRAVSOFT_FORALL=2\n"
"        MULTIPLEDOMAINS=4\n"
"        DM_SIDM=2\n"
"        DM_SIDM_AREPO\n"
"        DM_NGB_SORT\n"
"        HAVE_HDF5\n"
"        USE_FFTW3\n"
"\n");
}
