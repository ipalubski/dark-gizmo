This is a patched version of the publicly-available GIZMO (http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html). DARK-GIZMO reimplemented SIDM scattering with a smoothed-particle hydrodynamics (SPH) method. It supports isotropic and anisotropic (Yukawa-like scattering potential) and isotropic resonance SIDM cross sections. Additionally, an option to differentiate particle force softening and particle smoothing is added as well as an evolving 3-component Baryon Disk Potential is incorporated for DM only sims. \n
\n
New public compilation-time flags for the SIDM scattering are:
\n
DM_SIDM_AREPO # Use the SPH scattering method used in Arepo \n
DM_SIDM_GIZMO # Use the original GIZMO scattering method \n
SIDM_ISOTROPIC # Set scattering to be isotropic \n
SIDM_ANISOTROPIC # Set scattering to be anisotropic according to the cross section parameters. Only works for the Yukawa scattering potential. \n
PMAXLOW / PMAXHIGH $ Set the SIDM timestep accuracy parameter kappa to 0.002 / 0.02 (Default) This may improve simulation accuracy for large scattering scross sections. \n

These need to be set in addition to DM_SIDM. \n

For Use, Authorship & Citation Requirements see section 3 of the GIZMO documentation: http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html \n
