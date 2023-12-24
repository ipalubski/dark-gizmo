This is a patch to the current publicly-available version of GIZMO (http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html). It adds isotropic and anisotropic (Yukawa-like scattering potential) and isotropic resonance SIDM cross sections as well as an evolving 3-component Baryon Disk Potential.

New public compilation-time flags for the SIDM scattering are:

DM_SIDM_AREPO # Use the SPH scattering method used in Arepo
DM_SIDM_GIZMO # Use the original GIZMO scattering method
SIDM_ISOTROPIC # Set scattering to be isotropic
SIDM_ANISOTROPIC # Set scattering to be anisotropic according to the cross section parameters. Only works for the Yukawa scattering potential.
PMAXLOW / PMAXHIGH $ Set the SIDM timestep accuracy parameter kappa to 0.002 / 0.02 (Default) This may improve simulation accuracy for large scattering scross sections.

These need to be set in addition to DM_SIDM.

Installation: 
Simply copy and replace these files in the main GIZMO directory.

For Use, Authorship & Citation Requirements see section 3 of the GIZMO documentation: http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html
