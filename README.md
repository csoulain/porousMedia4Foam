# porousMedia4Foam

**************************************************************************
************* porousMedia4Foam (PM4F) for OpenFOAM ********************
**************************************************************************

* General Informations

Generic library to solve flow and transport in porous media using OpenFOAM.
The package includes:

    - a porous media library for
        * single phase flow
        * multiphase flow (relative permeability and capillary pressure)
        * reactive transport

    - solvers created using the porous media library
        * darcyFoam : saturated flow in porous media (Darcy's law)
        * impesFoam : unsaturated flow with Darcy's law for two-phase
        * dbsFoam   : hybrid-scale modeling for single-phase flow using
                      Darcy-Brinkman-Stokes equation

    - tutorials

* History

Some of the codes and tutorial have been adapted from the porousMultiphaseFoam
toolbox developed and maintained by Pierre Horgue.

* Install

- This toolbox has been tested on OpenFOAM v11 only

- It only needs a standard OpenFOAM installation from www.openfoam.org


- Read the COPYING_OPENFOAM file for information about OpenFOAM and this
  toolbox Copyrights.

* Installation instructions :

- First, source the OpenFOAM configuration file, i.e. (example for ubuntu
  version) :

  	  source /opt/openfoamv7/etc/bashrc

- then in the "porousMedia4Foam" directory, run :

       	  ./Allwmake

  to install the package.

- Dynamic libraries are compiled in the standard OpenFOAM user directory :

  $FOAM_USER_LIBBIN

- The executable solver "impesFoam" is placed in the standard OpenFOAM user
  directory $FOAM_USER_APPBIN.

- Each tutorial directory contains "run" and "clean" files to test installation
  and validate the solver.

- A python script runTutorials.py can be used to test all components.

- To remove compilation and temporary files, run :

	./Allwclean

- see the ReleaseNotes.txt file for detailed information about the toolbox.

********************************************************************************
