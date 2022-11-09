/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::twoPhaseMixture

Description
    densityModel module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef densityModel_H
#define densityModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class densityModel
{

  protected:

    // Protected data

    const fvMesh& mesh_;
//    const volScalarField &eps_;

    volScalarField rho_;

  // Private Member Functions

    //- Disallow copy construct
    densityModel(const densityModel&);


public:

  //- Runtime type information
    TypeName("densityModel");


    // Declare run-time constructor selection table

      declareRunTimeSelectionTable
      (
          autoPtr,
          densityModel,
          dictionary,
          (
              const fvMesh& mesh,
              const dictionary& dict
          ),
          (mesh, dict)
      );


  // Selectors

      //- Return a reference to the selected viscosity model
      static autoPtr<densityModel> New
      (
          const fvMesh& mesh,
          const dictionary& dict
      );

  // Constructors

    //- Construct from components

    densityModel
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

//    densityModel();

    //- Destructor
    virtual ~densityModel()
    {}


    // Member Functions


      virtual tmp<volScalarField> rho() const =0;


      virtual void updateDensity() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
