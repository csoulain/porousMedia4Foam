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
    PhreeqcRM module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef fromGeoChemicalPackage_H
#define fromGeoChemicalPackage_H

#include "volFields.H"
#include "speciesTable.H"
//#include "basicMultiComponentMixture.H"
//#include "runTimeSelectionTables.H"

#include "densityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace densityModels
{
/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class fromGeoChemicalPackage : public densityModel
{

  protected:

    // Protected data

    dictionary fromGeoChemicalPackageDict_;

    dimensionedScalar rho0_;

    volScalarField rho_;

  // Private Member Functions

    //- Disallow copy construct
    fromGeoChemicalPackage(const fromGeoChemicalPackage&);


public:

  //- Runtime type information
    TypeName("fromGeoChemicalPackage");


  // Declare run-time constructor selection table


  // Selectors

    //- Construct from components
    fromGeoChemicalPackage
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~fromGeoChemicalPackage()
    {}


    // Member Functions
    virtual tmp<volScalarField> rho() const;

    virtual void updateDensity();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace densityModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
