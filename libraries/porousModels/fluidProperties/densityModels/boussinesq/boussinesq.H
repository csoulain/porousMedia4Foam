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

#ifndef boussinesq_H
#define boussinesq_H

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

class boussinesq : public densityModel
{

  protected:

    // Protected data

    dictionary boussinesqDict_;

    dimensionedScalar rho0_;

    dimensionedScalar beta_;

    dimensionedScalar T0_;

    word TName_;

    volScalarField rho_;

  // Private Member Functions

    //- Disallow copy construct
    boussinesq(const boussinesq&);


public:

  //- Runtime type information
    TypeName("boussinesq");


  // Declare run-time constructor selection table


  // Selectors

    //- Construct from components
    boussinesq
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~boussinesq()
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
