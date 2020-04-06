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
    A general model to account for effectiveDispersion. This class includes
    both molecular diffusion and hydrodynamic dispersion.

    General dispersion model:
    Deff = porosity*molecularDiffusionCoeff*(1+hydrodynamicDispersionCoeff/molecularDiffusionCoeff*|velocity|)

SourceFiles
    linearDispersion.C

\*---------------------------------------------------------------------------*/

#ifndef linearDispersion_H
#define linearDispersion_H

#include "volFields.H"
#include "speciesTable.H"

#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{
/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class linearDispersion : public dispersionModel
{
  protected:

    dictionary linearDispersionDict_;

    // Protected data

    // Name of the porosity field
    word epsName_;

    // Name of the velocity field
    word UName_;

    dimensionedScalar Di_;              //molecular diffusion coefficient - m2/s
    dimensionedScalar alphaL_;          //hydrodynamic dispersion coeff - m/s
    volScalarField Deff_;                //effectiveDispersion
    const volScalarField & eps_;        //porosity
    const volVectorField & U_;          //velocity
    scalar n_;                          //for porosity


  // Private Member Functions

    //- Disallow copy construct
    linearDispersion(const linearDispersion&);


public:

  //- Runtime type information
    TypeName("linearDispersion");

  // Declare run-time constructor selection table

  // Selectors

    //- Construct from components
    linearDispersion
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~linearDispersion()
    {}


    // Member Functions
    virtual tmp<volScalarField> effectiveDispersion() const;

    virtual void updateDispersion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace linearDispersion
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
