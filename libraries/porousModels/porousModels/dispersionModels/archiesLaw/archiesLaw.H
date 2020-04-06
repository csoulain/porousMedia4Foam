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
    Archies model for tortuosity. Feed to update molecular diffusion
    and consequently effectiveDispersion.

    General dispersion model:
    Deff = porosity^{n}*molecularDiffusionCoeff

SourceFiles
    archiesLaw.C

Authors
    - Saideep Pavuluri

\*---------------------------------------------------------------------------*/

#ifndef archiesLaw_H
#define archiesLaw_H

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

class archiesLaw : public dispersionModel
{
  protected:

    dictionary archiesLawDict_;

    // Protected data

    // Name of the porosity field
    word epsName_;


    dimensionedScalar Di_;              //molecular diffusion coefficient - m2/s
    volScalarField Deff_;                //effectiveDispersion
    const volScalarField & eps_;        //porosity
    scalar n_;                          //for porosity

  // Private Member Functions

    //- Disallow copy construct
    archiesLaw(const archiesLaw&);


public:

  //- Runtime type information
    TypeName("archiesLaw");

  // Declare run-time constructor selection table

  // Selectors

    //- Construct from components
    archiesLaw
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~archiesLaw()
    {}


    // Member Functions
    virtual tmp<volScalarField> effectiveDispersion() const;

    virtual void updateDispersion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace dispersionModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
