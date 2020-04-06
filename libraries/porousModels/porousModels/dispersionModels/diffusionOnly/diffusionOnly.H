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
    A dispersion model that considers diffusion only. There are no tortuosity
    nor hydrodynamic dispersion effects.

    diffusionOnly model:
    Deff = Di

SourceFiles
    diffusionOnly.C

Authors
    C. Soulaine 6/4/2020 : creation

\*---------------------------------------------------------------------------*/

#ifndef diffusionOnly_H
#define diffusionOnly_H

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

class diffusionOnly : public dispersionModel
{
  protected:

    dictionary diffusionOnlyDict_;

    // Protected data


    dimensionedScalar Di_;              //molecular diffusion coefficient - m2/s
    volScalarField Deff_;                //effectiveDispersion


  // Private Member Functions

    //- Disallow copy construct
    diffusionOnly(const diffusionOnly&);


public:

  //- Runtime type information
    TypeName("diffusionOnly");

  // Declare run-time constructor selection table

  // Selectors

    //- Construct from components
    diffusionOnly
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~diffusionOnly()
    {}


    // Member Functions
    virtual tmp<volScalarField> effectiveDispersion() const;

    virtual void updateDispersion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace diffusionOnly
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
