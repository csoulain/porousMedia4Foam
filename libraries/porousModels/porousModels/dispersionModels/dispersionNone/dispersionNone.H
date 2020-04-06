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

    dispersionNone model:
    Deff = Di

SourceFiles
    dispersionNone.C

Authors
    C. Soulaine 6/4/2020 : creation

\*---------------------------------------------------------------------------*/

#ifndef dispersionNone_H
#define dispersionNone_H

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

class dispersionNone : public dispersionModel
{
  protected:

//    dictionary dispersionNoneDict_;

    // Protected data


    volScalarField Deff_;                //effectiveDispersion


  // Private Member Functions

    //- Disallow copy construct
    dispersionNone(const dispersionNone&);


public:

  //- Runtime type information
    TypeName("none");

  // Declare run-time constructor selection table

  // Selectors

    //- Construct from components
    dispersionNone
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~dispersionNone()
    {}


    // Member Functions
    virtual tmp<volScalarField> effectiveDispersion() const;

    virtual void updateDispersion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace dispersionNone
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
