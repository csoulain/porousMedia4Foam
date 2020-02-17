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
    Foam::surfaceAreaModels::constantSurfaceArea

Description
    constant specific surface area

    Ae = A0

SourceFiles
    constantSurfaceArea.C

Authors
    17/02/2020 : C. Soulaine
\*---------------------------------------------------------------------------*/

#ifndef constantSurfaceArea_H
#define constantSurfaceArea_H

#include "volFields.H"
#include "surfaceAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceAreaModels
{
/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class constantSurfaceArea : public surfaceAreaModel
{

    dictionary constantSurfaceAreaDict_;

    dimensionedScalar A0_;

    volScalarField Ae_;


  protected:

    // Protected data

  // Private Member Functions

    //- Disallow copy construct
    constantSurfaceArea(const constantSurfaceArea&);


public:

  //- Runtime type information
    TypeName("constantSurfaceArea");


  // Declare run-time constructor selection table


  // Selectors

    //- Construct from components
    constantSurfaceArea
    (
      const fvMesh& mesh,
      const volScalarField& eps,
      const dictionary& dict
    );


    //- Destructor
    virtual ~constantSurfaceArea()
    {}


    // Member Functions
    virtual tmp<volScalarField> surfaceArea() const;

    virtual void updateSurfaceArea();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace absolutePermeabilityModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
