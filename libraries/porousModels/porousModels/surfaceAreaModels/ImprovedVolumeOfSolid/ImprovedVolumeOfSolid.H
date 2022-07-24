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
    Foam::surfaceAreaModels::improvedVolumeOfSolid

Description
    Compute the surface area from the solid phase indicator function.
    Used for solving reactive transport problems at the pore-scale

    Ae = mag(grad(Ys))

Reference
    Soulaine et al. (2017) Mineral dissolution and wormholing from a pore-scale
    perspective. Journal of Fluid Mechanics 827, 457-483

SourceFiles
    improvedVolumeOfSolid.C

Authors
    17/02/2020 : C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef improvedVolumeOfSolid_H
#define improvedVolumeOfSolid_H

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

class improvedVolumeOfSolid : public surfaceAreaModel
{


  protected:

    // Protected data

    const volScalarField & Ys_;

    volScalarField Ae_;


  // Private Member Functions

    //- Disallow copy construct
    improvedVolumeOfSolid(const improvedVolumeOfSolid&);


public:

  //- Runtime type information
    TypeName("improvedVolumeOfSolid");


  // Declare run-time constructor selection table


  // Selectors

    //- Construct from components
    improvedVolumeOfSolid
    (
      const fvMesh& mesh,
      const volScalarField& Ys,
      const dictionary& dict
    );


    //- Destructor
    virtual ~improvedVolumeOfSolid()
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
