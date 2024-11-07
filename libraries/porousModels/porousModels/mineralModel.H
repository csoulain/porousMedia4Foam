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
    Foam::mineralModel

Description
    mineralModel module for geometric porous media properties including
      - porosity
      - absolute permeability models
      - surface area models


SourceFiles
    mineralModel.C

Authors
    - 4/2/2020 C. Soulaine

\*---------------------------------------------------------------------------*/

#ifndef mineralModel_H
#define mineralModel_H

#include "volFields.H"
#include "surfaceInterpolate.H"
#include "surfaceAreaModel.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class mineralModel
{


  protected:

    // Protected data

    const fvMesh &mesh_;

    word mineralName_;

    dictionary mineralDict_;


    autoPtr<surfaceAreaModel> surfaceAreaModelPtr_;


/*
    tmp<volScalarField> invK_;
*/
  // Private Member Functions

    //- Disallow copy construct
    mineralModel(const mineralModel&);


  public:

  // Constructors

    //- Construct from components

    mineralModel
    (
        const fvMesh& mesh,
        const word & name,
        const volScalarField& Ys,
        const dictionary& dict
    );

    //- Destructor
    ~mineralModel()
    {}


    // Member Functions

    //- Return the composition of the multi-component mixture
//    virtual basicMultiComponentMixture& composition() = 0;

    //- Return the composition of the multi-component mixture
//    virtual const basicMultiComponentMixture& composition() const = 0;


    tmp<volScalarField> surfaceArea() const
    {
        return surfaceAreaModelPtr_->surfaceArea();
    }

    void updateSurfaceArea()
    {
        surfaceAreaModelPtr_->updateSurfaceArea();
//        return surfaceAreaModelPtr_->updateSurfaceArea();
    }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
