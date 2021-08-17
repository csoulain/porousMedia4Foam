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
    Taylor-Aris dispersion model. Note that this model is specific for
    pipe flows. Therefore the porosity is equal to 1.

    Dispersion model:
    Deff = molecularDiffusionCoeff*(1+channelradius^2*|velocity|^2/(48*molecularDiffusionCoeff^2))

SourceFiles
    dispersionTensorNone.C

\*---------------------------------------------------------------------------*/

#ifndef dispersionTensorNone_H
#define dispersionTensorNone_H

#include "volFields.H"
#include "speciesTable.H"

#include "dispersionTensorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionTensorModels
{
/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class dispersionTensorNone : public dispersionTensorModel
{
  protected:

    dictionary dispersionTensorNoneDict_;

    // Protected data

    // Name of the velocity field
    word UName_;

    volTensorField Deff_;                //effectiveDispersion

  // Private Member Functions

    //- Disallow copy construct
    dispersionTensorNone(const dispersionTensorNone&);

public:

  //- Runtime type information
    TypeName("none");

  // Declare run-time constructor selection table

  // Selectors

    //- Construct from components
    dispersionTensorNone
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~dispersionTensorNone()
    {}


    // Member Functions
    virtual tmp<volTensorField> effectiveDispersionTensor() const;

    virtual void updateDispersionTensor();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace dispersionTensorModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
