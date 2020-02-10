/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::capillarityModels::linear

Description
    Standard linear model for capillary pressure

SourceFiles
    pcLinear.C

\*---------------------------------------------------------------------------*/

#ifndef pcLinear_H
#define pcLinear_H

#include "capillarityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
    Class pcLinear Declaration
    \*---------------------------------------------------------------------------*/

class pcLinear
:
    public capillarityModel
{
    // pcVan Genuchten coefficients
    dictionary pcLinearCoeffs_;

    volScalarField pc0_;
    volScalarField pcMax_;

public:

    //- Runtime type information
    TypeName("linear");

    // Constructors

    //- Construct from components
    pcLinear
    (
        const word& name,
        const dictionary& dict,
        const volScalarField& Sb,
        const autoPtr<reducedSaturationModel> & reducedSaturationModelPtr
    );

    //- Destructor
    ~pcLinear()
    {}

    // Member Functions

    //- Correct the capillarity
    void correct()
    {
        pc_ = pc0_ + (1 - Se_) * (pcMax_ - pc0_);
        dpcdS_ =  pc0_ - pcMax_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
