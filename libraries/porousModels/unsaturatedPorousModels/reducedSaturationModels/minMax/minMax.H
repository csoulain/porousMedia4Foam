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
    Foam::reducedSaturationModels::minMax

Description
    Calculation of the effective saturation to compute kr and pc using
    Se = (S-Smin)/(Smax-Smin)

SourceFiles
    minMax.C

\*---------------------------------------------------------------------------*/

#ifndef minMax_H
#define minMax_H

#include "reducedSaturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reducedSaturationModels
{

/*---------------------------------------------------------------------------*\
    Class minMax Declaration
    \*---------------------------------------------------------------------------*/

class minMax
:
    public reducedSaturationModel
{

    //- Private data
    dictionary minMaxCoeffs_;


    volScalarField Smin_;
    volScalarField Smax_;


public:

    //- Runtime type information
    TypeName("minMax");

    // Constructors

    //- Construct from components
    minMax
    (
        const dictionary& dict,
        const volScalarField& Sb
    );

    //- Destructor
    ~minMax()
    {}

    // Member Functions

    //- Compute the reduced saturation
    virtual void updateReducedSaturation()
    {
        Se_=(Sb_+SMALL-Smin_)/(Smax_-Smin_);
        dSedS_ = 1./(Smax_-Smin_);
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reducedSaturationModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
