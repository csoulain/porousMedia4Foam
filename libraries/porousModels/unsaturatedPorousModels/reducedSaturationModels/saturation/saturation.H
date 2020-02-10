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
    Foam::reducedSaturationModels::saturation

Description
    Calculation of the effective saturation to compute kr and pc using
    Se = (S-Smin)/(Smax-Smin)

SourceFiles
    saturation.C

\*---------------------------------------------------------------------------*/

#ifndef saturation_H
#define saturation_H

#include "reducedSaturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reducedSaturationModels
{

/*---------------------------------------------------------------------------*\
    Class saturation Declaration
    \*---------------------------------------------------------------------------*/

class saturation
:
    public reducedSaturationModel
{

    //- Private data


public:

    //- Runtime type information
    TypeName("saturation");

    // Constructors

    //- Construct from components
    saturation
    (
        const dictionary& dict,
        const volScalarField& Sb
    );

    //- Destructor
    ~saturation()
    {}

    // Member Functions

    //- Compute the reduced saturation
    virtual void updateReducedSaturation()
    {
        Se_=Sb_;
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reducedSaturationModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
