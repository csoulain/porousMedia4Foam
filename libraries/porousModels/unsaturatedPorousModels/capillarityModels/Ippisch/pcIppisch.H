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
    Foam::capillarityModels::Ippisch

Description
    Ippisch-Vogel-Bastian capillary pressure model (AWR 2006)

SourceFiles
    pcIppisch.C

\*---------------------------------------------------------------------------*/

#ifndef pcIppisch_H
#define pcIppisch_H

#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
                    Class pcIppisch Declaration
\*---------------------------------------------------------------------------*/

class pcIppisch
:
    public capillarityModel
{
    //- Ippisch coefficients
    dictionary pcIppischCoeffs_;
    volScalarField m_,n_,alpha_,tau_,he_;

    //- cut-off point saturation (computed with entry pressure)
    volScalarField Sc_;

public:

    //- Runtime type information
    TypeName("Ippisch");

    // Constructors

    //- Construct from components
    pcIppisch
    (
        const word& name,
        const dictionary& dict,
        const volScalarField& Sb,
        const autoPtr<reducedSaturationModel> & reducedSaturationModelPtr
    );

    //- Destructor
    ~pcIppisch()
    {}

    // Member Functions

    //- Correct the capillary pressure
    void correct()
    {
        FatalErrorIn("void correct()")
            << "Not implemented in pcIppish.H " << abort(FatalError);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
