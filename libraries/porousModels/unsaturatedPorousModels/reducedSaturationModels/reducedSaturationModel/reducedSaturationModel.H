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

Namespace
    Foam::reducedSaturationModel

Description
    A namespace for various capillarity Effect implementations.

Class
    Foam::reducedSaturationModel

Description
    An abstract base class for computing the reduced saturation.

SourceFiles
    reducedSaturationModel.C
    reducedSaturationModelNew.C

Authors
    06/02/2020: C.Soulaine

\*---------------------------------------------------------------------------*/

#ifndef reducedSaturationModel_H
#define reducedSaturationModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class reducedSaturationModel Declaration
\*---------------------------------------------------------------------------*/

class reducedSaturationModel
{

protected:

    // Protected data
    dictionary reducedSaturationProperties_;

    // reference to the saturation field
    const volScalarField& Sb_;

    // Effective saturation field
    volScalarField Se_;

    // Derivative of the effective saturation
    volScalarField dSedS_;

    // Private Member Functions

    //- Disallow copy construct
    reducedSaturationModel(const reducedSaturationModel&);

    //- Disallow default bitwise assignment
    void operator=(const reducedSaturationModel&);

public:

    //- Runtime type information
    TypeName("reducedSaturationModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        reducedSaturationModel,
        dictionary,
        (
            const dictionary& reducedSaturationProperties,
            const volScalarField& Sb
        ),
        (reducedSaturationProperties, Sb)
    );


    // Selectors

    //- Return a reference to the selected capillarity model
    static autoPtr<reducedSaturationModel> New
    (
        const dictionary& reducedSaturationProperties,
        const volScalarField& Sb
    );


    // Constructors

    //- Construct from components
    reducedSaturationModel
    (
        const dictionary& reducedSaturationProperties,
        const volScalarField& Sb
    );


    //- Destructor
    virtual ~reducedSaturationModel() {}


    // Member Functions

    //- Return effective saturation
    const volScalarField & Se() const
    {
        return Se_;
    }

    //- Return the derivative of Se
    const volScalarField & dSedS() const
    {
        return dSedS_;
    }

    // update effective saturation
    virtual void updateReducedSaturation() = 0;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
