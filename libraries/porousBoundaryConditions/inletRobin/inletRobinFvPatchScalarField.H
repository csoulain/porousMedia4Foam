/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::inletRobinFvPatchScalarField

Description
    This boundary condition provides a free-stream condition for pressure.

    It is an outlet-inlet condition that uses the velocity orientation to
    continuously blend between zero gradient for normal inlet and fixed value
    for normal outlet flow.

Usage
    \table
        Property        | Description                | Required | Default value
        U               | velocity field name        | no       | U
        kineticValue | freestream pressure        | yes      |
        supersonic      | Switch for supersonic flow | no       | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            inletRobin;
        kineticValue uniform 1e5;
    }
    \endverbatim

    Note:
        This condition is designed to operate with a freestreamVelocity
        condition

See also
    Foam::mixedFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    inletRobinFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef inletRobinFvPatchScalarField_H
#define inletRobinFvPatchScalarField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class inletRobinFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class inletRobinFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of the velocity field
        word diffusivityName_;

        word phiName_;


        word saturationName_;

        word porosityName_;


        scalar velocityValue_;

        scalar diffusivityValue_;

        scalar inletValue_;


public:

    //- Runtime type information
    TypeName("inletRobin");


    // Constructors

        //- Construct from patch, internal field and dictionary
        inletRobinFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given inletRobinFvPatchScalarField onto
        //  a new patch
        inletRobinFvPatchScalarField
        (
            const inletRobinFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        inletRobinFvPatchScalarField
        (
            const inletRobinFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        inletRobinFvPatchScalarField
        (
            const inletRobinFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new inletRobinFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

/*
            const scalarField& kineticValue() const
            {
                return refValue();
            }

            scalarField& kineticValue()
            {
                return refValue();
            }
*/

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
