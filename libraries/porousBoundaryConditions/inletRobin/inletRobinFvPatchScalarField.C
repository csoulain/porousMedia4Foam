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

\*---------------------------------------------------------------------------*/

#include "inletRobinFvPatchScalarField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletRobinFvPatchScalarField::
inletRobinFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    diffusivityName_(dict.lookup<word>("diffusivity")),
    phiName_(dict.lookup<word>("phi")),
    saturationName_(dict.lookup<word>("S")),
    porosityName_(dict.lookup<word>("porosity")),
    velocityValue_(dict.lookup<scalar>("velocityValue")),
    diffusivityValue_(dict.lookup<scalar>("diffusivityValue")),
    inletValue_(dict.lookup<scalar>("inletValue"))

{

//    velocityValue() = scalarField("velocityValue", dict, p.size());

    if (dict.found("value"))
    {
        Info << "test if" <<nl <<endl;

        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        Info << "test else" <<nl <<endl;
        //fvPatchScalarField::operator=(Zero);
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchScalarField::operator=(patchInternalField());
    }


    refValue() = inletValue_; //Zero;
    refGrad() = Zero;
    valueFraction() = 0;
}


Foam::inletRobinFvPatchScalarField::
inletRobinFvPatchScalarField
(
    const inletRobinFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    diffusivityName_(psf.diffusivityName_),
    phiName_(psf.phiName_),
    saturationName_(psf.saturationName_),
    porosityName_(psf.porosityName_),
    velocityValue_(psf.velocityValue_),
    diffusivityValue_(psf.diffusivityValue_),
    inletValue_(psf.inletValue_)
{}


Foam::inletRobinFvPatchScalarField::
inletRobinFvPatchScalarField
(
    const inletRobinFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    diffusivityName_(psf.diffusivityName_),
    phiName_(psf.phiName_),
    saturationName_(psf.saturationName_),
    porosityName_(psf.porosityName_),
    velocityValue_(psf.velocityValue_),
    diffusivityValue_(psf.diffusivityValue_),
    inletValue_(psf.inletValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inletRobinFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& porosity=
        patch().lookupPatchField<volScalarField, scalar>(porosityName_);

    const fvPatchField<scalar>& saturation=
        patch().lookupPatchField<volScalarField, scalar>(saturationName_);

    const fvsPatchField<scalar>& phi=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);


    this->refValue() = inletValue_; //Zero;

//    this->valueFraction() 
//        = 1./(1.+diffusivityValue_/(velocityValue_*this->patch().deltaCoeffs()));

    this->valueFraction() 
//        = 1./(1.+(porosity*saturation*diffusivityValue_)/(phi*this->patch().deltaCoeffs()));
        = 1./(1.+(porosity*saturation*diffusivityValue_)/(SMALL+phi*this->patch().deltaCoeffs()));


    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::inletRobinFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", diffusivityName_);
    writeEntry(os, "velocityValue", velocityValue_);
    writeEntry(os, "diffusivityValue", diffusivityValue_);    
    writeEntry(os, "inletValue", inletValue_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inletRobinFvPatchScalarField
    );
}

// ************************************************************************* //
