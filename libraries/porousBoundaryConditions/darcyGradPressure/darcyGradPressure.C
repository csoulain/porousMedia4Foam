/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "darcyGradPressure.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
  //  curTimeIndex_(-1)
    MfName_("Mf"),
    phiName_("phi"),
    phiGfName_("phiG"),
    phiPcName_("phiPc")
{}


Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict, false),
//    curTimeIndex_(-1)
    MfName_(dict.lookupOrDefault<word>("Mf", "Mf")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    phiGfName_(dict.lookupOrDefault<word>("phiG","phiG")),
    phiPcName_(dict.lookupOrDefault<word>("phiPc","phiPc"))
{
//    if (dict.found("value") && dict.found("gradient"))
//    {
//        fvPatchField<scalar>::operator=(scalarField("value", dict, p.size()));
//        gradient() = scalarField("gradient", dict, p.size());
//    }
//    else
//    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = Zero;
//    }
}


Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
(
    const darcyGradPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
//    curTimeIndex_(-1)
    MfName_(ptf.MfName_),
    phiName_(ptf.phiName_),
    phiGfName_(ptf.phiGfName_),
    phiPcName_(ptf.phiPcName_)
{}


Foam::darcyGradPressureFvPatchScalarField::darcyGradPressureFvPatchScalarField
(
    const darcyGradPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
//    curTimeIndex_(-1)
    MfName_(wbppsf.MfName_),
    phiName_(wbppsf.phiName_),
    phiGfName_(wbppsf.phiGfName_),
    phiPcName_(wbppsf.phiPcName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void Foam::darcyGradPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& snGradp
)
{
    if (updated())
    {
        return;
    }

    curTimeIndex_ = this->db().time().timeIndex();

    gradient() = snGradp;
    fixedGradientFvPatchScalarField::updateCoeffs();
}
*/

void Foam::darcyGradPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& Mf=
        patch().lookupPatchField<surfaceScalarField, scalar>(MfName_);

    const fvsPatchField<scalar>& phi=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvsPatchField<scalar>& phiGf=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiGfName_);

    const fvsPatchField<scalar>& phiPc=
        patch().lookupPatchField<surfaceScalarField, scalar>(phiPcName_);


    gradient() = - (phi-phiGf-phiPc)/Mf/(patch().magSf());

    fixedGradientFvPatchScalarField::updateCoeffs();

/*

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        FatalErrorInFunction
            << "updateCoeffs(const scalarField& snGradp) MUST be called before"
               " updateCoeffs() or evaluate() to set the boundary gradient."
            << exit(FatalError);
    }

    */
}


void Foam::darcyGradPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "Mf", "Mf", MfName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "phiG", "phiG", phiGfName_);
    writeEntryIfDifferent<word>(os, "phiPc", "phiPc", phiPcName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        darcyGradPressureFvPatchScalarField
    );
}


// ************************************************************************* //
