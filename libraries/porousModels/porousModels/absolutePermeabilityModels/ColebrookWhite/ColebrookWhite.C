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

\*---------------------------------------------------------------------------*/

#include "ColebrookWhite.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace absolutePermeabilityModels
{
    defineTypeNameAndDebug(ColebrookWhite, 0);

    addToRunTimeSelectionTable
    (
        absolutePermeabilityModel,
        ColebrookWhite,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absolutePermeabilityModels::ColebrookWhite::ColebrookWhite
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    absolutePermeabilityModel(mesh, dict),
    ColebrookWhiteDict_(dict.subDict(typeName+"Coeffs")),
    epsName_(ColebrookWhiteDict_.lookupOrDefault<word>("eps", "eps")),
    UName_(ColebrookWhiteDict_.lookupOrDefault<word>("U", "U")),
    rhoName_(ColebrookWhiteDict_.lookupOrDefault<word>("rho", "rho")),
    muName_(ColebrookWhiteDict_.lookupOrDefault<word>("mu", "mu")),
    Rw_(ColebrookWhiteDict_.lookupOrDefault
    (
        "Rw",
        dimensionedScalar("Rw",dimensionSet(0,1,0,0,0,0,0),SMALL))
    ),
    wallRoughness_(ColebrookWhiteDict_.lookupOrDefault
    (
        "wallRoughness",
        dimensionedScalar("wallRoughness",dimensionSet(0,1,0,0,0,0,0),0.0))
    ),
    K_
    (
        IOobject
        (
          "K",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        mesh_,
        Rw_*Rw_/8.,
        //calculatedFvPatchScalarField::typeName
        "zeroGradient"
    ),
    invK_
    (
        IOobject
        (
          "invK",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        1./K_,
        "zeroGradient"
    ),
    Kf_("Kf", fvc::interpolate(K_,"K")),
    U_(mesh.lookupObject<volVectorField>(UName_)),
    rho_(mesh.lookupObject<volScalarField>(rhoName_)),
    mu_(mesh.lookupObject<volScalarField>(muName_))
{

//    updatePermeability();

}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::ColebrookWhite::absolutePermeability() const
{
      return K_;
}

Foam::tmp<Foam::volScalarField>
Foam::absolutePermeabilityModels::ColebrookWhite::inversePermeability() const
{
      return invK_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::absolutePermeabilityModels::ColebrookWhite::Kf() const
{
      return Kf_;
}


void Foam::absolutePermeabilityModels::ColebrookWhite::updatePermeability()
{
    volScalarField Re ("Re", rho_*mag(U_)*2.*Rw_/mu_);

    volScalarField f ("f", 0*Re);

    forAll(Re,cellID)
    {


//        Info << "Re["<< cellID <<"] = " << Re[cellID] << nl <<endl;
        if(Re[cellID] <= 2400)
        {
            f[cellID]=64./(Re[cellID]+SMALL);
        }
        else
        {

            // Newton-Raphson algorithm to calculate the Colebrook-White
            // friction term in presence of turbulence after GSWELL

            scalar F1, F2, F3, fOld = 0;

            scalar fErr = 1;

            scalar A1  = wallRoughness_.value()/(2.*Rw_.value()*3.7);

            scalar A2 = 2.51/Re[cellID];

            scalar A3 = log(10.);

            f[cellID] = 0.02;

            while (fErr >= 0.01)
            {
                fOld = mag(f[cellID]);

                F1=1./sqrt(fOld) + 2.*log10(A1 + A2/sqrt(fOld));

                F2 = 2.*A2/(A3*(A1+A2/sqrt(fOld)));

                F3 = -0.5*pow(fOld,-1.5)*(1.+F2);

                f[cellID] = fOld - F1/F3;

                fErr = mag((f[cellID]-fOld)/f[cellID]);

        //        Info << "f[" << cellID << "]=" << f[cellID] << " fErr = " << fErr <<endl;
            }
        }
    }

    f.correctBoundaryConditions();

    invK_= f*Re/(8.*Rw_*Rw_);

    invK_.max(0.0);

    invK_.correctBoundaryConditions();

    dimensionedScalar smallInvK_("smallInvK",dimensionSet(0,-2,0,0,0),SMALL);

    K_ = 1./(invK_+smallInvK_);
    K_.correctBoundaryConditions();

    //Kf_ = fvc::interpolate(K_,"harmonic");

    Kf_ = fvc::interpolate(K_);

}
// -------------------------------------------------------------------------//
