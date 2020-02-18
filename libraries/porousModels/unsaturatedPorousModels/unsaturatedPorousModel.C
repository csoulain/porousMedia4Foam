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

#include "unsaturatedPorousModel.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unsaturatedPorousModel::unsaturatedPorousModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField &Sb,
    const incompressiblePhase &phasea,
    const incompressiblePhase &phaseb
)
:
        porousModel(mesh,dict),
        mesh_(mesh),
        porousMediaDict_(dict.subDict("porousMediaProperties")),
        Sb_(Sb),
        phasea_(phasea),
        phaseb_(phaseb),
        reducedSaturationModelPtr_
        (
            reducedSaturationModel::New(porousMediaDict_,Sb)
        ),
        relativePermeabilityModelPtr_
        (
            relativePermeabilityModel::New
            (
              "krModel",
              porousMediaDict_,
              Sb,
              reducedSaturationModelPtr_
            )
        ),
        capillarityModelPtr_
        (
            capillarityModel::New
            (
              "pcModel",
              porousMediaDict_,
              Sb,
              reducedSaturationModelPtr_
            )
        ),
        Maf_
        (
          IOobject
          (
              "Faf",
              Sb_.time().timeName(),
              Sb_.db(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          Sb.mesh(),
          dimensionedScalar("Maf",dimensionSet(-1,3,1,0,0),0)
        ),
        Mbf_("Mbf",0.0*Maf_),
        Mf_("Mf",0.0*Maf_),
        Laf_
        (
          IOobject
          (
              "Faf",
              Sb_.time().timeName(),
              Sb_.db(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          Sb.mesh(),
          dimensionedScalar("Maf",dimensionSet(0,0,1,0,0),0)
        ),
        Lbf_("Lbf",0.0*Laf_),
        Lf_("Lf",0.0*Laf_)
      //  Lf_("Lf",Laf_+Lbf_)
        /*,
        Faf_
        (
          IOobject
          (
              "Faf",
              Sb_.time().timeName(),
              Sb_.db(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          Sb.mesh(),
          dimensionedScalar("Faf",dimless,0)
      ),
      Fbf_
      (
        IOobject
        (
            "Fbf",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("Fbf",dimless,0)
    )
    */
{
    this->update();
}

// -------------------------------------------------------------------------//


// ************************************************************************* //
