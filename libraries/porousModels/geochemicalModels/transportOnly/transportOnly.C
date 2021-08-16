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

#include "transportOnly.H"
#include "addToRunTimeSelectionTable.H"

#include "fvMatrix.H"
//#include "fvmDdt.H"
//#include "fvmDiv.H"
//#include "fvmLaplacian.H"
#include "fvm.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace geochemicalModels
    {
        defineTypeNameAndDebug(transportOnly, 0);

        addToRunTimeSelectionTable
        (
            basicGeochemicalModel,
            transportOnly,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geochemicalModels::transportOnly::transportOnly
(
  const fvMesh& mesh,
  const dictionary& dict
)
:
    basicGeochemicalModel(mesh, dict),
//    transportOnlyDict_(dict.subDict(typeName)),
    transportPropertiesDict_(dict)
{
    Info << "initialization of the transportOnly calculation ....";
    Y_.resize(1);

    Y_.set
    (
      0,
      new volScalarField
      (
        IOobject
        (
          "Ci",
          mesh_.time().timeName(),
          mesh_,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
        ),
        mesh_
      )
    );

    Info<< "OK" << nl << endl;
}


// -------------------------------------------------------------------------//

void Foam::geochemicalModels::transportOnly::updateFluidComposition()
{

  //  Info << " Update fluid composition with transportOnly" << endl;

    word divPhiYiScheme = "div(phi,Yi)";

//    const volScalarField &Deff = effectiveDispersion();
    const volTensorField &DispT =  effectiveDispersionTensor();

//    Ak_ = 2*(1.-eps_)*Ak_;

    forAll(Y_,i)
    {
      //        if(Y[i].name() != inertSpecies)
        volScalarField& Yi = Y_[i];
      //        dimensionedScalar& Di = D[i];


        tmp<fvScalarMatrix> YiEqn
        (
                  fvm::ddt(eps_,Yi) + fvm::div(phi_,Yi,divPhiYiScheme)
                - fvm::laplacian(DispT,Yi,"laplacian(Di,Yi)")
//                - fvm::laplacian(Deff,Yi,"laplacian(Di,Yi)")
        );

        YiEqn.ref().relax();
        solve(YiEqn);

      //        Yi.max(0.0);
      //        Yi.min(1.0);
    }
  //  Info<<"Ok" << endl;

/*

//volScalarField m_s ("m_s", -stoec*ae*4.0*eps*epsSolid*McaCo3*alphai*(Ceq-Ci));

volScalarField m_s ("m_s", stoec*ae*McaCo3*alphai*Ci/Ceq);

*/

}

void Foam::geochemicalModels::transportOnly::updateMineralDistribution()
{}
// -------------------------------------------------------------------------//

/*
Foam::volScalarField Foam::transportOnly::dMl() const
{

    volScalarField dMl_(0.0*fvc::ddt(Y_[0])/this->rhol());
    forAll(Y_,s)
    {
        dMl_ = dMl_ + fvc::ddt(Y_[s])/this->rhol();
    }

    return dMl_;
}
*/

// ************************************************************************* //
