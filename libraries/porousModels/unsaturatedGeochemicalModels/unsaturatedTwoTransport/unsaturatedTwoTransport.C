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

#include "unsaturatedTwoTransport.H"
#include "addToRunTimeSelectionTable.H"

#include "fvMatrix.H"
//#include "fvmDdt.H"
//#include "fvmDiv.H"
//#include "fvmLaplacian.H"
#include "fvm.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace unsaturatedGeochemicalModels
    {
        defineTypeNameAndDebug(unsaturatedTwoTransport, 0);

        addToRunTimeSelectionTable
        (
            basicUnsaturatedGeochemicalModel,
            unsaturatedTwoTransport,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unsaturatedGeochemicalModels::unsaturatedTwoTransport::unsaturatedTwoTransport
(
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField &Sb,
    const incompressiblePhase &phasea,
    const incompressiblePhase &phaseb
)
:
    basicUnsaturatedGeochemicalModel(mesh, dict, Sb, phasea, phaseb),
//    unsaturatedTwoTransportDict_(dict.subDict(typeName)),
    geochemicalPropertiesDict_(dict.subDict("geochemicalProperties")),
    Sb_(Sb),
    Sa_(1.-Sb),
    massExchangeCoeff_(geochemicalPropertiesDict_.lookup("massExchangeCoefficient")),
    He_(geochemicalPropertiesDict_.lookup<scalar>("He"))
{
    Info << "initialization of the unsaturatedTwoTransport calculation ....";
    Yb_.resize(1);
    Ya_.resize(1);


    Yb_.set
    (
      0,
      new volScalarField
      (
        IOobject
        (
          "Cb",
          mesh_.time().timeName(),
          mesh_,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
        ),
        mesh_
      )
    );

    Ya_.set
    (
      0,
      new volScalarField
      (
        IOobject
        (
          "Ca",
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

void Foam::unsaturatedGeochemicalModels::unsaturatedTwoTransport::updateFluidComposition()
{

  //  Info << " Update fluid composition with unsaturatedTwoTransport" << endl;

    word divPhiYbiScheme = "div(phib,Ybi)";
    word divPhiYaiScheme = "div(phia,Yai)";



//    const volScalarField &Deff = effectiveDispersion();
    const volTensorField &DispT =  effectiveDispersionTensor();

//    Ak_ = 2*(1.-eps_)*Ak_;

    forAll(Yb_,i)
    {
      //        if(Y[i].name() != inertSpecies)
        volScalarField& Ybi = Yb_[i];
 
        volScalarField& Yai = Ya_[i];


        tmp<fvScalarMatrix> YbiEqn
        (
                  fvm::ddt(eps_*Sb_,Ybi) + fvm::div(phib_,Ybi,divPhiYbiScheme)
                - fvm::laplacian(eps_*Sb_*DispT,Ybi,"laplacian(Di,Yi)")
                ==
                  fvm::Sp(massExchangeCoeff_,Ybi)-massExchangeCoeff_*He_*Yai
        );

        YbiEqn.ref().relax();
        solve(YbiEqn);


        tmp<fvScalarMatrix> YaiEqn
        (
                  fvm::ddt(eps_*(1-Sb_),Yai) + fvm::div(phia_,Yai,divPhiYaiScheme)
                - fvm::laplacian(eps_*(1.-Sb_)*DispT,Yai,"laplacian(Di,Yi)")
                ==
                  -massExchangeCoeff_*Ybi+fvm::Sp(massExchangeCoeff_*He_,Yai)
        );

        YaiEqn.ref().relax();
        solve(YaiEqn);


    }
  //  Info<<"Ok" << endl;

/*

//volScalarField m_s ("m_s", -stoec*ae*4.0*eps*epsSolid*McaCo3*alphai*(Ceq-Ci));

volScalarField m_s ("m_s", stoec*ae*McaCo3*alphai*Ci/Ceq);

*/

}

void Foam::unsaturatedGeochemicalModels::unsaturatedTwoTransport::updateMineralDistribution()
{}
// -------------------------------------------------------------------------//

/*
Foam::volScalarField Foam::unsaturatedTwoTransport::dMl() const
{

    volScalarField dMl_(0.0*fvc::ddt(Y_[0])/this->rhol());
    forAll(Y_,s)
    {
        dMl_ = dMl_ + fvc::ddt(Y_[s])/this->rhol();
    }

    return dMl_;
}
*/


void Foam::unsaturatedGeochemicalModels::unsaturatedTwoTransport::updatedGasvdRho()
{
    dGasvdRho_ = 0.0*dGasvdRho_;
 //   forAll(mineralList_,s)
    {
        //dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./rhol_-1./rhos_[s]);
//        dMinvdRho_+= -rhos_[s]*fvc::ddt(Ys_[s])*(1./this->rho()-1./rhos_[s]);
    }
    dGasvdRho_.correctBoundaryConditions(); //necessary??
}


// ************************************************************************* //
