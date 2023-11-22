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

#include "simpleFirstOrderKineticMole.H"
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
        defineTypeNameAndDebug(simpleFirstOrderKineticMole, 0);

        addToRunTimeSelectionTable
        (
            basicGeochemicalModel,
            simpleFirstOrderKineticMole,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geochemicalModels::simpleFirstOrderKineticMole::simpleFirstOrderKineticMole
(
  const fvMesh& mesh,
  const dictionary& dict
)
:
    basicGeochemicalModel(mesh, dict),
    geochemicalModelDict_(dict.subDict("geochemicalProperties")),
    mineralSubDict_( mineralList_.size() ),
    Vm_( mineralList_.size() ),
    ki_( mineralList_.size() ),
    Acti_( mineralList_.size() )
{
    Info << "initialization of the simpleFirstOrderKineticMole calculation ....";
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

    readMineralProperties();
    Info<< "OK" << nl << endl;
}


// -------------------------------------------------------------------------//

void Foam::geochemicalModels::simpleFirstOrderKineticMole::readMineralProperties()
{

  forAll(mineralList_,s)
	{
		word currentMineral = mineralList_[s];
		Info << " Doing stuff for mineral: " << currentMineral << endl;

    mineralSubDict_.set
    (
        s,
        new dictionary
        (
              geochemicalModelDict_.subDict(currentMineral+"Properties")
        )
    );

		Vm_.set
		(
  			s,
  			new dimensionedScalar
  			(
  					   mineralSubDict_[s].lookup("Vm")
  			)
		);

    ki_.set
    (
        s,
        new dimensionedScalar
        (
               mineralSubDict_[s].lookup("ki")
        )
    );

    Acti_.set
    (
        s,
        new dimensionedScalar
        (
               mineralSubDict_[s].lookup("Acti")
        )
    );

  }

}

void Foam::geochemicalModels::simpleFirstOrderKineticMole::updateFluidComposition()
{

  //  Info << " Update fluid composition with simpleFirstOrderKineticMole" << endl;

    word divPhiYiScheme = "div(phi,Yi)";

    const volTensorField &Deff = effectiveDispersionTensor();


  //  tmp<volScalarField> Aee_ (this->surfaceArea());
    volScalarField Ak_
    (
        "Ak",
        0.0*porousMedia_[0].surfaceArea()*ki_[0]*Acti_[0]
    );

    forAll(mineralList_,s)
  	{
        Ak_ += porousMedia_[s].surfaceArea()*ki_[s]*Acti_[s];
    }

//    Ak_ = 2*(1.-eps_)*Ak_;

    forAll(Y_,i)
    {
      //        if(Y[i].name() != inertSpecies)
        volScalarField& Yi = Y_[i];
      //        dimensionedScalar& Di = D[i];

        tmp<fvScalarMatrix> YiEqn
        (
                  fvm::ddt(eps_,Yi) + fvm::div(phi_,Yi,divPhiYiScheme)
                - fvm::laplacian(Deff,Yi,"laplacian(Di,Yi)")
          ==
            fvm::Sp(Ak_,Yi)
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

void Foam::geochemicalModels::simpleFirstOrderKineticMole::updateMineralDistribution()
{
    forAll(Ys_,s)
    {
        volScalarField dMs_
        (
            "dMs",
            porousMedia_[s].surfaceArea()*ki_[s]*Vm_[s]*Acti_[s]*Y_[s]/(Ys_[s]+SMALL)
        );

        solve
        (
          fvm::ddt(Ys_[s]) == fvm::Sp(dMs_,Ys_[s])
        );

      //    Ys_[s].max(0.0);
      //    Ys_[s].min(0.999);
    }
}
// -------------------------------------------------------------------------//

/*
Foam::volScalarField Foam::simpleFirstOrderKineticMole::dMl() const
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
