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

#include "simpleFirstOrderKinetic.H"
#include "addToRunTimeSelectionTable.H"

#include "fvMatrix.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvm.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace geochemicalModels
    {
        defineTypeNameAndDebug(simpleFirstOrderKinetic, 0);

        addToRunTimeSelectionTable
        (
            basicGeochemicalModel,
            simpleFirstOrderKinetic,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geochemicalModels::simpleFirstOrderKinetic::simpleFirstOrderKinetic
(
  const fvMesh& mesh,
  const dictionary& dict
)
:
    basicGeochemicalModel(mesh, dict),
    simpleFirstOrderKineticDict_(dict.subDict(typeName)),
    transportPropertiesDict_(dict),
    mineralSubDict_( mineralList_.size() ),
    Vm_( mineralList_.size() ),
    ki_( mineralList_.size() ),
    Ceq_( mineralList_.size() ),
    Di_( simpleFirstOrderKineticDict_.lookup("Di") ) //,
//    alphaL_( simpleFirstOrderKineticDict_.lookup("alphaL") )
{
    Info << "initialization of the simpleFirstOrderKinetic calculation ....";
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

void Foam::geochemicalModels::simpleFirstOrderKinetic::readMineralProperties()
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
              transportPropertiesDict_.subDict(currentMineral+"Properties")
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

    Ceq_.set
    (
        s,
        new dimensionedScalar
        (
               mineralSubDict_[s].lookup("Ceq")
        )
    );

  }

}

void Foam::geochemicalModels::simpleFirstOrderKinetic::updateFluidComposition()
{

  //  Info << " Update fluid composition with simpleFirstOrderKinetic" << endl;

    word divPhiYiScheme = "div(phi,Yi)";

//    volScalarField Deff("Deff", Foam::pow(eps_,1)*Di_*(1.+alphaL_/Di_*mag(U_)));
//    volScalarField Deff("Deff", Foam::pow(eps_,2)*Di_);
    const volTensorField &Deff = effectiveDispersionTensor();


  //  tmp<volScalarField> Aee_ (this->surfaceArea());
    volScalarField Ak_
    (
        "Ak",
        0.0*porousMedia_[0].surfaceArea()*ki_[0]/Ceq_[0]
    );

    forAll(mineralList_,s)
  	{
        Ak_ += porousMedia_[s].surfaceArea()*ki_[s]/Ceq_[s];
    }

    forAll(Y_,i)
    {
      //        if(Y[i].name() != inertSpecies)
        volScalarField& Yi = Y_[i];
      //        dimensionedScalar& Di = D[i];
      //        dimensionedScalar& Di = D[i];



        fvScalarMatrix YiEqn
        (
                  fvm::ddt(eps_,Yi) + fvm::div(phi_,Yi,divPhiYiScheme)
                - fvm::laplacian(Deff,Yi,"laplacian(Di,Yi)")
          ==
            fvm::Sp(Ak_,Yi)
        );

        YiEqn.solve();

      //        Yi.max(0.0);
      //        Yi.min(1.0);
    }
  //  Info<<"Ok" << endl;

/*

//volScalarField m_s ("m_s", -stoec*ae*4.0*eps*epsSolid*McaCo3*alphai*(Ceq-Ci));

volScalarField m_s ("m_s", stoec*ae*McaCo3*alphai*Ci/Ceq);

*/

}

void Foam::geochemicalModels::simpleFirstOrderKinetic::updateMineralDistribution()
{
    forAll(Ys_,s)
    {
        volScalarField dMs_
        (
            "dMs",
            porousMedia_[s].surfaceArea()*ki_[s]*Vm_[s]/Ceq_[s]*Y_[s]/(Ys_[s]+SMALL)
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
Foam::volScalarField Foam::simpleFirstOrderKinetic::dMl() const
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
