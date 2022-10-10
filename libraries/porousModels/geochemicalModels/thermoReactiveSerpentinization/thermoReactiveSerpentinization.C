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

#include "thermoReactiveSerpentinization.H"
#include "addToRunTimeSelectionTable.H"

#include "fvMatrix.H"
//#include "fvmDdt.H"
//#include "fvmDiv.H"
//#include "fvmLaplacian.H"
#include "fvm.H"
#include "fvcDdt.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace geochemicalModels
    {
        defineTypeNameAndDebug(thermoReactiveSerpentinization, 0);

        addToRunTimeSelectionTable
        (
            basicGeochemicalModel,
            thermoReactiveSerpentinization,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geochemicalModels::thermoReactiveSerpentinization::thermoReactiveSerpentinization
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
basicGeochemicalModel(mesh, dict),
geochemicalModelDict_(dict.subDict("geochemicalProperties")),
thermoReactiveSerpentinizationDict_(geochemicalModelDict_.subDict(typeName)),
mineralSubDict_( mineralList_.size() ),
Vm_( mineralList_.size() ),
ki_( mineralList_.size() ),
Acti_( mineralList_.size() ),
aKin_( mineralList_.size() ),
bKin_( mineralList_.size() ),
cKin_( mineralList_.size() ),
T_
(
    IOobject
    (
      "T",
      mesh_.time().timeName(),
      mesh_,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
    ),
    mesh
),
Cpf_(thermoReactiveSerpentinizationDict_.lookup("Cpf")),
Cps_(thermoReactiveSerpentinizationDict_.lookup("Cps")),
lambda_(thermoReactiveSerpentinizationDict_.lookup("lambda")),
latentHeat_(thermoReactiveSerpentinizationDict_.lookup("latentHeat")),
initializeTemperatureGradient_
(
    thermoReactiveSerpentinizationDict_.lookupOrDefault("initializeTemperatureGradient",true)
),
solveTemperature_
(
    thermoReactiveSerpentinizationDict_.lookupOrDefault("solveTemperature",true)
),
solveChemistry_
(
    thermoReactiveSerpentinizationDict_.lookupOrDefault("solveChemistry",true)
)
{
    Info << "initialization of the thermoReactiveSerpentinization calculation ....";
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


    if (initializeTemperatureGradient_)
    {
        volScalarField Tini
        (
            IOobject
            (
              "T0",
              mesh_.time().timeName(),
              mesh_,
              IOobject::MUST_READ,
              IOobject::NO_WRITE
            ),
            mesh
        );

        solve(fvm::laplacian(Tini),"Tini");
        T_ = Tini;
        T_.correctBoundaryConditions();
    }
}


// -------------------------------------------------------------------------//

void Foam::geochemicalModels::thermoReactiveSerpentinization::readMineralProperties()
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

        aKin_.set
        (
            s,
            new dimensionedScalar
            (
                mineralSubDict_[s].lookup("a")
            )
        );
        bKin_.set
        (
            s,
            new dimensionedScalar
            (
                mineralSubDict_[s].lookup("b")
            )
        );
        cKin_.set
        (
            s,
            new dimensionedScalar
            (
                mineralSubDict_[s].lookup("c")
            )
        );
    }

}

void Foam::geochemicalModels::thermoReactiveSerpentinization::updateFluidComposition()
{
    /*
    if(solveChemistry_)
    {
        //  Info << " Update fluid composition with thermoReactiveSerpentinization" << endl;

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
    }
    */

    if(solveTemperature_)
    {
        updateTemperature();
    }
}


void Foam::geochemicalModels::thermoReactiveSerpentinization::updateTemperature()
{

    surfaceScalarField phiRhoCp ("phiRhoCp", phi_*fvc::interpolate(this->rho())*Cpf_);

    tmp<volScalarField> rhoCp = (1.-eps_)*rhos_[0]*Cps_+eps_*this->rho()*Cpf_;

    tmp<fvScalarMatrix> TEqn
    (
        rhoCp*fvm::ddt(T_) + fvm::div(phiRhoCp,T_) //- T*fvc::div(phiRhoCp)
        - fvm::laplacian(lambda_,T_) +rhos_[0]*latentHeat_*fvc::ddt(eps_)
//        ==
//        fvm::Sp(Ak_,Yi)
    );

/*
    surfaceScalarField phiCp ("phiCp", phi_*Cpf_);

    tmp<volScalarField> Cp = (1.-eps_)*Cps_+eps_*Cpf_;

    tmp<fvScalarMatrix> TEqn
    (
        fvm::ddt(Cp,T_) + fvm::div(phiCp,T_)
        - fvm::laplacian(lambda_,T_)+latentHeat_*fvc::ddt(eps_)
//        ==
//        fvm::Sp(Ak_,Yi)
    );
*/

    solve(TEqn);

}


void Foam::geochemicalModels::thermoReactiveSerpentinization::updateMineralDistribution()
{
    if(solveChemistry_)
    {
        forAll(Ys_,s)
        {

            tmp<volScalarField> kk =
                aKin_[s]*Foam::exp(-bKin_[s]*Foam::pow((T_-cKin_[s]),2));
//                eps_*aKin_[s];

            volScalarField dMs_
            (
                "dMs",
//                porousMedia_[s].surfaceArea()*ki_[s]*Vm_[s]*Acti_[s]*Y_[s]/(Ys_[s]+SMALL)
//                -porousMedia_[s].surfaceArea()*kk*Vm_[s]*Y_[s]/(Ys_[s]+SMALL)
                (eps_-0.1)*kk // /(Ys_[s]+SMALL)
            );

            solve
            (
                fvm::ddt(Ys_[s]) == fvm::Sp(dMs_,Ys_[s]), "Ys"
            );

/*
            tmp<fvScalarMatrix> YsEqn
            (
                fvm::ddt(Ys_[s]) == fvm::Sp(dMs_,Ys_[s]), "Ys"
            );
            solve(YsEqn);
            */
            //    Ys_[s].max(0.0);
            //    Ys_[s].min(0.999);
        }
    }
}
// -------------------------------------------------------------------------//

/*
Foam::volScalarField Foam::thermoReactiveSerpentinization::dMl() const
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
