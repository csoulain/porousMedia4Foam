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

#include "taylorArisTensor.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionTensorModels
{
    defineTypeNameAndDebug(taylorArisTensor, 0);

    addToRunTimeSelectionTable
    (
        dispersionTensorModel,
        taylorArisTensor,
        dictionary
    );
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dispersionTensorModels::taylorArisTensor::taylorArisTensor
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dispersionTensorModel(mesh, dict),
    taylorArisTensorDict_(dict.subDict(typeName+"Coeffs")),
    UName_(taylorArisTensorDict_.lookupOrDefault<word>("U", "U")),
    Di_( taylorArisTensorDict_.lookup("Di") ),
    Deff_
     (
      IOobject
      (
       "Deff",
       mesh_.time().timeName(),
       mesh_,
       IOobject::READ_IF_PRESENT,
       IOobject::NO_WRITE
      ),
      mesh_,
      dimensionedTensor("Deff",dimensionSet(0, 2, -1, 0, 0),tensor::zero),
      "zeroGradient"
    ),
    U_(mesh.lookupObject<volVectorField>(UName_)),
    d_( taylorArisTensorDict_.lookup("channelDiameter") )

{
}

// * * * * * * * * * * * * * * member functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField>
Foam::dispersionTensorModels::taylorArisTensor::effectiveDispersionTensor() const
{
      return Deff_;
}

void Foam::dispersionTensorModels::taylorArisTensor::updateDispersionTensor()
{
//      Deff_= Di_*(1 + Foam::pow(d_*0.5,2)*Foam::pow(mag(U_),2)/(48*Foam::pow(Di_,2)));


//      Info << "max U_ = " << max(U_) << nl << endl;


//      Info << "Mean U_ = " << fvc::domainIntegrate(U_)/sum(mesh_.V()) << nl << endl;

//      Info << "channelDiameter = " << d_.value() << nl << endl;

//      Info << "Di = " <<  Di_.value() << nl << endl;


      volScalarField Pe  ("Pe", 0.5*d_*mag(U_)/Di_);


      dimensionedScalar time = U_.mesh().time();

      volScalarField dispCorr ("dispCorr", Pe*0.0);
      scalar pi = 3.141592653589793;


      for(int i=1 ; i <= 3 ; i++)
      {
          dispCorr += 18./pow(i*pi,6)*Foam::exp(-Foam::pow(i*pi/(0.5*d_),2)*Di_*time);
      }

  //    Deff_=
  //      tensor(1,0,0,0,0,0,0,0,0)*Di_*(1. + Foam::pow(Pe,2)*mag(2./105-dispCorr))
  //    + tensor(0,0,0,0,1,0,0,0,0)*Di_*(1. + Foam::pow(Pe,2)*mag(2./105-dispCorr))

       Deff_= tensor(1,0,0,0,1,0,0,0,0)*Di_ + Foam::pow(0.5*d_,2)/Di_*U_*U_*mag(2./105-0.0*dispCorr)

      ;
}
// -------------------------------------------------------------------------//
