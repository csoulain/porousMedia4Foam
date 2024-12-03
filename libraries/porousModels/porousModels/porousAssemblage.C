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

#include "porousAssemblage.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousAssemblage::porousAssemblage
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
        mesh_(mesh),
        mineralList_(dict.lookup("mineral")),
        Ys_(mineralList_.size() ),
        inertMineral_
        (
            IOobject
            (
                "inertMineral",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("inertMineral",dimless,0),
            "zeroGradient"
        ),
        eps_
        (
            IOobject
            (
                "eps",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("eps",dimless,1.0),
            "zeroGradient"
        ),
        activateUpdatePorosity_(dict.lookupOrDefault("activateUpdatePorosity",true)),
        rhos_(mineralList_.size() ),
        mineralSurfaceArea_(mineralList_.size() ),
        mineral_(mineralList_.size()),
        absolutePermeabilityModelPtr_
        (
            absolutePermeabilityModel::New(mesh, dict)
        ),
        dispersionTensorModelPtr_
        (
            dispersionTensorModel::New(mesh, dict)
        )
{
    forAll(mineralList_,s)
    {
      word currentMineral = mineralList_[s];
      Info << " Doing stuff for mineral: " << currentMineral << endl;

      Ys_.set
      (
        s,
        new volScalarField
        (
          IOobject
          (
            "Ys."+mineralList_[s],
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
          ),
          mesh_ 
        )
      );
      //Ys_[s].write();

      rhos_.set // a bouger dans mineralModel
      (
          s,
          new dimensionedScalar
          (
              dict.subDict(currentMineral+"Properties").lookup("rhos")
          )
      );

      mineral_.set
      (
        s,
        new mineralModel
        (
          mesh,
          mineralList_[s],
          Ys_[s],
          dict
        )
      );

    }
    updatePorosity();
}


// -------------------------------------------------------------------------//

void Foam::porousAssemblage::updatePorosity()
{
    Info << "update dans porous assemblage" << endl;
    eps_ = 0.0*eps_;
    forAll(mineralList_,s)
    {
        eps_+=Ys_[s];
    }
    eps_ = 1.-eps_-inertMineral_;
    eps_.correctBoundaryConditions(); //necessary??
}



// ************************************************************************* //
