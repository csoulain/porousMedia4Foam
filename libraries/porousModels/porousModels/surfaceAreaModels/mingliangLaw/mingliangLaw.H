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

Class
    Foam::surfaceAreaModels::mingliangLaw

Description
    power law function to compute surface area

    Ae = A0*(Ys/Y0)^n 

SourceFiles
    mingliangLaw.C

Authors
    05/03/2020 : S. Pavuluri
\*---------------------------------------------------------------------------*/

#ifndef mingliangLaw_H
#define mingliangLaw_H

#include "volFields.H"
#include "surfaceAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceAreaModels
{
/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class mingliangLaw : public surfaceAreaModel
{


  protected:

    // Protected data
    dictionary mingliangLawDict_;

    dimensionedScalar A0_;

    scalar n_;
    scalar initVolFrac_;
    volScalarField Ae_;

    const volScalarField & Ys_;



  // Private Member Functions

    //- Disallow copy construct
    mingliangLaw(const mingliangLaw&);


public:

  //- Runtime type information
    TypeName("mingliangLaw");


  // Declare run-time constructor selection table


  // Selectors

    //- Construct from components
    mingliangLaw
    (
      const fvMesh& mesh,
      const volScalarField& Ys,
      const dictionary& dict
    );


    //- Destructor
    virtual ~mingliangLaw()
    {}


    // Member Functions
    virtual tmp<volScalarField> surfaceArea() const;

    virtual void updateSurfaceArea();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace absolutePermeabilityModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
