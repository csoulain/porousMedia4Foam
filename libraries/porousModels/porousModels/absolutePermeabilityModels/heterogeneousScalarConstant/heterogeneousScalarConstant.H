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
    Foam::twoPhaseMixture

Description
    PhreeqcRM module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousScalarConstant_H
#define heterogeneousScalarConstant_H

#include "volFields.H"
#include "speciesTable.H"
//#include "basicMultiComponentMixture.H"
//#include "runTimeSelectionTables.H"

#include "absolutePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace absolutePermeabilityModels
{
/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousScalarConstant : public absolutePermeabilityModel
{

  protected:

    // Protected data

    dictionary heterogeneousScalarConstantDict_;

    dimensionedScalar K0_;

    volScalarField K_;

    volScalarField invK_;

    surfaceScalarField Kf_;


  // Private Member Functions

    //- Disallow copy construct
    heterogeneousScalarConstant(const heterogeneousScalarConstant&);


public:

  //- Runtime type information
    TypeName("heterogeneousScalarConstant");


  // Declare run-time constructor selection table


  // Selectors

    //- Construct from components
    heterogeneousScalarConstant
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~heterogeneousScalarConstant()
    {}


    // Member Functions
    virtual tmp<volScalarField> absolutePermeability() const;

    virtual tmp<volScalarField> inversePermeability() const;

    virtual tmp<surfaceScalarField> Kf() const;

    virtual void updatePermeability();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace absolutePermeabilityModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
