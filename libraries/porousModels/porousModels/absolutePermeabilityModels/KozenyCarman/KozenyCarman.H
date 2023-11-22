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
    Several versions of the Kozeny-Carman equation are available in literature. We use a switch to select the suitable variant

    1. K_(t+dt) = K(t)*( (eps(t+dt)^3)/((1-eps(t+dt))^2) )*( ((1-eps(t))^2)/(eps(t)^3) )

    Ref: Implementation and evaluation of permeability-porosity and tortuosity-porosity relationships linked to mineral dissolution-precipitation. Mingliang Xie et.al. (2015) -> Eq. 2

    2. K_(t+dt) = K(0)*( ((1-eps(0))/(1-eps))^2*(eps/eps(0))^3 )

    Ref.: Benchmarking of reactive transport codes for 2D simulations with mineral dissolution-precipitation reactions and feedback on transport parameters. J. Poonoosamy et.al.(2018) -> Eq. 10
    In code: use switch 'updateFromInitialPoroPerm_'

    3. K_(t+dt) = K(0)*(eps(t)/eps(0))^3

    Ref.: Benchmarking of reactive transport codes for 2D simulations w
ith mineral dissolution-precipitation reactions and feedback on transpo
rt parameters. J. Poonoosamy et.al.(2018) -> Eq. 11
    In code: use switch 'modifiedKozenyCarman_'

SourceFiles
    KozenyCarman.C

\*---------------------------------------------------------------------------*/

#ifndef KozenyCarman_H
#define KozenyCarman_H

#include "volFields.H"
#include "speciesTable.H"

#include "absolutePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace absolutePermeabilityModels
{
/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class KozenyCarman : public absolutePermeabilityModel
{


  protected:

    dictionary KozenyCarmanDict_;

    // Protected data

    // Name of the porosity field
    word epsName_;

    dimensionedScalar K0_;

    Switch updateFromInitialValue_;
    Switch modifiedKozenyCarman_;


    volScalarField K_;
    volScalarField initK_;
    volScalarField invK_;

    surfaceScalarField Kf_;

    const volScalarField & eps_;
    volScalarField eps0_;

  // Private Member Functions

    //- Disallow copy construct
    KozenyCarman(const KozenyCarman&);


public:

  //- Runtime type information
    TypeName("KozenyCarman");


  // Declare run-time constructor selection table


  // Selectors

    //- Construct from components
    KozenyCarman
    (
      const fvMesh& mesh,
      const dictionary& dict
    );


    //- Destructor
    virtual ~KozenyCarman()
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
