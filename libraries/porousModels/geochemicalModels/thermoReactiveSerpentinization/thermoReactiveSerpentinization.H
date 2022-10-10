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
    thermoReactiveSerpentinization module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef thermoReactiveSerpentinization_H
#define thermoReactiveSerpentinization_H

//#include "dictionary.H"
//#include "volFieldsFwd.H"
#include "volFields.H"
#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "basicGeochemicalModel.H"
//#include "autoPtr.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geochemicalModels
{

/*---------------------------------------------------------------------------*\
                      Class twoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermoReactiveSerpentinization : public basicGeochemicalModel
{


    //  dictionary  thermoReactiveSerpentinizationDict_;


protected:

    // Protected data

    word mineral_;

    dictionary geochemicalModelDict_;

    dictionary thermoReactiveSerpentinizationDict_;

    PtrList<dictionary> mineralSubDict_;

    PtrList<dimensionedScalar> Vm_;

    PtrList<dimensionedScalar> ki_;
    PtrList<dimensionedScalar> Acti_;
    PtrList<dimensionedScalar> aKin_;
    PtrList<dimensionedScalar> bKin_;
    PtrList<dimensionedScalar> cKin_;

    volScalarField T_;
    dimensionedScalar Cpf_;
    dimensionedScalar Cps_;
    dimensionedScalar lambda_;
    dimensionedScalar latentHeat_;

    Switch initializeTemperatureGradient_;
    Switch solveTemperature_;
    Switch solveChemistry_;

    // Private Member Functions

    void readMineralProperties();

    //- Disallow copy construct
    thermoReactiveSerpentinization(const thermoReactiveSerpentinization&);


public:

    //- Runtime type information
    TypeName("thermoReactiveSerpentinization");


    // Declare run-time constructor selection table

    // Constructors

    //- Construct from components

    thermoReactiveSerpentinization
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~thermoReactiveSerpentinization()
    {}


        // Member Functions

        void updateTemperature();

        virtual void updateFluidComposition();

        virtual void updateMineralDistribution();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace geochemicalModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
