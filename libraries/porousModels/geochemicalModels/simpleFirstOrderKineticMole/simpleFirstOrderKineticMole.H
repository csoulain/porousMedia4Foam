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
    simpleFirstOrderKineticMole module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef simpleFirstOrderKineticMole_H
#define simpleFirstOrderKineticMole_H

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

class simpleFirstOrderKineticMole : public basicGeochemicalModel
{


//  dictionary  simpleFirstOrderKineticMoleDict_;


  protected:

    // Protected data

    word mineral_;

    dictionary geochemicalModelDict_;

    PtrList<dictionary> mineralSubDict_;

    PtrList<dimensionedScalar> Vm_;

    PtrList<dimensionedScalar> ki_;
    PtrList<dimensionedScalar> Acti_;



  // Private Member Functions

    void readMineralProperties();

    //- Disallow copy construct
    simpleFirstOrderKineticMole(const simpleFirstOrderKineticMole&);


  public:

    //- Runtime type information
      TypeName("simpleFirstOrderKineticMole");


    // Declare run-time constructor selection table

  // Constructors

    //- Construct from components

    simpleFirstOrderKineticMole
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~simpleFirstOrderKineticMole()
    {}


    // Member Functions

    virtual void updateFluidComposition();

    virtual void updateMineralDistribution();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace geochemicalModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
