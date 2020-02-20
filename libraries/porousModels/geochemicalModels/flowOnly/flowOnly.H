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
    flowOnly module

SourceFiles
    phreeqc.C

\*---------------------------------------------------------------------------*/

#ifndef flowOnly_H
#define flowOnly_H

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

class flowOnly : public basicGeochemicalModel
{


  protected:

    // Protected data



  // Private Member Functions

    //- Disallow copy construct
    flowOnly(const flowOnly&);


  public:

    //- Runtime type information
      TypeName("flowOnly");


    // Declare run-time constructor selection table

  // Constructors

    //- Construct from components

    flowOnly
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~flowOnly()
    {}


    // Member Functions

    virtual void updateFluidComposition()
    {}

    virtual void updateMineralDistribution()
    {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace geochemicalModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
