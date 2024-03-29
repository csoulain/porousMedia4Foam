/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::capillarityModel> Foam::capillarityModel::New
(
 const word& name,
 const dictionary& capillarityProperties,
 const volScalarField& Sb,
 const autoPtr<reducedSaturationModel> & reducedSaturationModelPtr
 )
{
  const word modelType(capillarityProperties.lookup("capillarityModel"));

  Info<< "Selecting capillarity model => " << modelType << "\n" << endl;

  dictionaryConstructorTable::iterator cstrIter =
    dictionaryConstructorTablePtr_->find(modelType);

  if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
      FatalErrorIn
	(
	 "capillarityModel::New"
	 "("
	 "const word& name,"
	 "const dictionary& capillarityProperties,"
	 "const volScalarField& Sb,"
   "const autoPtr<reducedSaturationModel> & reducedSaturationModelPtr"
	 ")"
	 )
	<< "Unknown capillarityModel type "
	<< modelType << nl << nl
	<< "Valid capillarityModels are : " << endl
	<< dictionaryConstructorTablePtr_->sortedToc()
	<< exit(FatalError);
    }

  return autoPtr<capillarityModel>
    (cstrIter()(name, capillarityProperties,Sb,reducedSaturationModelPtr));
}


// ************************************************************************* //
