
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

mesh.schemes().setFluxRequired(p.name());


Info<< "Reading field sourceTerm\n" << endl;
volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sourceTerm",dimensionSet(0,0,-1,0,0,0,0),0.0)
);



/*
Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.name(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));
*/


IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

geochemicalModel geochemistry(mesh,physicalProperties);


Info << nl << "Reading/generating porosity field eps" << endl;
const volScalarField & eps = geochemistry.porosity();

Info << nl << "Reading/generating permeability field invK" << endl;
const volScalarField & invK = geochemistry.inversePermeability();


Switch activateHeleShaw
(
    physicalProperties.lookupOrDefault<Switch>("activateHeleShaw",false)
);


volScalarField HeleShawCorrection ("HeleShawCorrection", 0.0*invK);

if (activateHeleShaw)
{
    Foam::absolutePermeabilityModels::HeleShaw HeleShaw(mesh,physicalProperties);
    HeleShawCorrection = HeleShaw.inversePermeability();
    Info << "HeleShaw correction is activated" << nl << endl;
}

const volScalarField & nu = geochemistry.nu();


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
