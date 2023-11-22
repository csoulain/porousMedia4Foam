#include "createRDeltaT.H"

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

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

pressureReference pressureReference(p, pimple.dict());

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


IOdictionary transportProperties
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

//autoPtr<geochemicalModel> geochemistry
//  = geochemicalModel::New(mesh,transportProperties);

geochemicalModel geochemistry(mesh,transportProperties);


Info << nl << "Reading/generating porosity field eps" << endl;
const volScalarField & eps = geochemistry.porosity();

Info << nl << "Reading/generating permeability field invK" << endl;
const volScalarField & invK = geochemistry.inversePermeability();


Switch activateHeleShaw
(
    transportProperties.lookupOrDefault<Switch>("activateHeleShaw",false)
);


volScalarField HeleShawCorrection ("HeleShawCorrection", 0.0*invK);

if (activateHeleShaw)
{
    Foam::absolutePermeabilityModels::HeleShaw HeleShaw(mesh,transportProperties);
    HeleShawCorrection = HeleShaw.inversePermeability();
    Info << "HeleShaw correction is activated" << nl << endl;
}

const volScalarField & nu = geochemistry.nu();


/*
Info << nl << "Reading/generating momentumTransportModel";
//autoPtr<viscosityModel> viscosity(viscosityModel::New(mesh));


autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, geochemistry.viscosity())
);
Info << " ... OK" << endl;
*/


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
