//word phaseName("");
//if (args.optionFound("phase")) phaseName = args.option("phase");

Info << nl << "Reading field p\n" << endl;
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


Info << nl << "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U",dimensionSet(0,2,-1,0,0,0,0),vector::zero)
);


Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


dimensionedScalar rho (transportProperties.lookup("rho"));

dimensionedScalar mu (transportProperties.lookup("mu"));


//autoPtr<incompressiblePhase> phase = incompressiblePhase::New(mesh,transportProperties,phaseName);
//volVectorField& UphaseName = phase->U();
//const dimensionedScalar& rho = phase->rho();
//const dimensionedScalar& mu = phase->mu();
//volVectorField U("U",UphaseName);

#include "createPhi.H"


geochemicalModel porousMedia(mesh,transportProperties);


Info << "Reading field K\n" << endl;
const surfaceScalarField& Kf = porousMedia.Kf();

surfaceScalarField Mf ("Mf",Kf/mu);
surfaceScalarField phiG("phiG",(rho * Mf * g) & mesh.Sf());
surfaceScalarField phiPc("phiPc", 0 * phiG); //- required for porousBoundaryConditions


//- For constant or event source injection/extraction
volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimless/dimTime,0)
);
