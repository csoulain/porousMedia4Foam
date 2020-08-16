//word phaseName("");
//if (args.optionFound("phase")) phaseName = args.option("phase");

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

#include "createPhi.H"


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




geochemicalModel porousMedia(mesh,transportProperties);
