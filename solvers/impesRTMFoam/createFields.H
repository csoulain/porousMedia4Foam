///////////////////////////////////////////////////////////////
////////////////////// TIMESTEP CONTROL ///////////////////////
///////////////////////////////////////////////////////////////

//word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Courant");
word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
if (adjustTimeStep) adjustTimeStep=true; // to remove warnings at compilation
int CFLint = -1;
if (CFL == "Coats") CFLint = 0;
else if (CFL == "Todd") CFLint = 1;
else if (CFL == "Courant") CFLint = 2;
else
{
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "CFL condition unrecongnized : Coats, Todd and Courant available"
            << exit(FatalError);
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/*
Info << nl << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U",dimensionSet(0,1,-1,0,0),vector::zero)
);

#include "createPhi.H"
surfaceScalarField phiP = phi;
*/

//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

Info << "Reading pressure field p" << endl;
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

Info << nl << "Reading saturation field Sb" << endl;
volScalarField Sb
(
    IOobject
    (
        "Sb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",1.));

//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading transportProperties" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//- list that receives event files of event-based boundary conditions
//List<patchEventFile*> patchEventList;
//eventInfiltration::setEventFileRegistry(&patchEventList, Sb.name());

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//Info<< "Reading physicalProperties\n" << endl;
//immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

//non-wetting
autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(mesh,physicalProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

// wetting
autoPtr<incompressiblePhase> phaseb = incompressiblePhase::New(mesh,physicalProperties,"b");
volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& mub = phaseb->mu();

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

unsaturatedGeochemicalModel porousMedia(mesh,physicalProperties,Sb,phasea,phaseb);

// Porosity
Info << nl << "Reading porosity field eps" << endl;
const volScalarField & eps = porousMedia.porosity();


Info << nl << "Reading mobilities" << endl;
// absolute permeability is included in mobilities

const surfaceScalarField &Mbf = porousMedia.Mbf();
const surfaceScalarField &Lbf = porousMedia.Lbf();

const surfaceScalarField &Mf = porousMedia.Mf();
const surfaceScalarField &Lf = porousMedia.Lf();




/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
Info << nl << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Ua + Ub
);

#include "createPhi.H"
surfaceScalarField phiP = phi;

//U= Ua+Ub;
///////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS //////////////////////////
///////////////////////////////////////////////////////////////////

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
