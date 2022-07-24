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





/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info<< "Reading wetting fluid velocity field Uwetting\n" << endl;
/*
volVectorField Uwetting
(
    IOobject
    (
        "Uwetting",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading wetting fluid velocity field UnonWetting\n" << endl;
volVectorField UnonWetting
(
    IOobject
    (
        "UnonWetting",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/
Info<< "Reading field U\n" << endl;
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
    dimensionedVector("U",dimensionSet(0,1,-1,0,0,0,0),vector::zero)
//    UnonWetting + Uwetting
);

#include "createPhi.H"

Info<< "Reading Solid Fraction Indicator\n" << endl;
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimensionSet(0,0,0,0,0,0,0), 0)
);

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//non-wetting
Info<< "Reading transportProperties\n" << endl;

immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

//const dimensionedScalar& rho1 = mixture.rho1();
//const dimensionedScalar& rho2 = mixture.rho2();

/*
autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(mesh,transportProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

// wetting
autoPtr<incompressiblePhase> phaseb = incompressiblePhase::New(mesh,transportProperties,"b");
volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& mub = phaseb->mu();
*/


//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

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

//dictionary wetting (transportProperties.subDict(mixture.phase1Name()));
//dictionary nonWetting (transportProperties.subDict(mixture.phase2Name()));

autoPtr<incompressiblePhase> phase1
  = incompressiblePhase::New(mesh,transportProperties,mixture.phase1Name());

Info << "TESTT after wetting 2" << nl << endl;


volVectorField& Uwetting = phase1->U();
surfaceScalarField& phiw = phase1->phi();
const dimensionedScalar& rho1 = phase1->rho();
const dimensionedScalar& mu1 = phase1->mu();
dimensionedScalar nu1("nu1",mu1/rho1);


autoPtr<incompressiblePhase> phase2
  = incompressiblePhase::New(mesh,transportProperties,mixture.phase2Name());

  Info << "TESTT after non wetting" << nl << endl;


  volVectorField& UnonWetting = phase2->U();
  surfaceScalarField& phinw = phase2->phi();
  const dimensionedScalar& rho2 = phase2->rho();
  const dimensionedScalar& mu2 = phase2->mu();
  dimensionedScalar nu2("nu2",mu2/rho2);

Info << "OK" << nl << endl;

U = UnonWetting + Uwetting;

//surfaceScalarField phiw ("phiw", linearInterpolate(Uwetting)&mesh.Sf());
//surfaceScalarField phinw ("phinw", linearInterpolate(UnonWetting)&mesh.Sf());


/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

//unsaturatedPorousModel porousMedia(mesh,transportProperties,Sb,phasea,phaseb);

unsaturatedPorousModel porousMedia(mesh,transportProperties,alpha1,phase1,phase2);


// Porosity
Info << nl << "Reading porosity field eps" << endl;
const volScalarField & eps = porousMedia.porosity();


Info << nl << "Reading mobilities" << endl;
// absolute permeability is included in mobilities

/*
const surfaceScalarField &Mbf = porousMedia.Mbf();
const surfaceScalarField &Lbf = porousMedia.Lbf();

const surfaceScalarField &Mf = porousMedia.Mf();
const surfaceScalarField &Lf = porousMedia.Lf();
*/

/*
Info<< "Reading fluid fraction field eps\n" << endl;
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/



#include "createPorousMediaFields.H"



// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

Info<< "Calculating averaged mass flux\n" << endl;
surfaceScalarField rhoPhiByEps
(
        IOobject
        (
        "rhoPhiByEps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    rhoPhi*linearInterpolate(1/(eps+SMALL))
);

// Construct incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"
