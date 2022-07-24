volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

/*
surfaceScalarField phiDragU
(
    fvc::interpolate(laminarTransport.nu()*invK,"invK")*phi
);
*/

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

//phiHbyA += phiDragU;

if (p.needReference())
{
    adjustPhi(phiHbyA, U, p);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU());

// Source term in the continuity equation due to phase change
tmp<volScalarField> source(geochemistry.dMinvdRho());


// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) - sourceTerm*0.001//- source
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = fvc::reconstruct(phi);
//U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
