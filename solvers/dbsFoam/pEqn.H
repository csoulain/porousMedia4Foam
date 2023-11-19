//rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

tmp<volScalarField> rAtU(rAU);

{
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(HbyA)
//          + fvc::interpolate(rAU)*fvc::ddtCorr(rho, U, phi)
          + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAtU(), p)
          - sourceTerm
//         ==
//            betav*fvModels.source(psi, p, rho.name())
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

//#include "rhoEqn.H"
#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

//U = HbyA - (invA & (betav*fvc::grad(p)));
U = fvc::reconstruct(phi);
U.correctBoundaryConditions();
fvConstraints.constrain(U);
//K = 0.5*magSqr(U);

//if (thermo.dpdt())
//{
//    dpdt = fvc::ddt(p);
//}
