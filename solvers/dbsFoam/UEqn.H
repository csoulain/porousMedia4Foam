    fvVectorMatrix UEqn
    (
        1./(eps+SMALL)*(fvm::ddt(U) + fvm::div(phi, U))
    //  - fvm::laplacian(turbulence->nu()/(eps+SMALL),U, "laplacian(nu,U)")
    //  + fvm::Sp(turbulence->nu()*invK,U)
    - fvm::laplacian(nu/(eps+SMALL),U, "laplacian(nu,U)")
    + fvm::Sp(nu*invK,U)
    ==
        fvModels.source(U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);
    }