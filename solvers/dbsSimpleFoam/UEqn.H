    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        1./(eps+SMALL)*fvm::div(phi, U)
      - fvm::laplacian(turbulence->nu()/(eps+SMALL),U)
      + fvm::Sp(turbulence->nu()*invK,U)
      + fvm::Sp(turbulence->nu()*HeleShawCorrection,U)
     ==
    fvModels.source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);
    }
