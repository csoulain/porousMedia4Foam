    // Momentum predictor

    tmp<fvVectorMatrix> tUEqn
    (
        1./(eps+SMALL)*fvm::div(phi, U)
      - fvm::laplacian(laminarTransport.nu()/(eps+SMALL),U)
      + fvm::Sp(laminarTransport.nu()*invK,U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();


    UEqn.relax();

    fvOptions.constrain(UEqn);


    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
