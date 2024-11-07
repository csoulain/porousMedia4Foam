{

    surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());
    surfaceScalarField phiPc
    (
        "phiPc",
         Mbf * fvc::interpolate(porousMedia.dpcdS(),"pc")* fvc::snGrad(Sb) * mesh.magSf()
    );


    fvScalarMatrix pEqn
    (
        fvm::laplacian(-Mf, p) + fvc::div(phiG)
        // capillary term
        + fvc::div(phiPc)
        ==
        // event source terms
        - sourceTerm
    );

    pEqn.solve();

    phiP = pEqn.flux();

    phi = phiP+phiG+phiPc;

//    phib = Fbf*phiP + (Lbf/Lf)*phiG + phiPc;
    phib = (Mbf/Mf)*phiP + (Lbf/Lf)*phiG + phiPc;
    phia = phi - phib;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Ub = fvc::reconstruct(phib);
    Ua = U-Ub;

    Ub.correctBoundaryConditions();
    Ua.correctBoundaryConditions();

}
