// --- Coats (2003) IMPES STABILITY: CFL LIMIT
{
    dimensionedScalar Mmu(mub/mua);//viscosity ratio

    const volScalarField &K = porousMedia.absolutePermeability();
    // permeability interpolation
    const surfaceScalarField & Kf = porousMedia.Kf();

    const volScalarField& kra = porousMedia.kra();
    const volScalarField& krb = porousMedia.krb();

    const volScalarField& dkradS = porousMedia.dkradS();
    const volScalarField& dkrbdS = porousMedia.dkrbdS();


    porousMedia.update();

    // - inertial part
    volScalarField dfw
    (
        "dfw",
        (dkrbdS*kra - dkradS*krb)
        /(Mmu*Foam::pow(kra,2)+2*kra*krb+1./Mmu*Foam::pow(krb,2))
    );
    dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);

    // - gravity part
    dfw -= K*(rhoa-rhob)*fvc::surfaceSum(mag(mesh.Sf() & g))
        /fvc::surfaceSum(mag(phi)+smallRate)
        * (Foam::pow(kra,2)*dkrbdS/mua
        + Foam::pow(krb,2)*dkradS/mub)/(Mmu*Foam::pow(kra,2)
        +2*kra*krb+1/Mmu*Foam::pow(krb,2));


    scalarField CFLCoats((runTime.deltaT()/eps)*dfw*fvc::surfaceSum(mag(phi)));

    // - capillarity part
    if(porousMedia.activateCapillarity())
    {
        CFLCoats +=
        (runTime.deltaT()/eps)*2*mag(porousMedia.dpcdS())
        *fvc::surfaceSum(Kf*mesh.magSf()/mag(mesh.delta()))
        *(kra*krb/(mub*kra+mua*krb));
    }

    CFLCoats /= mesh.V();
    CFLUse = gMax(CFLCoats);
    maxDeltaTFact = maxCo/(CFLUse + SMALL);

    Info<< "Coats CFL Number mean: " << gAverage(CFLCoats) << " max: " << CFLUse << endl;
}
