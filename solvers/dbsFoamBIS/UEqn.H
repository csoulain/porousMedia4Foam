// Solve the Momentum equation
laminarTransport.correct();

tmp<fvVectorMatrix> tUEqn
(
    1./(eps+SMALL)*(fvm::ddt(U) + fvm::div(phi, U))
  - fvm::laplacian(laminarTransport.nu()/(eps+SMALL),U)
  + fvm::Sp(laminarTransport.nu()*invK,U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
