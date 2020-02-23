// Solve the Momentum equation
laminarTransport.correct();

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  - fvm::laplacian(laminarTransport.nu(),U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
