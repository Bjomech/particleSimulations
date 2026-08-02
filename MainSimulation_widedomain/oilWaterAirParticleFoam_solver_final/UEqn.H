// UEqn.H — momentum equation with surface tension
MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(mixture.mu(), U)
 ==
    fvOptions(rho, U)
  + fvc::reconstruct(mixture.surfaceTensionForce() * mesh.magSf())
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p_rgh));
    fvOptions.correct(U);
}
