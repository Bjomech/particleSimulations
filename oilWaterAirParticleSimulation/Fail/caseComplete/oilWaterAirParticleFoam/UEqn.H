// UEqn.H  (FIXED VERSION)
//
// Momentum equation in stock interFoam/interIsoFoam form.
//
// What changed and why
// ────────────────────
//  1. Surface tension, buoyancy and the pressure gradient are now ALL applied
//     as face-based (snGrad) terms inside one fvc::reconstruct(...) call.
//     The old version mixed a reconstructed capillary force with a cell-based
//     fvc::grad(p_rgh) and left buoyancy out of the predictor completely.
//     Mixing cell-gradient and face-gradient discretisations of terms that
//     must cancel each other (Laplace pressure jump vs pressure gradient) is
//     exactly what produces spurious interface currents strong enough to stop
//     a small bubble from rising.
//  2. Added the missing deviatoric viscous term. With mu ratios of ~50 (water
//     vs oil) and ~55 (water vs air) the div(mu*dev2(T(grad(U)))) term is not
//     negligible; stock interFoam always includes it via divDevRhoReff.
//     NOTE: add this to system/fvSchemes divSchemes:
//         div((mu*dev2(T(grad(U))))) Gauss linear;
//  3. fvOptions.constrain(UEqn) added (was missing).

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(mixture.mu(), U)
  - fvc::div(mixture.mu()*dev2(T(fvc::grad(U))))
 ==
    fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
