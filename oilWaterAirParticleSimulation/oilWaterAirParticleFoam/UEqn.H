/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/UEqn.H
    Copyright (C) OpenFOAM Foundation and OpenCFD Ltd.; see the upstream
    file for the exact original copyright years.
    Developed and ported mainly against OpenFOAM v2512; the exact upstream
    version originally copied has not been verified.

Modifications
    Copyright (C) 2026 Bjoern Klepsch
    Last modified: 2026-09-08 (see git history for the detailed change log)
-------------------------------------------------------------------------------
License
    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See the LICENSE file in the repository root.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
\*---------------------------------------------------------------------------*/

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
//  4. Two-way particle coupling: "+ particles.SU(U)" added on the RHS.
//     KinematicCloud::SU(U) returns an implicit fvVectorMatrix sink/source
//     representing the reaction (Newton's third law) of the drag +
//     pressure-gradient + virtual-mass force the particles felt during the
//     PRECEDING particles.evolve() call (standard lagged/semi-implicit
//     coupling, exactly as in stock kinematicParcelFoam's UEqn.H). Only has
//     an effect while constant/kinematicCloudProperties has "coupled yes"
//     (with "coupled no" SU(U) returns a zero matrix, so this term is a
//     harmless no-op and the file stays correct either way).

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(mixture.mu(), U)
  - fvc::div(mixture.mu()*dev2(T(fvc::grad(U))))
 ==
    particles.SU(U)
  + fvOptions(rho, U)
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
