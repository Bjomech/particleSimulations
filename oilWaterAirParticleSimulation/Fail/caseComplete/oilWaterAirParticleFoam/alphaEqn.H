/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/alphaEqn.H
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

// alphaEqn.H  (FIXED VERSION)
//
// Three-phase isoAdvector advection, called ONCE per time step before the
// PIMPLE loop.
//
// What changed and why
// ────────────────────
//  1. rhoPhi no longer double-counts the water flux.
//     isoAdvection::getRhoPhi(rhoA, rhoB) returns
//
//         (rhoA - rhoB)*dVf/deltaT + rhoB*phi
//
//     so the old
//
//         rhoPhi = advector1.getRhoPhi(rho1, rho3)
//                + advector2.getRhoPhi(rho2, rho3)
//
//     contained rho3*phi TWICE, i.e. the density flux in fvm::div(rhoPhi, U)
//     was about 2000 kg/m3 instead of 1000 kg/m3 in the whole water bulk while
//     fvm::ddt(rho, U) still used 1000. That inconsistency alone makes the
//     momentum balance (and therefore the rise velocity) wrong. Correct form:
//
//         rhoPhi = rho3*phi
//                + (rho1 - rho3)*alphaPhi1
//                + (rho2 - rho3)*alphaPhi2
//
//     which is what the expression below evaluates.
//
//  2. The hand-written "mass-conservative scale-back" of alpha1/alpha2 is gone.
//     isoAdvector already bounds and conserves alpha internally (see the
//     "isoAdvection: After conservative bounding" lines in the log). The old
//     code multiplied the WHOLE alpha field by sumBefore/sumAfter, which moves
//     volume between unrelated regions of the domain, and it clamped
//     alpha.oil to (1 - alpha.air) on every step, which deleted oil wherever
//     the bubble's interface band touched a droplet. Measured in the delivered
//     run of oilWaterAirParticle_large: oil volume dropped from 60.6 mm3 at
//     t = 0.002 s to 44.4 mm3 at t = 0.318 s, i.e. -27 %.
//
//     Replacement: clamp to [0,1], then in the (rare) cells where
//     alpha.air + alpha.oil > 1 scale BOTH by 1/(alpha.air + alpha.oil).
//     That keeps the local air/oil ratio, touches only genuinely overlapping
//     cells and leaves every other cell untouched.
//
//  3. No sub-cycling. isoAdvector is not a MULES limiter; it wants one
//     geometric advection step per deltaT with a bounded interface Courant
//     number (see maxAlphaCo in system/controlDict). Wrapping only alpha1 in a
//     subCycle also left alpha2's old-time levels inconsistent.
//     nAlphaSubCycles / nAlphaCorr / cAlpha in fvSolution are ignored here.

{
    // ── advect air ────────────────────────────────────────────────────────
    {
        volScalarField::Internal Sp1
        (
            IOobject("Sp1", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        );
        volScalarField::Internal Su1
        (
            IOobject("Su1", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        );

        advector1.advect(Sp1, Su1);
    }

    // ── advect oil ────────────────────────────────────────────────────────
    {
        volScalarField::Internal Sp2
        (
            IOobject("Sp2", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        );
        volScalarField::Internal Su2
        (
            IOobject("Su2", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        );

        advector2.advect(Sp2, Su2);
    }

    // ── consistent density flux (no duplicated rho3*phi) ──────────────────
    rhoPhi =
        advector1.getRhoPhi(rho1, rho3)
      + advector2.getRhoPhi(rho2, rho3)
      - rho3*phi;

    // ── bounding: only overlapping cells are renormalised ─────────────────
    alpha1 = max(min(alpha1, one_), zero_);
    alpha2 = max(min(alpha2, one_), zero_);

    {
        // scale = 1 wherever alpha1 + alpha2 <= 1
        const volScalarField scale(one_/max(alpha1 + alpha2, one_));
        alpha1 *= scale;
        alpha2 *= scale;
    }

    alpha3 = one_ - alpha1 - alpha2;
}

// ── conservation monitor ─────────────────────────────────────────────────
// Watch these two numbers: they must stay constant to within a fraction of a
// percent for the whole run. If the oil volume drifts down, the air/oil
// overlap treatment above is being triggered too often (mesh too coarse for
// the droplet, or deltaT too large).
{
    const scalar Vair  = gSum(alpha1.primitiveField()*mesh.V().field());
    const scalar Voil  = gSum(alpha2.primitiveField()*mesh.V().field());

    Info<< "Phase fractions:"
        << "  alpha.air = "   << alpha1.weightedAverage(mesh.V()).value()
        << " (V = " << Vair*1e9 << " mm^3)"
        << "  alpha.oil = "   << alpha2.weightedAverage(mesh.V()).value()
        << " (V = " << Voil*1e9 << " mm^3)"
        << "  min/max alpha.water = "
        << min(alpha3).value() << "/" << max(alpha3).value()
        << endl;
}
