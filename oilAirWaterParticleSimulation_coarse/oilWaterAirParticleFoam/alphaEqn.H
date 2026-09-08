// alphaEqn.H — three-phase sequential constrained isoAdvection
//
// Key fixes vs previous version:
//   1. Hard clamp replaced by mass-conservative scale-back for alpha1
//      AND alpha2 (oil gets the same treatment as air — see STEP 2).
//      Previously alpha2 was only hard-clamped to [0, 1-alpha1] with no
//      rescale, which silently deleted oil volume every time the air
//      bubble's footprint overlapped a cell that held oil. Over a full
//      run this looked like oil droplets numerically "vanishing" instead
//      of physically interacting with the bubble. Fixed by mirroring
//      alpha1's conservative rescale, but bounded by the oil's actual
//      available room in a cell (1 - alpha1) instead of 1.
//   2. rhoPhi built from both advectors (air + oil contributions)
//   3. This file must be included ONCE per timestep, BEFORE pimple.loop()
//
// Strategy:
//   STEP 1 — advect alpha.air; scale back any overshoot (conserves mass)
//   STEP 2 — pre-cap alpha.oil, advect, post-cap; rescale back after each
//            cap so oil volume is conserved (oil still can't enter
//            cells already occupied by air — that constraint is physical)
//   STEP 3 — alpha.water = 1 - alpha.air - alpha.oil  (algebraic remainder)

{
    dimensionedScalar totalDeltaT = runTime.deltaT();

    surfaceScalarField rhoPhiSum
    (
        IOobject("rhoPhiSum", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        // ── STEP 1: advect alpha.air ──────────────────────────────────────
        {
            volScalarField::Internal Sp1
            (
                IOobject("Sp1", runTime.timeName(), mesh),
                mesh, dimensionedScalar(dimless/dimTime, Zero)
            );
            volScalarField::Internal Su1
            (
                IOobject("Su1", runTime.timeName(), mesh),
                mesh, dimensionedScalar(dimless/dimTime, Zero)
            );
            advector1.advect(Sp1, Su1);

            // Mass-conservative scale-back:
            // isoAdvector may produce alpha1 slightly > 1 at the interface.
            // Hard clipping destroys that mass; instead we measure the total
            // before and after clamping and rescale to return the lost mass
            // proportionally to all cells that still hold air.
            scalar sumBefore =
                gSum(alpha1.primitiveField() * mesh.V().field());

            alpha1 = max(min(alpha1, one_), zero_);

            scalar sumAfter =
                gSum(alpha1.primitiveField() * mesh.V().field());

            if (sumAfter > SMALL && mag(sumBefore - sumAfter) > SMALL)
            {
                alpha1 *= sumBefore / sumAfter;
                alpha1 = min(alpha1, one_);   // safety — stays bounded
            }
        }

        // ── STEP 2: advect alpha.oil (constrained, mass-conservative) ─────
        // Oil's available room in a cell is (1 - alpha1), not 1, since air
        // may already occupy part of it. Same rescale-back logic as
        // alpha1's STEP 1: measure total oil volume before/after each clamp
        // and rescale to recover mass a hard clamp would otherwise delete
        // (e.g. when the rising air bubble sweeps into a cell holding oil).
        // Oil still physically cannot enter a cell already fully occupied
        // by air — only the *loss of total oil volume* is being fixed here,
        // not the exclusion rule itself.
        {
            scalar sumBefore2 =
                gSum(alpha2.primitiveField() * mesh.V().field());

            alpha2 = min(max(alpha2, zero_), one_ - alpha1);

            scalar sumAfter2 =
                gSum(alpha2.primitiveField() * mesh.V().field());

            if (sumAfter2 > SMALL && mag(sumBefore2 - sumAfter2) > SMALL)
            {
                alpha2 *= sumBefore2 / sumAfter2;
                alpha2 = min(alpha2, one_ - alpha1);   // safety — stays bounded
            }
        }
        {
            volScalarField::Internal Sp2
            (
                IOobject("Sp2", runTime.timeName(), mesh),
                mesh, dimensionedScalar(dimless/dimTime, Zero)
            );
            volScalarField::Internal Su2
            (
                IOobject("Su2", runTime.timeName(), mesh),
                mesh, dimensionedScalar(dimless/dimTime, Zero)
            );
            advector2.advect(Sp2, Su2);

            scalar sumBefore2 =
                gSum(alpha2.primitiveField() * mesh.V().field());

            alpha2 = min(max(alpha2, zero_), one_ - alpha1);

            scalar sumAfter2 =
                gSum(alpha2.primitiveField() * mesh.V().field());

            if (sumAfter2 > SMALL && mag(sumBefore2 - sumAfter2) > SMALL)
            {
                alpha2 *= sumBefore2 / sumAfter2;
                alpha2 = min(alpha2, one_ - alpha1);   // safety — stays bounded
            }
        }

        // ── STEP 3: alpha.water = remainder ──────────────────────────────
        alpha3 = max(one_ - alpha1 - alpha2, zero_);

        // ── rhoPhi: blend fluxes from both advectors ──────────────────────
        // advector1 carries the air/water interface flux (rho1 vs rho3)
        // advector2 carries the oil/water interface flux (rho2 vs rho3)
        // Summing both gives a consistent density flux for the momentum eqn.
        // Direct initialization (not copy-initialization with '=') avoids an
        // ambiguous tmp<surfaceScalarField> -> surfaceScalarField conversion
        // that newer Clang (macOS) rejects but GCC (Linux) silently accepts.
        // Semantics are identical on both compilers/platforms.
        surfaceScalarField rhoPhiAlpha1(advector1.getRhoPhi(rho1, rho3));
        surfaceScalarField rhoPhiAlpha2(advector2.getRhoPhi(rho2, rho3));

        rhoPhiSum += (runTime.deltaT()/totalDeltaT)
                   * (rhoPhiAlpha1 + rhoPhiAlpha2);
    }

    rhoPhi = rhoPhiSum;
}

Info<< "Phase fractions:"
    << "  alpha.air = "
    << alpha1.weightedAverage(mesh.V()).value()
    << "  Min=" << min(alpha1).value()
    << "  Max=" << max(alpha1).value()
    << "  alpha.oil = "
    << alpha2.weightedAverage(mesh.V()).value()
    << "  Max=" << max(alpha2).value()
    << endl;
