// alphaEqn.H — three-phase sequential constrained isoAdvection
//
// Key fixes vs previous version:
//   1. Hard clamp replaced by mass-conservative scale-back for alpha1
//   2. rhoPhi built from both advectors (air + oil contributions)
//   3. This file must be included ONCE per timestep, BEFORE pimple.loop()
//
// Strategy:
//   STEP 1 — advect alpha.air; scale back any overshoot (conserves mass)
//   STEP 2 — pre-cap alpha.oil, advect, post-cap (oil can't enter pure-air cells)
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

        // ── STEP 2: advect alpha.oil (constrained) ────────────────────────
        alpha2 = min(max(alpha2, zero_), one_ - alpha1);
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
            alpha2 = min(max(alpha2, zero_), one_ - alpha1);
        }

        // ── STEP 3: alpha.water = remainder ──────────────────────────────
        alpha3 = max(one_ - alpha1 - alpha2, zero_);

        // ── rhoPhi: blend fluxes from both advectors ──────────────────────
        // advector1 carries the air/water interface flux (rho1 vs rho3)
        // advector2 carries the oil/water interface flux (rho2 vs rho3)
        // Summing both gives a consistent density flux for the momentum eqn.
        surfaceScalarField rhoPhiAlpha1 = advector1.getRhoPhi(rho1, rho3);
        surfaceScalarField rhoPhiAlpha2 = advector2.getRhoPhi(rho2, rho3);

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
