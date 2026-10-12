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

// alphaEqn.H  (CONSERVATIVE MULES VERSION — supersedes the isoAdvector fix)
//
// Why this file was rewritten again
// ──────────────────────────────────
// The previous fix (two independent isoAdvector instances + a clamp/rescale
// "scale = 1/max(alpha1+alpha2,1)" step) bounds each phase individually but
// NOT their sum. Whenever the two independently-reconstructed geometric
// interfaces (bubble vs. droplet) happen to overlap in a cell, the clamp
// deletes volume from BOTH phases symmetrically, converting it into
// "water". This is exactly what produced the measured -8.76 % air / -7.99 %
// oil volume drift in the run that triggered this rewrite, concentrated at
// the timesteps where the rising bubble first touches an oil droplet -
// i.e. precisely the physics this thesis is about.
//
// isoAdvector is a *geometric* PLIC method; there is no widely used,
// verified "isoAdvector for N immiscible phases with joint boundedness" in
// mainline OpenFOAM (only the two-phase interIsoFoam uses it). What OpenFOAM
// itself does for N immiscible phases is an *algebraic* MULES scheme with a
// JOINT flux limiter across all phases at once. This file ports that
// scheme, verified against the real ESI/OpenCFD (openfoam.com, the fork
// v2512 belongs to) N-phase source:
//
//   Foam::multiphaseMixture::solveAlphas()
//   applications/solvers/multiphase/multiphaseInterFoam/multiphaseMixture/multiphaseMixture.C
//   https://www.openfoam.com/documentation/guides/latest/api/multiphaseMixture_8C_source.html
//
// NOTE (ESI vs Foundation API - this caused real compile errors on v2512):
// an earlier revision of this file was ported from the FOUNDATION
// (openfoam.org) reference (incompressibleMultiphaseVoF/alphaPredictor.C),
// whose MULES API differs: Foundation has a MULES::control class and a
// limitSum(alphas, alphaPhis, phi) overload; ESI v2512 has NEITHER. The ESI
// algorithm instead: (1) MULES::limit(...) with returnCorr = true, so each
// phase's flux is reduced to its limited CORRECTION; (2) a single-argument
// MULES::limitSum(alphaPhiCorrs) that jointly rescales the three correction
// fluxes; (3) each phase's bounded upwind flux is added back
// (upwind<scalar>(mesh, phi).flux(alpha)) before MULES::explicitSolve. Same
// physics/boundedness, different decomposition. Limiter tuning in ESI is
// read internally by MULES from fvSolution's "alpha.*" dict (nLimiterIter,
// default 3) - there is no control object to construct here.
//
// How it works, in one sentence: instead of advecting each phase then
// clamping the result, build all three phases' compression-corrected fluxes
// FIRST, bound each individually with MULES::limit, then call
// MULES::limitSum() to jointly rescale the three fluxes so their SUM can
// never transport more than one full cell volume in a single step - i.e.
// the three phases cannot geometrically over-claim a cell in the first
// place. This is exactly analogous to what MULES::limitSum was written for
// (see src/finiteVolume/fvMatrices/solvers/MULES/MULES.H).
//
// Trade-off you are accepting with this fix
// ──────────────────────────────────────────
// MULES with interface compression is algebraic/flux-based, not geometric.
// It is slightly more diffusive than isoAdvector's PLIC reconstruction, i.e.
// the air/oil/water interfaces will be a cell or two less sharp than
// before. This is the same trade-off stock OpenFOAM's own N-phase immiscible
// solver makes, so it is a defensible, mainstream choice for a case with
// this many simultaneous interfaces - not something invented for this
// project.
//
// What changed relative to the isoAdvector version
// ──────────────────────────────────────────────────
//  * isoAdvector (advector1/advector2, createFields.H) is no longer used at
//    all - removed from createFields.H and oilWaterAirParticleFoam.C.
//  * All three phases (air/oil/water) are now explicitly MULES-solved (the
//    old code derived alpha3 = 1 - alpha1 - alpha2 algebraically; the
//    reference N-phase algorithm solves alpha3 the same way as the other
//    two and only uses a tiny defensive renormalisation at the end).
//  * fvSchemes needs one new entry: div(phirb,alpha) Gauss interfaceCompression;
//    (added). fvSolution's "alpha.*" subdict lost its isoAdvector-only keys
//    (reconstructionScheme, isoFaceTol, ...) and kept cAlpha.
//  * No sub-cycling is implemented at all (nAlphaSubCycles in fvSolution is
//    documentation-only and not read here) - deliberately, to keep this
//    rewrite's OpenFOAM-API surface as small and close to the verified
//    reference as possible, since none of this can be compiled/tested in
//    this environment. If you need sub-cycling at a coarser mesh, say so
//    and it can be added.

{
    const dictionary& alphaControls = mesh.solverDict(alpha1.name());
    const scalar cAlpha(alphaControls.get<scalar>("cAlpha"));

    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    // Compression velocity magnitude, capped at the local face velocity -
    // identical to the reference algorithm and to stock two-phase interFoam.
    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(cAlpha*phic, max(phic));

    // Non-owning references to the three phase fields / the three phases
    // are fixed (air, oil, water); a UPtrList just lets the compression-flux
    // loop below iterate "for every OTHER phase" generically, mirroring
    // the reference's forAllIter(phases_) structure without needing the
    // phase/PtrDictionary machinery multiphaseInterFoam uses internally.
    UPtrList<volScalarField> alphaFields(3);
    alphaFields.set(0, &alpha1);
    alphaFields.set(1, &alpha2);
    alphaFields.set(2, &alpha3);

    // PtrList (owning) because MULES::limitSum needs to hold onto all three
    // correction fluxes at once - built with `new`, exactly as the
    // reference solveAlphas() does.
    PtrList<surfaceScalarField> alphaPhis(3);

    for (label i = 0; i < 3; i++)
    {
        volScalarField& alphai = alphaFields[i];

        alphaPhis.set
        (
            i,
            new surfaceScalarField
            (
                "phi" + alphai.name() + "Corr",
                fvc::flux(phi, alphai, alphaScheme)
            )
        );

        surfaceScalarField& alphaPhi = alphaPhis[i];

        // Compression term, summed pairwise over the two OTHER phases -
        // this is what actively sharpens/holds each interface and is the
        // piece that was missing entirely from a naive per-phase MULES
        // solve.
        for (label j = 0; j < 3; j++)
        {
            if (j == i) continue;

            volScalarField& alphaj = alphaFields[j];

            surfaceScalarField phir(phic*mixture.nHatf(alphai, alphaj));

            alphaPhi += fvc::flux
            (
                -fvc::flux(-phir, alphaj, alpharScheme),
                alphai,
                alpharScheme
            );
        }

        // Per-phase bound: keeps this one phase in [0,1] on its own.
        // returnCorr = true (last argument): on return alphaPhi holds only
        // the limited CORRECTION flux (full flux minus bounded upwind part),
        // which is what the ESI limitSum below expects; the upwind part is
        // added back before explicitSolve - exactly as in solveAlphas().
        MULES::limit
        (
            1.0/mesh.time().deltaTValue(),
            geometricOneField(),
            alphai,
            phi,
            alphaPhi,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField(),
            true
        );
    }

    // ── The actual fix: joint bound across ALL THREE phases at once ───────
    // Jointly rescales the three limited correction fluxes so their sum
    // cannot transport more than one cell volume of "stuff" per step, i.e.
    // it is no longer possible for air+oil to geometrically over-claim a
    // cell the way the old isoAdvector clamp discovered only AFTER the
    // fact. (ESI single-argument overload; operates on corrections, the
    // bounded upwind parts are consistent by construction.)
    MULES::limitSum(alphaPhis);

    volScalarField sumAlpha
    (
        IOobject("sumAlpha", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    for (label i = 0; i < 3; i++)
    {
        surfaceScalarField& alphaPhi = alphaPhis[i];

        // Rebuild the full phase flux: jointly-limited correction + bounded
        // upwind flux (identical to solveAlphas(); upwind.H is included by
        // oilWaterAirParticleFoam.C).
        alphaPhi += upwind<scalar>(mesh, phi).flux(alphaFields[i]);

        MULES::explicitSolve(geometricOneField(), alphaFields[i], alphaPhi);
        sumAlpha += alphaFields[i];
    }

    // Consistent density flux built directly from the three full, jointly
    // bounded phase fluxes (correction + upwind, as rebuilt above) -
    // replaces the old getRhoPhi()-based expression.
    rhoPhi = alphaPhis[0]*rho1 + alphaPhis[1]*rho2 + alphaPhis[2]*rho3;

    // Defensive renormalisation of the tiny (sub-percent, floating-point /
    // finite-MULES-tolerance level) residual left after the joint bound -
    // ported verbatim from the reference alphaSolve(), NOT the old
    // "delete-from-both-phases" clamp. With limitSum already enforcing the
    // bound above, sumAlpha should already be within ~1e-6 of 1 everywhere;
    // this only mops up what's left.
    volScalarField sumCorr(1.0 - sumAlpha);
    alpha1 += alpha1*sumCorr;
    alpha2 += alpha2*sumCorr;
    alpha3 += alpha3*sumCorr;
}

// ── conservation monitor ─────────────────────────────────────────────────
// Watch these two numbers: with the joint MULES bound above they should now
// stay flat to within a fraction of a percent for the whole run, including
// through bubble/droplet contact events.
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
