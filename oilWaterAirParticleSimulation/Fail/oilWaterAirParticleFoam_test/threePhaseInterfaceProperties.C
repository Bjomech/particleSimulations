/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/interMixingFoam/threePhaseInterfaceProperties/threePhaseInterfaceProperties.C
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

#include "threePhaseInterfaceProperties.H"
#include "fvc.H"
#include "fvCFD.H"

namespace Foam
{

threePhaseInterfaceProperties::threePhaseInterfaceProperties
(
    const incompressibleThreePhaseMixture& mixture
)
:
    mixture_(mixture),
    sigma12_
    (
        "sigma12",
        dimensionSet(1,0,-2,0,0),
        mixture_.lookup("sigma12")
    ),
    sigma13_
    (
        "sigma13",
        dimensionSet(1,0,-2,0,0),
        mixture_.lookup("sigma13")
    ),
    sigma23_
    (
        "sigma23",
        dimensionSet(1,0,-2,0,0),
        mixture_.lookup("sigma23")
    ),
    K1_
    (
        IOobject
        (
            "K1",
            mixture_.alpha1().time().timeName(),
            mixture_.alpha1().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mixture_.alpha1().mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    K2_(K1_),
    K3_(K1_)
{}

tmp<surfaceScalarField> threePhaseInterfaceProperties::surfaceTensionForce() const
{
    // Each term is (interpolated sigma*curvature) * snGrad(alpha) for that
    // phase, summed to give the total surface-tension pressure-equation source.
    // Coefficient pairing (matches which two phases border each interface):
    //   phase 1 (air)  term uses sigma13_ (air-water) with K1_ (air curvature)
    //   phase 2 (oil)  term uses sigma23_ (oil-water)  with K2_ (oil curvature)
    //   phase 3 (water) term uses sigma23_ (oil-water) with K3_ (water curvature)
    // sigma12_ (air-oil) is not used directly here; its effect enters through
    // how alpha2 is derived (alpha2 = 1 - alpha1 - alpha3) and how the phases
    // are ordered in incompressibleThreePhaseMixture. This is the same
    // sigma/curvature pairing scheme as OpenFOAM's stock threePhaseInterfaceProperties.
    return
        fvc::interpolate(sigma13_*K1_)*fvc::snGrad(mixture_.alpha1())
      + fvc::interpolate(sigma23_*K2_)*fvc::snGrad(mixture_.alpha2())
      + fvc::interpolate(sigma23_*K3_)*fvc::snGrad(mixture_.alpha3());
}

void threePhaseInterfaceProperties::correct()
{
    // Compute interface curvature K = -div(nHat) for each phase.
    // nHat = grad(alpha) / (|grad(alpha)| + delta)  — stabilised unit normal.
    // The small delta prevents division by zero in single-phase cells where
    // grad(alpha) = 0.

    // FIXED: deltaN was hardcoded to a bare 1e-8 (units 1/m), completely
    // unscaled by mesh resolution. Stock OpenFOAM (interfaceProperties.C)
    // instead scales it as 1e-8/(mean cell size), i.e.
    // 1e-8/pow(average(mesh.V()), 1/3). For this mesh's ~14um cells that
    // gives deltaN ~= 7e-4 (1/m) -- about 70,000x LARGER than the bare 1e-8
    // that was here. With deltaN far too small, magGradAlpha barely above
    // machine noise near the tail of the (sharp, isoAdvector) interface was
    // still enough to dominate the denominator, so near-noise gradients got
    // normalised into full-length, effectively random unit normals -- these
    // feed straight into curvature (K = -div(nHat)) and then into
    // surfaceTensionForce(), injecting large, chaotic, sign-flipping forces
    // into the momentum equation ('spurious currents', a well-known VOF/CSF
    // artifact, worsened here by the small bubble radius + realistic
    // air-water sigma13=0.07 N/m -> small Eo -> surface tension already
    // dominates buoyancy, so this noise was easily strong enough to mask or
    // reverse the real (weak) buoyant rise).
    const scalar meshSize
    (
        Foam::pow(Foam::average(mixture_.alpha1().mesh().V()).value(), 1.0/3.0)
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        dimless/dimLength,
        1e-8/meshSize
    );

    auto computeK = [&](const volScalarField& alpha) -> volScalarField
    {
        // FIXED (v3): smooth alpha before differentiating it for curvature.
        //
        // v1/v2 masked the SYMPTOM (spurious curvature) after computing it
        // from the raw, isoAdvector-sharp alpha field. That reduced but did
        // not fix the underlying noise: interface-band curvature was still
        // both wrong-sign in a large fraction of cells (~40-60% split, for
        // what should be a single-sign curvature on a roughly circular
        // bubble) and several times too large in magnitude. A controlled
        // test with sigma set to 0 (no surface tension at all) showed the
        // bubble rises fast and strongly -- reaching ~94% of the domain
        // height in 0.18s -- while the same case with real sigma and v2
        // masking stayed completely frozen for the full 0.6s run. That
        // rules out "small real signal buried in noise": the buoyant force
        // is large, and the raw curvature-driven surface tension force was
        // effectively cancelling it outright, most likely via spatially
        // incoherent (sign-flipping cell-to-cell) forces rather than a
        // simple wrong magnitude.
        //
        // This applies nSmoothIter passes of simple topological averaging
        // (each cell's value becomes the face-interpolated average of
        // itself and its neighbours) to a LOCAL COPY of alpha, used only
        // for this curvature calculation. This is the same idea used by
        // "smoothed VOF"/CSF implementations to tame curvature noise from
        // sharp interface-capturing schemes: differentiating a smoother
        // field gives a smoother, more consistently-signed gradient/
        // curvature. It does NOT touch the real alpha field used for
        // advection, rho, or mu anywhere else in the solver -- alpha
        // itself (the function parameter) is never reassigned, only this
        // separate alphaSmoothed copy is.
        const label nSmoothIter = 2;

        volScalarField alphaSmoothed
        (
            IOobject
            (
                "alphaSmoothed",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha
        );

        for (label iter = 0; iter < nSmoothIter; ++iter)
        {
            alphaSmoothed = fvc::average(fvc::interpolate(alphaSmoothed));
            alphaSmoothed.correctBoundaryConditions();
        }

        volVectorField gradAlpha(fvc::grad(alphaSmoothed));

        volScalarField magGradAlpha
        (
            IOobject
            (
                "magGradAlpha",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mag(gradAlpha)
        );

        // Stabilised unit normal pointing out of the phase
        volVectorField nHat
        (
            IOobject
            (
                "nHat",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gradAlpha / (magGradAlpha + deltaN)
        );

        // Curvature = -divergence of unit normal
        volScalarField K
        (
            IOobject
            (
                "K_tmp",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -fvc::div(nHat)
        );

        // FIXED (v2, retained in v3): mask curvature to a narrow band around
        // the real interface. Bulk cells (alpha ~ 0 or alpha ~ 1, no real
        // interface present) should have exactly zero curvature -- any
        // nonzero K there is numerical noise. v2 found that a mask based on
        // alpha in [1e-6, 1-1e-6] was too permissive (still let through
        // near-pure cells carrying the largest spikes); [0.01, 0.99] fixed
        // that. In v3 the mask is evaluated on alphaSmoothed (not the raw
        // alpha) since curvature is now computed from the smoothed field --
        // smoothing spreads the transition over a few more cells than the
        // raw isoAdvector front, so the mask needs to track the same field
        // that actually has a nonzero gradient, or it would clip the smoothed
        // curvature right back down to a sharp, noisy-again band.
        const dimensionedScalar alphaMin("alphaMin", dimless, 0.01);
        const dimensionedScalar alphaMax("alphaMax", dimless, 0.99);

        volScalarField interfaceMask
        (
            IOobject
            (
                "interfaceMask",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pos0(alphaSmoothed - alphaMin) * pos0(alphaMax - alphaSmoothed)
        );

        K *= interfaceMask;

        return K;
    };

    K1_ = computeK(mixture_.alpha1());
    K2_ = computeK(mixture_.alpha2());
    K3_ = computeK(mixture_.alpha3());
}

} // End namespace Foam
