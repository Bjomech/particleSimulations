/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/interMixingFoam/threePhaseInterfaceProperties/threePhaseInterfaceProperties.H
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

#ifndef threePhaseInterfaceProperties_H
#define threePhaseInterfaceProperties_H

#include "incompressibleThreePhaseMixture.H"

// threePhaseInterfaceProperties  (FIXED VERSION)
//
// Description
//     Per-phase interface curvature (K.air, K.oil, K.water) and the combined
//     CSF surface-tension force for a three-phase (air/oil/water) VOF mixture.
//
//     Curvature is computed exactly the way stock OpenFOAM's two-phase
//     interfaceProperties::calculateK() does it, i.e. from the *face*
//     interpolated gradient:
//
//         gradAlphaf = interpolate(grad(alpha))
//         nHatfv     = gradAlphaf/(mag(gradAlphaf) + deltaN)
//         nHatf      = nHatfv & Sf
//         K          = -div(nHatf)
//
//     This is what makes the CSF force "well balanced" against the
//     face-based pressure gradient in pEqn.H. Computing K from a smoothed
//     *cell* gradient and then masking it (previous version) breaks that
//     balance and produces a large spurious net force on a closed interface.
//
//     Pairwise sigma handling: with three immiscible phases the CSF force is
//     written as sum_i sigma_i * K_i * grad(alpha_i), where the per-phase
//     coefficients follow from requiring sigma_i + sigma_j = sigma_ij:
//
//         sigma1 = (sigma12 + sigma13 - sigma23)/2      (air)
//         sigma2 = (sigma12 + sigma23 - sigma13)/2      (oil)
//         sigma3 = (sigma13 + sigma23 - sigma12)/2      (water)
//
//     At a pure air/water interface alpha.oil = 0, grad(alpha.air) =
//     -grad(alpha.water) and K.air = -K.water, so the sum collapses to
//     sigma13*K.air*grad(alpha.air) - exactly the two-phase result. The same
//     holds for the oil/water and air/oil interfaces, so sigma12 (air-oil) is
//     now actually used instead of being ignored.

namespace Foam
{

class threePhaseInterfaceProperties
{
private:

    const incompressibleThreePhaseMixture& mixture_;

    dimensionedScalar sigma12_;   // air-oil   surface tension [N/m]
    dimensionedScalar sigma13_;   // air-water surface tension [N/m]
    dimensionedScalar sigma23_;   // oil-water surface tension [N/m]

    // Per-phase CSF coefficients derived from the three pairwise values
    dimensionedScalar sigma1_;
    dimensionedScalar sigma2_;
    dimensionedScalar sigma3_;

    // Stabilisation of the unit normal, scaled with mesh size exactly as in
    // stock interfaceProperties: 1e-8/pow(average(mesh.V()), 1/3)
    dimensionedScalar deltaN_;

    volScalarField K1_;   // curvature of the air   interface   ("K.air")
    volScalarField K2_;   // curvature of the oil   interface   ("K.oil")
    volScalarField K3_;   // curvature of the water interface   ("K.water")

    //- Curvature of one phase, stock-style (face-based)
    void calculateK(const volScalarField& alpha, volScalarField& K) const;

public:

    threePhaseInterfaceProperties
    (
        const incompressibleThreePhaseMixture& mixture
    );

    //- Combined surface-tension force as a face field, used BOTH in the
    //  momentum predictor (UEqn.H) and in phiHbyA (pEqn.H)
    tmp<surfaceScalarField> surfaceTensionForce() const;

    //- Recompute K.air/K.oil/K.water from the current alpha fields
    void correct();

    // ── Pairwise interface-normal (added for the MULES-based alphaEqn.H) ──
    //
    // Face unit interface-normal VECTOR between two phase fractions alpha1
    // and alpha2, used as the compression direction in the conservative
    // MULES alpha advection. Identical formula to stock OpenFOAM's
    // multiphaseMixture::nHatfv(alpha1, alpha2) (verified against
    // OpenFOAM-10 source, applications/solvers/multiphase/multiphaseInterFoam/
    // multiphaseMixture/multiphaseMixture.C):
    //
    //     gradAlphaf = interpolate(alpha2)*interpolate(grad(alpha1))
    //                - interpolate(alpha1)*interpolate(grad(alpha2))
    //     nHatfv     = gradAlphaf/(mag(gradAlphaf) + deltaN)
    //
    // Note this is deliberately NOT the same as the single-phase nHatfv used
    // internally by calculateK() above (which is built from grad(alpha) of
    // ONE phase only, for the CSF curvature). This pairwise version gives
    // the normal of the INTERFACE BETWEEN two specific phases, which is what
    // the compression flux in alphaEqn.H needs.
    tmp<surfaceVectorField> nHatfv
    (
        const volScalarField& alpha1,
        const volScalarField& alpha2
    ) const;

    //- Face-normal flux of nHatfv (nHatfv & Sf). Used directly as the
    //  compression-flux direction in MULES::limit's phir term.
    tmp<surfaceScalarField> nHatf
    (
        const volScalarField& alpha1,
        const volScalarField& alpha2
    ) const;

    // Access
    const volScalarField& K1() const { return K1_; }
    const volScalarField& K2() const { return K2_; }
    const volScalarField& K3() const { return K3_; }

    // Pairwise surface tensions (added for InterfaceAttachmentForce).
    // NOTE: InterfaceAttachmentForce (a Lagrangian ParticleForce, see
    // InterfaceAttachmentForce.H) cannot actually call these: it is built
    // by OpenFOAM's dictionary-based run-time selection, and this
    // threePhaseInterfaceProperties object is a plain local C++ object in
    // oilWaterAirParticleFoam.C, not a registered IOobject in the mesh
    // database - a ParticleForce can only look up *registered* fields by
    // name (as it does for alpha.air/alpha.oil themselves). These getters
    // are kept public for consistency/debugging and possible future
    // refactoring (e.g. registering this object), but the force currently
    // reads its own sigma13/sigma23 dict entries, defaulted to exactly the
    // same physical values (see constant/kinematicCloudProperties).
    const dimensionedScalar& sigma12() const { return sigma12_; }
    const dimensionedScalar& sigma13() const { return sigma13_; }
    const dimensionedScalar& sigma23() const { return sigma23_; }
};

} // End namespace Foam

#endif
