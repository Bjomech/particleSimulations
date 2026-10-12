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

// threePhaseInterfaceProperties
//
// Description
//     Computes per-phase interface curvature (K1, K2, K3) and the resulting
//     surface tension force for a three-phase (air/oil/water) VOF mixture,
//     using the three pairwise surface tension coefficients sigma12, sigma13,
//     sigma23 read from incompressibleThreePhaseMixture (constant/transportProperties).
//     Modelled after interfaceProperties (stock two-phase interFoam) extended
//     to three phases: each phase gets its own curvature field, and
//     surfaceTensionForce() combines them with the sigma coefficients — see
//     the comment in threePhaseInterfaceProperties.C for how each sigma is
//     paired with which curvature term.

namespace Foam
{

class threePhaseInterfaceProperties
{
private:
    const incompressibleThreePhaseMixture& mixture_;

    dimensionedScalar sigma12_;   // air-oil surface tension [N/m]
    dimensionedScalar sigma13_;   // air-water surface tension [N/m]
    dimensionedScalar sigma23_;   // oil-water surface tension [N/m]

    volScalarField K1_;   // curvature of the air phase interface (-div(nHat1))
    volScalarField K2_;   // curvature of the oil phase interface
    volScalarField K3_;   // curvature of the water phase interface

public:
    threePhaseInterfaceProperties(const incompressibleThreePhaseMixture& mixture);

    // Combined surface tension force (face-interpolated), added to pEqn.H's
    // pressure equation as a source term.
    tmp<surfaceScalarField> surfaceTensionForce() const;

    // Recomputes K1_, K2_, K3_ from the current alpha fields; call once per
    // outer PIMPLE iteration after alpha is updated.
    void correct();
};

} // End namespace Foam

#endif
