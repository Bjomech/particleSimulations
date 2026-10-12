/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interIsoFoam/UEqn.H
    Copyright (C) OpenFOAM Foundation and OpenCFD Ltd.; see the upstream
    file for the exact original copyright years.
    Developed and ported mainly against OpenFOAM v2512; the exact upstream
    version originally copied has not been verified.

Modifications
    Copyright (C) 2026 Bjoern Klepsch
    Last modified: 2026-08-02 (see git history for the detailed change log)
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
