/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interIsoFoam/pEqn.H
    Copyright (C) OpenFOAM Foundation and OpenCFD Ltd.; see the upstream
    file for the exact original copyright years.
    Developed and ported mainly against OpenFOAM v2512; the exact upstream
    version originally copied has not been verified.

Modifications
    Copyright (C) 2026 Bjoern Klepsch
    Last modified: 2026-07-20 (see git history for the detailed change log)
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

// pEqn.H — pressure correction with buoyancy + surface tension (same as turbulentThreePhaseFoam's pEqn.H)
tmp<volScalarField> trAU(1.0/UEqn.A());
const volScalarField& rAU = trAU();

surfaceScalarField rAUf
(
    "rAUf",
    fvc::interpolate(rAU)
);

volVectorField HbyA
(
    constrainHbyA
    (
        rAU*UEqn.H(),
        U,
        p_rgh
    )
);

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi))
);

MRF.makeRelative(phiHbyA);

if (p_rgh.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p_rgh);
    fvc::makeAbsolute(phiHbyA, U);
}

surfaceScalarField phig
(
    (
        mixture.surfaceTensionForce()
      + ghf*fvc::snGrad(rho)
    )*rAUf*mesh.magSf()
);

phiHbyA += phig;

constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference
    (
        pRefCell,
        getRefCellValue(p_rgh, pRefCell)
    );

    p_rghEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - p_rghEqn.flux();

        p_rgh.relax();

        U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
        U.correctBoundaryConditions();

        fvOptions.correct(U);
    }
}

#include "continuityErrs.H"

fvc::makeRelative(phi, U);

p == p_rgh + rho*gh;

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

p_rgh == p - rho*gh;
