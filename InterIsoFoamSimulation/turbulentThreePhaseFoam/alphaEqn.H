/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interIsoFoam/alphaEqn.H
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

// alphaEqn.H — one geometric-VOF advection substep via isoAdvector.
// Unlike interFoam's alphaEqn.H, there is no matrix solve here: advector.advect()
// reconstructs the interface geometrically per cell and fluxes alpha1 by exact
// volume intersection. Sp/Su below are unused sink/source placeholders kept
// only for interface compatibility with fvOptions.
{
    volScalarField::Internal Sp
    (
        IOobject
        (
            "Sp",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    volScalarField::Internal Su
    (
        IOobject
        (
            "Su",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero)
    );

    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // Core isoAdvector step: reconstructs the water/oil interface as a plane
    // in each interfacial cell and computes the exact geometric flux of alpha1.
    advector.advect(Sp, Su);

    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    // rhoPhi built directly from the geometric alpha flux (not from fvc::interpolate(rho)*phi).
    rhoPhi = advector.getRhoPhi(rho1, rho2);

    alpha2 = 1.0 - alpha1;   // two phases only — oil fraction is just the complement
    mixture.correct();
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
