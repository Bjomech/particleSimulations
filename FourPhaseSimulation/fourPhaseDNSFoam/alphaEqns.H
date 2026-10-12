/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/ and interFoam/interMixingFoam/ (solver family)
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

// Alpha transport with MULES subcycling (2512 API)
{
    label nAlphaCorr
    (
        mesh.solution().subDict("PIMPLE").lookupOrDefault<label>("nAlphaCorr", 1)
    );

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(mixture.cAlpha()*phic, max(phic));

    surfaceScalarField phir(phic*mixture.nHatf());

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha1
        (
            fvc::flux
            (
                phi,
                alpha1,
                "div(phi,alpha)"
            )
          + fvc::flux
            (
                phir,
                scalar(1) - alpha1,
                "div(phirb,alpha)"
            )
        );

        MULES::explicitSolve(alpha1, phi, phiAlpha1, 1, 0);

        surfaceScalarField phiAlpha2
        (
            fvc::flux
            (
                phi,
                alpha2,
                "div(phi,alpha)"
            )
          + fvc::flux
            (
               -phir,
                alpha2,
                "div(phirb,alpha)"
            )
        );

        MULES::explicitSolve(alpha2, phi, phiAlpha2, 1, 0);

        alpha3 = scalar(1) - alpha1 - alpha2;
        alpha3.max(0);
        alpha3.min(1);
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
