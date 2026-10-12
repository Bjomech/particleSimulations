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

// Reading porosity properties from constant directory
IOdictionary porosityProperties
(
    IOobject
    (
        "porosityProperties",
        runTime.constant(),
        runTime,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const bool porosityEnabled
(
    porosityProperties.getOrDefault<bool>("porosityEnabled", false)
);

tmp<volScalarField> tporosity;
if (porosityEnabled)
{
    tporosity = tmp<volScalarField>::New
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
}
