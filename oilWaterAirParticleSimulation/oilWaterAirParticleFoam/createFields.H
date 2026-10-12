/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/createFields.H
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

// createFields.H — oilWaterAirParticleFoam
// Three-phase immiscible mixture + particle cloud fields. Phase advection
// itself (MULES-based, see alphaEqn.H) needs no extra objects constructed
// here beyond the mixture and the flux-required flags below.

#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U", runTime.timeName(), mesh,
        IOobject::MUST_READ, IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Constructing three-phase mixture\n" << endl;

// Reads phases (air oil water) from constant/transportProperties.
// alpha1 = alpha.air, alpha2 = alpha.oil, alpha3 = alpha.water
immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

volScalarField& alpha1 = mixture.alpha1();   // air
volScalarField& alpha2 = mixture.alpha2();   // oil
volScalarField& alpha3 = mixture.alpha3();   // water

const dimensionedScalar& rho1 = mixture.rho1();   // air density
const dimensionedScalar& rho2 = mixture.rho2();   // oil density
const dimensionedScalar& rho3 = mixture.rho3();   // water density

// Dimensioned scalars for unambiguous field arithmetic.
// Use  (one_ - alpha1)  not  (1.0 - alpha1) — the latter is ambiguous
// in OpenFOAM's operator overload set and causes compile errors.
const dimensionedScalar one_ ("one",  dimless, 1.0);
const dimensionedScalar zero_("zero", dimless, 0.0);

volScalarField rho
(
    IOobject
    (
        "rho", runTime.timeName(), mesh,
        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();

surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Turbulence model — commented out entirely, not just switched to laminar.
// This solver never constructs a turbulence model object at all (contrast
// with turbulentThreePhaseFoam/interisofoamSolver, which always build one
// and just set simulationType laminar in turbulenceProperties). UEqn.H
// instead calls mixture.mu() directly for the viscous term, so there's no
// turbulence->correct() to call later in the main .C file either.
// The lines below are dead/leftover code from an earlier RAS/LES-capable
// version of this solver (not from an "isobubbleFoam" API) and are kept
// commented out only for reference if turbulence modelling is added back.
//#include "turbulentTransportModel.H"

//autoPtr<incompressible::turbulenceModel> turbulence
//(
//    incompressible::turbulenceModel::New(U, phi, mixture)
//);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// Dynamic viscosity for particle drag calculation
volScalarField muc
(
    IOobject
    (
        "muc", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    max(mixture.mu(), dimensionedScalar("mucMin", dimDynamicViscosity, 1e-5))
);

volScalarField p
(
    IOobject
    (
        "p", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label  pRefCell  = 0;
scalar pRefValue = 0.0;
setRefCell(p, p_rgh, pimple.dict(), pRefCell, pRefValue);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p", p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());
mesh.setFluxRequired(alpha2.name());
mesh.setFluxRequired(alpha3.name());

#include "createMRF.H"
#include "createFvOptions.H"

// NOTE: this solver used to construct two isoAdvection objects here
// (advector1/advector2, one per transported phase, geometric PLIC
// advection). They have been removed: alphaEqn.H now uses a MULES-based
// joint-limiter scheme for ALL THREE phases instead (see alphaEqn.H's
// header comment for why). No isoAdvection/geometricVoF objects, headers,
// or libraries are needed any more - see Make/options.
