/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/interFoam.C
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

/*---------------------------------------------------------------------------*\
  oilWaterAirParticleFoam  (FIXED VERSION)
  Three-phase isoAdvector VOF + Lagrangian kinematic particle cloud.

  Physics goal
  ────────────
  Air bubble rises through water, entrains solid particles and carries them
  upward through a field of oil droplets.

  Structure (unchanged)
  ─────────────────────
  alphaEqn.H  : once per time step, BEFORE pimple.loop()
  particles.evolve() : once per time step, AFTER pimple.loop()

  What changed in this file
  ────────────────────────
  1. Interface (alpha) Courant number is now limited via the stock
     alphaCourantNo.H, so deltaT respects maxAlphaCo from system/controlDict.
     Previously alphaCoNum was hard-wired to 0 and maxAlphaCo to GREAT, i.e.
     the isoAdvector stability limit was never enforced.
  2. muc (the viscosity the cloud sees) is now simply the mixture viscosity.
     The old blend
         muc = alpha1*1.81e-5 + alpha2*50e-3 + alpha3*mixture.mu()
     was wrong twice: it hard-coded viscosities that disagree with
     constant/transportProperties (mu2 = 0.0425 there, 0.05 here), and the
     third term is quadratic in alpha3 because mixture.mu() already contains
     alpha1*mu1 + alpha2*mu2 + alpha3*mu3. mixture.mu() is exactly the
     intended three-phase blend.
  3. Optional two-way coupling note: with "coupled no" in
     constant/kinematicCloudProperties the cloud does not act back on the
     fluid - and UEqn.H has no particles.SU(U) term either, so momentum
     coupling would be silently ignored even with "coupled yes". If you want
     two-way coupling, set "coupled yes" and add "+ particles.SU(U)" to the
     UEqn matrix.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleThreePhaseMixture.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "basicKinematicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "oilWaterAirParticleFoam: three-phase isoAdvector VOF solver "
        "(air bubble rising through oil droplets in water) "
        "with a Lagrangian kinematic particle cloud."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    const word kinematicCloudName
    (
        args.getOrDefault<word>("cloud", "kinematicCloud")
    );

    Info<< "Constructing kinematic cloud: " << kinematicCloudName << nl;
    basicKinematicCloud particles(kinematicCloudName, rho, U, muc, g);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"   // needs maxAlphaCo in controlDict
            #include "setDeltaT.H"
        }

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ── Alpha advection: once per time step, before the PIMPLE loop ────
        // (no alphaControls.H: nAlphaCorr/nAlphaSubCycles/cAlpha are MULES
        //  parameters and are not used by the isoAdvector path)
        #include "alphaEqn.H"

        // Update density - do NOT overwrite rhoPhi (set in alphaEqn.H)
        rho = alpha1*rho1 + alpha2*rho2 + alpha3*rho3;

        // ── PIMPLE loop ───────────────────────────────────────────────────
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
            }

            if (mesh.changing())
            {
                gh  = (g & mesh.C())  - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;
                MRF.update();

                if (correctPhi)
                {
                    phi = mesh.Sf() & Uf();
                    #include "correctPhi.H"
                    fvc::makeRelative(phi, U);
                }
            }

            // Recompute the interface curvatures K.air/K.oil/K.water
            mixture.correct();

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }

        // ── Particle cloud: once per time step, after the flow solve ───────
        // muc = three-phase mixture viscosity (air inside the bubble, oil in
        // the droplets, water elsewhere), floored to avoid a zero denominator
        // in the drag correlation.
        muc = max
        (
            mixture.mu(),
            dimensionedScalar("mucMin", dimDynamicViscosity, 1e-6)
        );

        particles.evolve();

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
