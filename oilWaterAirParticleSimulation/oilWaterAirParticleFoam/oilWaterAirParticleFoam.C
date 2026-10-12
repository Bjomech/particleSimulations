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
  Three-phase MULES/interface-compression VOF + Lagrangian kinematic
  particle cloud. (Earlier revision used two independent isoAdvector
  instances; replaced because it advected air and oil geometrically but
  independently, so overlapping reconstructed interfaces at bubble/droplet
  contact caused an -8.76 %/-7.99 % volume drift. See alphaEqn.H's header
  comment for the full rationale and the real OpenFOAM source it is ported
  from.)

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
     the interface-advection stability limit was never enforced.
  2. muc (the viscosity the cloud sees) is now simply the mixture viscosity.
     The old blend
         muc = alpha1*1.81e-5 + alpha2*50e-3 + alpha3*mixture.mu()
     was wrong twice: it hard-coded viscosities that disagree with
     constant/transportProperties (mu2 = 0.0425 there, 0.05 here), and the
     third term is quadratic in alpha3 because mixture.mu() already contains
     alpha1*mu1 + alpha2*mu2 + alpha3*mu3. mixture.mu() is exactly the
     intended three-phase blend.
  3. Two-way coupling: ENABLED. UEqn.H now has "+ particles.SU(U)" in the
     momentum equation, and constant/kinematicCloudProperties has
     "coupled yes". SU(U) is KinematicCloud's own implicit source term - the
     Newton's-third-law reaction of the drag/pressure-gradient/virtual-mass
     force the fluid applied to the parcels during the PRECEDING
     particles.evolve() call (standard lagged coupling, as in stock
     kinematicParcelFoam). See kinematicCloudProperties' header comment for
     why this has negligible visible effect at this particle count/size (St
     ~ 5e-3, three 14 um parcels) despite being numerically active and
     physically correct - it is orthogonal to, and does not replace, a
     dedicated interface-attachment ParticleForce if you want visible
     trapping.
  4. Interface-attachment force: IMPLEMENTED. A new custom ParticleForce,
     InterfaceAttachmentForce (InterfaceAttachmentForce.H/.C, registered via
     makeOilWaterAirParticleForces.C, selected as "interfaceAttachment" in
     kinematicCloudProperties), gives particles the sub-grid capillary
     pinning behaviour that plain drag/gravity/pressureGradient/virtualMass
     cannot provide at St ~ 5e-3. This supersedes the note in (3)/README
     that such a force "needs to be added" - it now exists. See
     InterfaceAttachmentForce.H's header comment for the model, its units,
     and citations.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "MULES.H"
#include "upwind.H"       // upwind<scalar>(mesh, phi).flux(alpha) in alphaEqn.H
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
        "oilWaterAirParticleFoam: three-phase MULES VOF solver "
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
        // MULES-based, jointly bounded across all three phases - reads
        // cAlpha (and optional MULES::control tuning keys) directly from
        // fvSolution's "alpha.*" subdict inside alphaEqn.H itself.
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
