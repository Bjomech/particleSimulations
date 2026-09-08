/*---------------------------------------------------------------------------*\
  oilWaterAirParticleFoam
  Three-phase isoAdvector VOF + Lagrangian kinematic particle cloud.

  Physics goal
  ────────────
  Air bubble rises through water carrying a particle upward through oil
  droplets.  We want to see:
    • Sharp air|water interface (bubble shape preserved)
    • Sharp oil|water interfaces (droplets don't dissolve)
    • Particle captured inside bubble, rising with it
    • Bubble-oil and oil-particle surface interactions

  Key structural decisions
  ────────────────────────
  alphaEqn.H is included ONCE per timestep, BEFORE the pimple.loop().
  isoAdvector runs exactly once per dt — not once per outer corrector.

  rhoPhi is built inside alphaEqn.H and must NOT be overwritten after
  the include — it discards the geometrically-correct isoAdvector flux.

  particles.evolve() is called ONCE per timestep, AFTER pimple.loop().

  muc uses a full three-phase viscosity blend:
    alpha1 (air)   muAir  = 1.81e-5 Pa.s  low drag  → particle entrained
    alpha2 (oil)   muOil  = 50e-3   Pa.s  high drag → particle captured
    alpha3 (water) mixture.mu()            reference phase
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentTransportModel.H"
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
        "with Lagrangian kinematic particle cloud."
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

    // alpha3 (water) is declared in createFields.H as:
    //   volScalarField& alpha3 = mixture.alpha3();
    // No redeclaration needed here.

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
            scalar alphaCoNum = 0.0;
            scalar maxAlphaCo = GREAT;
            #include "setDeltaT.H"
        }

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ── Alpha advection ── ONCE per timestep, before PIMPLE loop ─────
        #include "alphaControls.H"
        #include "alphaEqn.H"

        // Update density — do NOT overwrite rhoPhi (set in alphaEqn.H).
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

            mixture.correct();

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }

            if (pimple.frozenFlow()) continue;

            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            // turbulence->correct() omitted: there is no turbulence model
            // object to correct in this solver at all (see createFields.H
            // — the turbulence model construction is commented out; UEqn.H
            // uses mixture.mu() directly, i.e. always laminar). This is not
            // caused by basicKinematicCloud, which has no bearing on
            // turbulence modelling. Add a turbulence model object back in
            // createFields.H (and call ->correct() here) if RAS/LES is needed.
        }

        // ── Particle cloud ── ONCE per timestep, after flow solve ─────────
        // Three-phase viscosity blend — alpha3 is mixture.alpha3() (water).
        muc = alpha1 * dimensionedScalar("muAir",  dimDynamicViscosity, 1.81e-5)
            + alpha2 * dimensionedScalar("muOil",  dimDynamicViscosity, 50e-3)
            + alpha3 * mixture.mu();
        muc = max
        (
            muc,
            dimensionedScalar("mucMin", dimDynamicViscosity, 1e-5)
        );

        particles.evolve();

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
