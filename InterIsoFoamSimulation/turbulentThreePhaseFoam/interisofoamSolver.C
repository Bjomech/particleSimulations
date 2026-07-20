/*---------------------------------------------------------------------------*\
  interisofoamSolver

  Same physical case as turbulentThreePhaseFoam (InterFoamSimulation): two
  phases (water/oil) + basicKinematicCloud particles, on the identical case
  geometry (oilWaterParticlesInterIsoFoam mirrors oilWaterParticlesInterFoam).
  The ONLY intended difference is the interface-capturing method:
  - turbulentThreePhaseFoam: implicit MULES + interface-compression (interFoam-style)
  - this solver: explicit geometric VOF via isoAdvector (interIsoFoam-style)
  Run side by side to compare the two VOF methods on the same setup.

  Turbulence/particle notes are identical to turbulentThreePhaseFoam — see
  that file's header (laminar in this case via turbulenceProperties; particle
  cloud active every timestep).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

#include "isoAdvection.H"           // geometric VOF interface reconstruction/advection
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"

#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

#include "basicKinematicCloud.H"   // Lagrangian particle cloud, same as turbulentThreePhaseFoam

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "VOF core of turbulentThreePhaseFoam ported to OpenFOAM v2512, "
        "using the isoAdvector geometric-VoF interface capturing scheme "
        "(interIsoFoam-based). Kinematic particle cloud included."
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

    // Reads constant/kinematicCloudProperties + kinematicCloudPositions,
    // identical files to the InterFoamSimulation case.
    const word kinematicCloudName
    (
        args.getOrDefault<word>("cloud", "kinematicCloud")
    );

    Info<< "Constructing kinematic cloud " << kinematicCloudName << endl;

    basicKinematicCloud particles
    (
        kinematicCloudName,
        rho,
        U,
        muc,
        g
    );

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
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
            }

            if (mesh.changing())
            {
                gh = (g & mesh.C()) - ghRef;
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

            // isoAdvector's own geometric advection sub-cycle, replacing the
            // MULES alpha solve used by turbulentThreePhaseFoam's pEqn-side loop.
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            rho = alpha1*rho1 + alpha2*rho2;
            rhoPhi = fvc::interpolate(rho)*phi;

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            // Laminar in this case (turbulenceProperties); see file header.
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            // Advance all Lagrangian parcels one flow timestep.
            particles.evolve();
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
