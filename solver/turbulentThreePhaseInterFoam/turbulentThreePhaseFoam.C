/*---------------------------------------------------------------------------*\
  turbulentThreePhaseFoam

  NOTE ON THE NAME: despite the name, this solver is currently TWO-PHASE,
  not three. It is stock OpenFOAM interFoam (MULES-based VOF, standard
  immiscibleIncompressibleTwoPhaseMixture) with a Lagrangian
  basicKinematicCloud bolted on for particle tracking. constant/
  transportProperties defines only "phases (water oil)" — there is no
  air/third phase in this case. The name was inherited from an earlier,
  more ambitious version of the case and never updated.

  Role in this thesis: this is the baseline oil-water VOF + particle setup,
  used to compare the standard MULES/interface-compression alpha advection
  (this solver) against the geometric isoAdvector method used by
  interisofoamSolver on the identical case geometry (oilWaterParticlesInterFoam
  vs oilWaterParticlesInterIsoFoam).

  Turbulence: an incompressibleInterPhaseTransportModel IS constructed
  (see createFields.H) and turbulence->correct() IS called every outer
  corrector below — but constant/turbulenceProperties sets
  "simulationType laminar;", so in this case it runs laminar. The model
  object supports switching to RAS/LES later just by changing that entry.

  Particle cloud: basicKinematicCloud is constructed once before the time
  loop and evolved once per PIMPLE iteration via particles.evolve() — it IS
  active in this file (see below), not disabled.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"

#include "immiscibleIncompressibleTwoPhaseMixture.H"   // two-phase (water/oil) VOF mixture
#include "incompressibleInterPhaseTransportModel.H"     // laminar/turbulent switch, see turbulenceProperties
#include "turbulentTransportModel.H"

#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"

#include "basicKinematicCloud.H"   // Lagrangian particle cloud, added on top of stock interFoam

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "VOF core of turbulentThreePhaseFoam ported to OpenFOAM v2512 "
        "(actually two-phase water/oil — see file header). "
        "Kinematic particle cloud is active and evolved every timestep."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    // Cloud reads its config from constant/kinematicCloudProperties and its
    // injection points from constant/kinematicCloudPositions.
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
        muc,   // carrier dynamic viscosity used for particle drag
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
                if (mesh.topoChanging())
                {
                    talphaPhi1Corr0.clear();
                }

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

            // Runs laminar in this case (turbulenceProperties: laminar);
            // kept so the case can switch to RAS/LES by editing that file only.
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            // Advance all Lagrangian parcels one flow timestep using the
            // just-solved U, rho and muc fields.
            particles.evolve();
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
