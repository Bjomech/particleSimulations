/*---------------------------------------------------------------------------*\
  fourPhaseDNSFoam — 4-phase DNS solver
  3-phase Eulerian VoF (interMixingFoam base) + Lagrangian particles/bubbles
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleThreePhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "basicKinematicCloud.H"
#include "attachmentModel.H"

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting 4-phase DNS simulation\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "alphaControls.H"
        #include "alphaEqnSubCycle.H"

        mixture.correct();

        particleCloud.evolve();
        bubbleCloud.evolve();

        attachment.solve(particleCloud, bubbleCloud, runTime.deltaT().value());

        #include "UEqn.H"

        while (pimple.loop())
        {
            #include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n" << endl;
    }

    return 0;
}
