/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/interMixingFoam/interMixingFoam.C
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
