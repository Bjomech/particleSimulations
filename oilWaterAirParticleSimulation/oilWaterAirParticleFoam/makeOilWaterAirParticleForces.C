/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    New file written for this project. Uses the OpenFOAM lagrangian/intermediate
    ParticleForce framework and follows the structure of the OpenFOAM
    PressureGradientForce<CloudType> class (Copyright OpenFOAM Foundation /
    OpenCFD Ltd.); see the file description for details.

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
  makeOilWaterAirParticleForces.C  (NEW)

  Registers InterfaceAttachmentForce for basicKinematicCloud so it becomes
  selectable in constant/kinematicCloudProperties's particleForces dict,
  alongside the stock sphereDrag/gravity/pressureGradient/virtualMass
  entries.

  IMPORTANT - only the *Type* macro is called here, never the base macro:

      makeParticleForceModelType(InterfaceAttachmentForce, basicKinematicCloud);   // done here
      makeParticleForceModel(basicKinematicCloud);                                // NOT done here

  makeParticleForceModel(CloudType) defines the run-time-selection TABLE
  ITSELF for ParticleForce<basicKinematicCloud> (defineTemplateTypeName
  AndDebug + defineTemplateRunTimeSelectionTable). liblagrangianIntermediate
  (already linked via Make/options: -llagrangianIntermediate) already calls
  this exactly once when it builds the stock forces for basicKinematicCloud.
  Calling it again from this solver's own object file would attempt to
  redefine those same symbols -> multiple-definition / RTS-table clashes at
  link time, exactly the failure mode reported by users who tried to
  recompile the whole intermediate library as a separate user library
  (https://www.cfd-online.com/Forums/openfoam-programming-development/235716-how-add-particle-forces-user-lib.html,
  https://www.cfd-online.com/Forums/openfoam-programming-development/156573-adding-new-drag-force-icouncoupledkinematicparcelfoam.html).

  makeParticleForceModelType(SS, CloudType) only INSERTS a new entry
  (InterfaceAttachmentForce) into the EXISTING table via
  adddictionaryConstructorToTable - safe to call from a brand-new
  translation unit for a brand-new type name, and exactly the supported
  mechanism for adding a custom ParticleForce from a solver's own Make/files
  without touching or recompiling liblagrangianIntermediate. Macro bodies
  verified against:
  https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-2.2.x/master/src/lagrangian/intermediate/submodels/Kinematic/ParticleForces/ParticleForce/ParticleForce.H
  https://cpp.openfoam.org/v13/ParticleForce_8H_source.html

  This file is listed directly in Make/files (like oilWaterAirParticleFoam.C
  itself); including InterfaceAttachmentForce.H here triggers the
  NoRepository compile-via-inclusion of InterfaceAttachmentForce.C for
  CloudType = basicKinematicCloud, instantiating the template exactly once.
\*---------------------------------------------------------------------------*/

#include "basicKinematicCloud.H"
#include "InterfaceAttachmentForce.H"

namespace Foam
{
    makeParticleForceModelType(InterfaceAttachmentForce, basicKinematicCloud);
}

// ************************************************************************* //
