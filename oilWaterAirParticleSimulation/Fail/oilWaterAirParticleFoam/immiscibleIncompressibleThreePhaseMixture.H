/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/interMixingFoam/immiscibleIncompressibleThreePhaseMixture/immiscibleIncompressibleThreePhaseMixture.H
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

#ifndef immiscibleIncompressibleThreePhaseMixture_H
#define immiscibleIncompressibleThreePhaseMixture_H

#include "incompressibleThreePhaseMixture.H"
#include "threePhaseInterfaceProperties.H"

// immiscibleIncompressibleThreePhaseMixture
//
// Description
//     Thin mixin/facade class combining incompressibleThreePhaseMixture
//     (phase fractions + densities/viscosities) with threePhaseInterfaceProperties
//     (interface curvature + surface tension force) into a single object,
//     mirroring the pattern of stock OpenFOAM's immiscibleIncompressibleTwoPhaseMixture
//     but extended to three phases (air/oil/water). This is the type actually
//     used as "mixture" in createFields.H/the main solver — it forwards calls
//     to whichever base class implements them (see "using" declarations below).

namespace Foam
{

class immiscibleIncompressibleThreePhaseMixture
:
    public incompressibleThreePhaseMixture,
    public threePhaseInterfaceProperties
{
public:
    immiscibleIncompressibleThreePhaseMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    )
    :
        incompressibleThreePhaseMixture(U, phi),
        // threePhaseInterfaceProperties needs a reference to the (already
        // constructed) incompressibleThreePhaseMixture base above, hence the cast.
        threePhaseInterfaceProperties(static_cast<const incompressibleThreePhaseMixture&>(*this))
    {}

    // Expose the phase-fraction/density/viscosity accessors unambiguously
    // (both base classes could otherwise clash via multiple inheritance).
    using incompressibleThreePhaseMixture::alpha1;
    using incompressibleThreePhaseMixture::alpha2;
    using incompressibleThreePhaseMixture::alpha3;
    using incompressibleThreePhaseMixture::rho1;
    using incompressibleThreePhaseMixture::rho2;
    using incompressibleThreePhaseMixture::rho3;
    using incompressibleThreePhaseMixture::mu;

    // Forwards to threePhaseInterfaceProperties (disambiguates which base to call).
    tmp<surfaceScalarField> surfaceTensionForce() const
    {
        return threePhaseInterfaceProperties::surfaceTensionForce();
    }

    // Updates both bases each outer iteration: mixture densities/viscosities
    // (currently a no-op) then interface curvature K1/K2/K3.
    void correct()
    {
        incompressibleThreePhaseMixture::correct();
        threePhaseInterfaceProperties::correct();
    }
};

} // End namespace Foam

#endif
