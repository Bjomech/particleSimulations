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
