#ifndef immiscibleIncompressibleThreePhaseMixture_H
#define immiscibleIncompressibleThreePhaseMixture_H

#include "incompressibleThreePhaseMixture.H"
#include "threePhaseInterfaceProperties.H"

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
        threePhaseInterfaceProperties(static_cast<const incompressibleThreePhaseMixture&>(*this))
    {}

    using incompressibleThreePhaseMixture::alpha1;
    using incompressibleThreePhaseMixture::alpha2;
    using incompressibleThreePhaseMixture::alpha3;
    using incompressibleThreePhaseMixture::rho1;
    using incompressibleThreePhaseMixture::rho2;
    using incompressibleThreePhaseMixture::rho3;
    using incompressibleThreePhaseMixture::mu;

    tmp<surfaceScalarField> surfaceTensionForce() const
    {
        return threePhaseInterfaceProperties::surfaceTensionForce();
    }

    void correct()
    {
        incompressibleThreePhaseMixture::correct();
        threePhaseInterfaceProperties::correct();
    }
};

} // End namespace Foam

#endif
