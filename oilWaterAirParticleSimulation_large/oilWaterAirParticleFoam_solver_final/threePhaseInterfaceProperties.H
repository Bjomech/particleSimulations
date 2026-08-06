#ifndef threePhaseInterfaceProperties_H
#define threePhaseInterfaceProperties_H

#include "incompressibleThreePhaseMixture.H"

namespace Foam
{

class threePhaseInterfaceProperties
{
private:
    const incompressibleThreePhaseMixture& mixture_;

    dimensionedScalar sigma12_;
    dimensionedScalar sigma13_;
    dimensionedScalar sigma23_;

    volScalarField K1_;
    volScalarField K2_;
    volScalarField K3_;

public:
    threePhaseInterfaceProperties(const incompressibleThreePhaseMixture& mixture);

    tmp<surfaceScalarField> surfaceTensionForce() const;
    void correct();
};

} // End namespace Foam

#endif
