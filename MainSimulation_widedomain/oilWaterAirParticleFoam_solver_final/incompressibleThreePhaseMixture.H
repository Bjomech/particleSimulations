#ifndef incompressibleThreePhaseMixture_H
#define incompressibleThreePhaseMixture_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"

// incompressibleThreePhaseMixture
//
// Description
//     Holds the three incompressible, immiscible phase fractions
//     (alpha1=air, alpha2=oil, alpha3=water — see the .C constructor for the
//     field-name mapping) plus their constant densities/viscosities, read
//     from constant/transportProperties (rho1/rho2/rho3, mu1/mu2/mu3). Also
//     an IOdictionary itself, so "mixture_.lookup(...)" (used by
//     threePhaseInterfaceProperties for sigma12/13/23) reads directly from
//     transportProperties. mu() below returns the linear alpha-weighted
//     mixture viscosity used by UEqn.H in place of a turbulence model.

namespace Foam
{

class incompressibleThreePhaseMixture
:
    public IOdictionary
{
protected:
    volScalarField alpha1_;   // "alpha.air"
    volScalarField alpha2_;   // "alpha.oil"
    volScalarField alpha3_;   // "alpha.water"

    dimensionedScalar rho1_;   // air density
    dimensionedScalar rho2_;   // oil density
    dimensionedScalar rho3_;   // water density

    dimensionedScalar mu1_;   // air dynamic viscosity
    dimensionedScalar mu2_;   // oil dynamic viscosity
    dimensionedScalar mu3_;   // water dynamic viscosity

public:
    // Constructor — reads constant/transportProperties and 0/alpha.{air,oil,water}
    incompressibleThreePhaseMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    // Access to phase fractions
    volScalarField& alpha1() { return alpha1_; }
    volScalarField& alpha2() { return alpha2_; }
    volScalarField& alpha3() { return alpha3_; }
    const volScalarField& alpha1() const { return alpha1_; }
    const volScalarField& alpha2() const { return alpha2_; }
    const volScalarField& alpha3() const { return alpha3_; }

    // Densities
    const dimensionedScalar& rho1() const { return rho1_; }
    const dimensionedScalar& rho2() const { return rho2_; }
    const dimensionedScalar& rho3() const { return rho3_; }

    // Viscosities
    const dimensionedScalar& mu1() const { return mu1_; }
    const dimensionedScalar& mu2() const { return mu2_; }
    const dimensionedScalar& mu3() const { return mu3_; }

    // Mixture viscosity (weighted by phase fraction): mu = alpha1*mu1 +
    // alpha2*mu2 + alpha3*mu3. Used directly in UEqn.H's viscous term since
    // no turbulence model is constructed in this solver (see createFields.H).
    tmp<volScalarField> mu() const;

    // Currently a no-op placeholder (densities/viscosities are constant in
    // this incompressible mixture) but kept so the class satisfies the same
    // "correct()" interface other mixture classes (e.g. two-phase interFoam's)
    // expose — called once per outer iteration from immiscibleIncompressible...::correct().
    void correct();
};

} // End namespace Foam

#endif
