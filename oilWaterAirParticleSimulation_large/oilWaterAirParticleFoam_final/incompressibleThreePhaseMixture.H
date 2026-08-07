#ifndef incompressibleThreePhaseMixture_H
#define incompressibleThreePhaseMixture_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"

namespace Foam
{

class incompressibleThreePhaseMixture
:
    public IOdictionary
{
protected:
    volScalarField alpha1_;
    volScalarField alpha2_;
    volScalarField alpha3_;

    dimensionedScalar rho1_;
    dimensionedScalar rho2_;
    dimensionedScalar rho3_;

    dimensionedScalar mu1_;
    dimensionedScalar mu2_;
    dimensionedScalar mu3_;

public:
    // Constructor
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

    // Mixture viscosity (weighted by phase fraction)
    tmp<volScalarField> mu() const;

    // Update (does nothing, but needed for interface)
    void correct();
};

} // End namespace Foam

#endif
