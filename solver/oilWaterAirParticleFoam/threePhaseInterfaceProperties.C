#include "threePhaseInterfaceProperties.H"
#include "fvc.H"
#include "fvCFD.H"

namespace Foam
{

threePhaseInterfaceProperties::threePhaseInterfaceProperties
(
    const incompressibleThreePhaseMixture& mixture
)
:
    mixture_(mixture),
    sigma12_
    (
        "sigma12",
        dimensionSet(1,0,-2,0,0),
        mixture_.lookup("sigma12")
    ),
    sigma13_
    (
        "sigma13",
        dimensionSet(1,0,-2,0,0),
        mixture_.lookup("sigma13")
    ),
    sigma23_
    (
        "sigma23",
        dimensionSet(1,0,-2,0,0),
        mixture_.lookup("sigma23")
    ),
    K1_
    (
        IOobject
        (
            "K1",
            mixture_.alpha1().time().timeName(),
            mixture_.alpha1().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mixture_.alpha1().mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    K2_(K1_),
    K3_(K1_)
{}

tmp<surfaceScalarField> threePhaseInterfaceProperties::surfaceTensionForce() const
{
    // Each term is (interpolated sigma*curvature) * snGrad(alpha) for that
    // phase, summed to give the total surface-tension pressure-equation source.
    // Coefficient pairing (matches which two phases border each interface):
    //   phase 1 (air)  term uses sigma13_ (air-water) with K1_ (air curvature)
    //   phase 2 (oil)  term uses sigma23_ (oil-water)  with K2_ (oil curvature)
    //   phase 3 (water) term uses sigma23_ (oil-water) with K3_ (water curvature)
    // sigma12_ (air-oil) is not used directly here; its effect enters through
    // how alpha2 is derived (alpha2 = 1 - alpha1 - alpha3) and how the phases
    // are ordered in incompressibleThreePhaseMixture. This is the same
    // sigma/curvature pairing scheme as OpenFOAM's stock threePhaseInterfaceProperties.
    return
        fvc::interpolate(sigma13_*K1_)*fvc::snGrad(mixture_.alpha1())
      + fvc::interpolate(sigma23_*K2_)*fvc::snGrad(mixture_.alpha2())
      + fvc::interpolate(sigma23_*K3_)*fvc::snGrad(mixture_.alpha3());
}

void threePhaseInterfaceProperties::correct()
{
    // Compute interface curvature K = -div(nHat) for each phase.
    // nHat = grad(alpha) / (|grad(alpha)| + delta)  — stabilised unit normal.
    // The small delta prevents division by zero in single-phase cells where
    // grad(alpha) = 0.

    const dimensionedScalar deltaN
    (
        "deltaN",
        dimless/dimLength,
        1e-8
    );

    auto computeK = [&](const volScalarField& alpha) -> volScalarField
    {
        volVectorField gradAlpha(fvc::grad(alpha));

        volScalarField magGradAlpha
        (
            IOobject
            (
                "magGradAlpha",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mag(gradAlpha)
        );

        // Stabilised unit normal pointing out of the phase
        volVectorField nHat
        (
            IOobject
            (
                "nHat",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gradAlpha / (magGradAlpha + deltaN)
        );

        // Curvature = -divergence of unit normal
        volScalarField K
        (
            IOobject
            (
                "K_tmp",
                alpha.time().timeName(),
                alpha.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -fvc::div(nHat)
        );

        return K;
    };

    K1_ = computeK(mixture_.alpha1());
    K2_ = computeK(mixture_.alpha2());
    K3_ = computeK(mixture_.alpha3());
}

} // End namespace Foam
