#include "threePhaseInterfaceProperties.H"
#include "fvCFD.H"
#include "fvc.H"

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
        dimensionSet(1, 0, -2, 0, 0),
        mixture_.get<scalar>("sigma12")
    ),
    sigma13_
    (
        "sigma13",
        dimensionSet(1, 0, -2, 0, 0),
        mixture_.get<scalar>("sigma13")
    ),
    sigma23_
    (
        "sigma23",
        dimensionSet(1, 0, -2, 0, 0),
        mixture_.get<scalar>("sigma23")
    ),

    // sigma_i + sigma_j = sigma_ij   =>   pairwise-consistent CSF coefficients
    sigma1_(0.5*(sigma12_ + sigma13_ - sigma23_)),
    sigma2_(0.5*(sigma12_ + sigma23_ - sigma13_)),
    sigma3_(0.5*(sigma13_ + sigma23_ - sigma12_)),

    // Mesh-size-scaled normal stabiliser (identical to stock interfaceProperties)
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mixture_.alpha1().mesh().V()), 1.0/3.0)
    ),

    K1_
    (
        IOobject
        (
            "K.air",
            mixture_.alpha1().mesh().time().timeName(),
            mixture_.alpha1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE          // set AUTO_WRITE if you want it in ParaView
        ),
        mixture_.alpha1().mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    K2_
    (
        IOobject
        (
            "K.oil",
            mixture_.alpha1().mesh().time().timeName(),
            mixture_.alpha1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture_.alpha1().mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    K3_
    (
        IOobject
        (
            "K.water",
            mixture_.alpha1().mesh().time().timeName(),
            mixture_.alpha1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture_.alpha1().mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    )
{
    Info<< "Three-phase CSF coefficients:"
        << " sigma.air = "   << sigma1_.value()
        << ", sigma.oil = "  << sigma2_.value()
        << ", sigma.water = "<< sigma3_.value()
        << " [N/m], deltaN = " << deltaN_.value() << " [1/m]" << endl;
}


void threePhaseInterfaceProperties::calculateK
(
    const volScalarField& alpha,
    volScalarField& K
) const
{
    const fvMesh& mesh = alpha.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha, "nHat"));

    // Interpolated face gradient of alpha
    const surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    const surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    // Face unit interface normal flux
    const surfaceScalarField nHatf(nHatfv & Sf);

    // Curvature
    K = -fvc::div(nHatf);
}


tmp<surfaceScalarField>
threePhaseInterfaceProperties::surfaceTensionForce() const
{
    return
        fvc::interpolate(sigma1_*K1_)*fvc::snGrad(mixture_.alpha1())
      + fvc::interpolate(sigma2_*K2_)*fvc::snGrad(mixture_.alpha2())
      + fvc::interpolate(sigma3_*K3_)*fvc::snGrad(mixture_.alpha3());
}


void threePhaseInterfaceProperties::correct()
{
    calculateK(mixture_.alpha1(), K1_);
    calculateK(mixture_.alpha2(), K2_);
    calculateK(mixture_.alpha3(), K3_);
}

} // End namespace Foam
