#include "incompressibleThreePhaseMixture.H"
#include "fvCFD.H"

namespace Foam
{

// Reads constant/transportProperties (as this object's own IOdictionary base)
// and the three phase-fraction fields from the current time directory.
incompressibleThreePhaseMixture::incompressibleThreePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    alpha1_
    (
        IOobject
        (
            "alpha.air",
            U.time().timeName(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    alpha2_
    (
        IOobject
        (
            "alpha.oil",
            U.time().timeName(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    alpha3_
    (
        IOobject
        (
            "alpha.water",
            U.time().timeName(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    rho1_
    (
        dimensionedScalar
        (
            "rho1",
            dimDensity,
            readScalar(lookup("rho1"))
        )
    ),
    rho2_
    (
        dimensionedScalar
        (
            "rho2",
            dimDensity,
            readScalar(lookup("rho2"))
        )
    ),
    rho3_
    (
        dimensionedScalar
        (
            "rho3",
            dimDensity,
            readScalar(lookup("rho3"))
        )
    ),
    mu1_
    (
        dimensionedScalar
        (
            "mu1",
            dimDynamicViscosity,
            readScalar(lookup("mu1"))
        )
    ),
    mu2_
    (
        dimensionedScalar
        (
            "mu2",
            dimDynamicViscosity,
            readScalar(lookup("mu2"))
        )
    ),
    mu3_
    (
        dimensionedScalar
        (
            "mu3",
            dimDynamicViscosity,
            readScalar(lookup("mu3"))
        )
    )
{}

// Linear alpha-weighted mixture viscosity (no sharpening/limiting) — simple
// arithmetic mean, not a harmonic mean, consistent with stock OpenFOAM's
// two-phase mixture classes.
tmp<volScalarField> incompressibleThreePhaseMixture::mu() const
{
    return volScalarField::New
    (
        "mu",
        alpha1_*mu1_ + alpha2_*mu2_ + alpha3_*mu3_
    );
}

// Indicator field: 1 in cells that contain an interface of ANY of the three
// phases, 0 in bulk cells. Three-phase generalisation of stock
// incompressibleTwoPhaseMixture::nearInterface(), which is required by
// $FOAM_SOLVERS/multiphase/VoF/alphaCourantNo.H to restrict the interface
// Courant number to interface cells only.
tmp<volScalarField> incompressibleThreePhaseMixture::nearInterface() const
{
    return volScalarField::New
    (
        "nearInterface",
        max
        (
            max
            (
                pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_),
                pos0(alpha2_ - 0.01)*pos0(0.99 - alpha2_)
            ),
            pos0(alpha3_ - 0.01)*pos0(0.99 - alpha3_)
        )
    );
}


void incompressibleThreePhaseMixture::correct()
{
    // No additional correction needed for now
}

} // End namespace Foam
