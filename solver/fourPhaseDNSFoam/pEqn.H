tmp<volScalarField> trAU(1.0/UEqn.A());
const volScalarField& rAU = trAU();

surfaceScalarField rAUf
(
    "rAUf",
    fvc::interpolate(rAU)
);

volVectorField HbyA
(
    constrainHbyA
    (
        rAU*UEqn.H(),
        U,
        p_rgh
    )
);

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi))
);

MRF.makeRelative(phiHbyA);

if (p_rgh.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p_rgh);
    fvc::makeAbsolute(phiHbyA, U);
}

surfaceScalarField phig
(
    (
        mixture.surfaceTensionForce()
      - ghf*fvc::snGrad(rho)
    )*rAUf*mesh.magSf()
);
// NOTE: sign of the ghf*snGrad(rho) buoyancy term must match UEqn.H's
// momentum predictor (also "-ghf*fvc::snGrad(rho)") and the standard
// OpenFOAM interFoam/interMixingFoam convention. A "+" here previously
// reversed the buoyancy-driven pressure-correction flux, which corrupts
// the carrier-fluid rise/fall behaviour that the particle/bubble clouds
// are advected/dragged by.

phiHbyA += phig;

constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference
    (
        pRefCell,
        p_rgh[pRefCell]
    );

    p_rghEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - p_rghEqn.flux();

        p_rgh.relax();

        U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
        U.correctBoundaryConditions();

        fvOptions.correct(U);
    }
}

#include "continuityErrs.H"

fvc::makeRelative(phi, U);

p == p_rgh + rho*gh;

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - p[pRefCell]
    );
}

p_rgh == p - rho*gh;
