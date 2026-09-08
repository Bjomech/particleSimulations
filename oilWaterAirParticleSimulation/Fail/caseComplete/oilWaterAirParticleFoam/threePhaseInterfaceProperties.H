#ifndef threePhaseInterfaceProperties_H
#define threePhaseInterfaceProperties_H

#include "incompressibleThreePhaseMixture.H"

// threePhaseInterfaceProperties  (FIXED VERSION)
//
// Description
//     Per-phase interface curvature (K.air, K.oil, K.water) and the combined
//     CSF surface-tension force for a three-phase (air/oil/water) VOF mixture.
//
//     Curvature is computed exactly the way stock OpenFOAM's two-phase
//     interfaceProperties::calculateK() does it, i.e. from the *face*
//     interpolated gradient:
//
//         gradAlphaf = interpolate(grad(alpha))
//         nHatfv     = gradAlphaf/(mag(gradAlphaf) + deltaN)
//         nHatf      = nHatfv & Sf
//         K          = -div(nHatf)
//
//     This is what makes the CSF force "well balanced" against the
//     face-based pressure gradient in pEqn.H. Computing K from a smoothed
//     *cell* gradient and then masking it (previous version) breaks that
//     balance and produces a large spurious net force on a closed interface.
//
//     Pairwise sigma handling: with three immiscible phases the CSF force is
//     written as sum_i sigma_i * K_i * grad(alpha_i), where the per-phase
//     coefficients follow from requiring sigma_i + sigma_j = sigma_ij:
//
//         sigma1 = (sigma12 + sigma13 - sigma23)/2      (air)
//         sigma2 = (sigma12 + sigma23 - sigma13)/2      (oil)
//         sigma3 = (sigma13 + sigma23 - sigma12)/2      (water)
//
//     At a pure air/water interface alpha.oil = 0, grad(alpha.air) =
//     -grad(alpha.water) and K.air = -K.water, so the sum collapses to
//     sigma13*K.air*grad(alpha.air) - exactly the two-phase result. The same
//     holds for the oil/water and air/oil interfaces, so sigma12 (air-oil) is
//     now actually used instead of being ignored.

namespace Foam
{

class threePhaseInterfaceProperties
{
private:

    const incompressibleThreePhaseMixture& mixture_;

    dimensionedScalar sigma12_;   // air-oil   surface tension [N/m]
    dimensionedScalar sigma13_;   // air-water surface tension [N/m]
    dimensionedScalar sigma23_;   // oil-water surface tension [N/m]

    // Per-phase CSF coefficients derived from the three pairwise values
    dimensionedScalar sigma1_;
    dimensionedScalar sigma2_;
    dimensionedScalar sigma3_;

    // Stabilisation of the unit normal, scaled with mesh size exactly as in
    // stock interfaceProperties: 1e-8/pow(average(mesh.V()), 1/3)
    dimensionedScalar deltaN_;

    volScalarField K1_;   // curvature of the air   interface   ("K.air")
    volScalarField K2_;   // curvature of the oil   interface   ("K.oil")
    volScalarField K3_;   // curvature of the water interface   ("K.water")

    //- Curvature of one phase, stock-style (face-based)
    void calculateK(const volScalarField& alpha, volScalarField& K) const;

public:

    threePhaseInterfaceProperties
    (
        const incompressibleThreePhaseMixture& mixture
    );

    //- Combined surface-tension force as a face field, used BOTH in the
    //  momentum predictor (UEqn.H) and in phiHbyA (pEqn.H)
    tmp<surfaceScalarField> surfaceTensionForce() const;

    //- Recompute K.air/K.oil/K.water from the current alpha fields
    void correct();

    // Access
    const volScalarField& K1() const { return K1_; }
    const volScalarField& K2() const { return K2_; }
    const volScalarField& K3() const { return K3_; }
};

} // End namespace Foam

#endif
