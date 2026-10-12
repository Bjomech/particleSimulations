/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    New file written for this project. Uses the OpenFOAM lagrangian/intermediate
    ParticleForce framework and follows the structure of the OpenFOAM
    PressureGradientForce<CloudType> class (Copyright OpenFOAM Foundation /
    OpenCFD Ltd.); see the file description for details.

Modifications
    Copyright (C) 2026 Bjoern Klepsch
    Last modified: 2026-09-08 (see git history for the detailed change log)
-------------------------------------------------------------------------------
License
    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See the LICENSE file in the repository root.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  InterfaceAttachmentForce  (NEW - custom, non-stock ParticleForce)

  Description
  ────────────
  Sub-grid capillary-attachment / interface-pinning force for Lagrangian
  point particles in the three-phase (air/oil/water) VOF field of this
  solver. NOT a stock OpenFOAM force model - it is written from scratch for
  this thesis case, following the structural pattern of the real, stock
  PressureGradientForce<CloudType> (verified against actual OpenFOAM
  source, see citations below) and registered into the SAME run-time
  selection table stock forces (sphereDrag, gravity, pressureGradient,
  virtualMass) already use, so it is selectable in kinematicCloudProperties
  exactly like them:

      subModels { particleForces { interfaceAttachment { c 1.0; } } }

  Why this force exists (motivation, see also constant/kinematicCloudProperties)
  ──────────────────────────────────────────────────────────────────────────────
  At the particle size forced by the point-particle constraint on this mesh
  (d_p = 14 um = 1 cell, rho0 = 2500 kg/m3) the Stokes number is
  St = tau_p*U_b/d_b ~ 5e-3: the particles are near-passive tracers. Drag,
  gravity, pressureGradient and virtualMass (the stock forces already
  enabled) cannot, by themselves, make a particle attach to and travel with
  a bubble/droplet interface - a tracer just follows the local flow and
  crosses interfaces freely. Real particle-laden multiphase flows show
  particles CAN be captured and held at a fluid-fluid interface, because
  crossing the interface costs surface energy (this is the same physics
  behind flotation, Pickering emulsions, etc.) - i.e. there is a genuine
  physical capillary "pinning" effect the stock point-particle drag/gravity
  forces do not capture at all. This class is a deliberately simple,
  heuristic sub-grid MODEL of that pinning effect, not a rigorously derived
  formula - that is stated here explicitly so it is not mistaken for a
  validated stock capability.

  Model
  ─────
  For each of the two phase interfaces the particle can be pinned to - the
  air/water interface (alpha.air, tension sigma13) and the oil/water
  interface (alpha.oil, tension sigma23) - define at the particle position:

      nHat_i  = grad(alpha_i) / ( |grad(alpha_i)| + deltaN )     [unit vector]
      L_i     = clamp( 4*alpha_i*(1 - alpha_i), 0, 1 )           [0 in bulk, 1 at alpha_i=0.5]

  and the force contribution from interface i:

      F_i = -c * sigma_i * d_p * L_i * (alpha_i - 0.5) * nHat_i

  Total force: F = F_air + F_oil (forceSuSp::Su(), fully explicit, Sp = 0).

  Why this form (dimensions + physical behaviour):
    - sigma_i [N/m] * d_p [m] = [N], already a force; L_i, (alpha_i-0.5) and
      nHat_i are all dimensionless, so F is correctly in Newtons - unlike a
      naive F ~ sigma*d_p*grad(alpha) (grad(alpha) alone has units 1/m, one
      length short of a force). This refines the rough F ∝ -c*sigma*d_p*
      grad(alpha_air) sketch noted earlier in this solver's README to be
      dimensionally correct while keeping the same physical ingredients.
    - (alpha_i - 0.5) changes sign exactly at the interface: a particle
      sitting deeper in phase i (alpha_i > 0.5) is pushed back OUT along
      -nHat_i (toward the interface); a particle on the other side
      (alpha_i < 0.5) is pushed IN along +nHat_i (also toward the
      interface). This is a genuine two-sided restoring/pinning force
      around alpha_i = 0.5, not a one-directional push - the defining
      property of "attachment".
    - L_i and nHat_i both vanish smoothly away from the interface (nHat_i
      because grad(alpha_i) -> 0 in bulk phases, regularised by deltaN
      exactly as threePhaseInterfaceProperties::nHatfv already does, so it
      -> 0, not NaN), so the force is genuinely localised to the interface
      region (a few cells), consistent with it being a SUB-GRID model of an
      interfacial phenomenon rather than a bulk body force.
    - c is a dimensionless, user-tunable strength; c = O(1) makes the force
      comparable in scale to sigma_i*d_p (an actual capillary force scale
      for a particle of this size), which is deliberately large relative to
      this particle's tiny inertia/drag/gravity so it can visibly override
      the tracer-like behaviour noted above. Tune c down if it makes
      particles stick unrealistically hard; tune it up if St ~ 5e-3 still
      dominates and no attachment is visible.

  sigma13 (air-water) and sigma23 (oil-water) default to 0.07 and 0.05 N/m
  respectively - the SAME physical pairwise tensions already implied by
  threePhaseInterfaceProperties's per-phase CSF coefficients (sigma.air =
  0.025, sigma.oil = 0.005, sigma.water = 0.045 N/m in constant/
  transportProperties, via sigma13 = sigma.air + sigma.water,
  sigma23 = sigma.oil + sigma.water). They are re-specified here (not
  looked up from the mixture object) because a ParticleForce is constructed
  by dictionary-based run-time selection and can only look up objects
  registered in the mesh database; threePhaseInterfaceProperties is a plain
  local object in oilWaterAirParticleFoam.C, not registered, so it cannot be
  referenced from here directly (see the NOTE next to
  threePhaseInterfaceProperties::sigma13()/sigma23() in
  threePhaseInterfaceProperties.H). Keep the two values consistent by hand
  if constant/transportProperties is ever changed.

  Structural reference / citations (verified against real OpenFOAM source)
  ──────────────────────────────────────────────────────────────────────────
    PRIMARY: every signature below was checked against the actual ESI
    OpenFOAM-v2512 source tree (https://dl.openfoam.com/source/v2512/),
    src/lagrangian/intermediate/submodels/Kinematic/ParticleForces/:
      - ParticleForce/ParticleForce.H: 5-arg base constructor
        (owner, mesh, dict, forceType, readCoeffs); calcCoupled/
        calcNonCoupled DO take a trackingData parameter in ESI too;
        RTS parameter list (owner, mesh, dict); macro bodies.
      - PressureGradient/PressureGradientForce.C and
        Interface/InterfaceForce.C: cacheFields getObjectPtr/store()/
        checkOut() pattern; interpolation at the parcel position is
        interp.interpolate(p.coordinates(), p.currentTetIndices()) -
        currentTetIndices() takes NO argument in ESI (the particle holds
        its own mesh reference; the Foundation fork passes td.mesh instead -
        using that here was a real v2512 compile error, since fixed).
      - NB: ESI v2512 ships a stock single-interface 'InterfaceForce'
        (coeffs C + alpha, non-coupled, Su = C*m*grad(alpha*(1-alpha)));
        this class is its dimensional, two-interface, coupled analogue.
    Secondary (design provenance; Foundation v13 pages used for the first
    draft - correct for the overall structure, NOT for ESI signatures):
      https://cpp.openfoam.org/v13/ParticleForce_8H_source.html
      https://cpp.openfoam.org/v13/PressureGradientForce_8C_source.html
    - makeParticleForceModel / makeParticleForceModelType macro bodies (why
      only the *Type* macro is invoked in makeOilWaterAirParticleForces.C,
      not the base macro - see the comment there):
      https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-2.2.x/master/src/lagrangian/intermediate/submodels/Kinematic/ParticleForces/ParticleForce/ParticleForce.H
    - forceSuSp container (Su/Sp semantics, F = Sp*(U-Up) + Su):
      https://cpp.openfoam.org/v13/forceSuSp_8H_source.html
    - Extensibility pattern for adding a new ParticleForce type without
      recompiling liblagrangianIntermediate:
      https://www.cfd-online.com/Forums/openfoam-programming-development/235716-how-add-particle-forces-user-lib.html

  SourceFiles
      InterfaceAttachmentForce.C
\*---------------------------------------------------------------------------*/

#ifndef InterfaceAttachmentForce_H
#define InterfaceAttachmentForce_H

#include "ParticleForce.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class InterfaceAttachmentForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InterfaceAttachmentForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected data

        //- Name of the air volume-fraction field (default "alpha.air")
        const word alphaAirName_;

        //- Name of the oil volume-fraction field (default "alpha.oil")
        const word alphaOilName_;

        //- Dimensionless strength coefficient c (default 1.0)
        const scalar c_;

        //- Air-water surface tension sigma13 [N/m] (default 0.07)
        const scalar sigma13_;

        //- Oil-water surface tension sigma23 [N/m] (default 0.05)
        const scalar sigma23_;

        //- Regularisation length scale [1/m] for the unit interface
        //  normal, same formula/role as
        //  threePhaseInterfaceProperties::deltaN_: 1e-8/(avg cell size)
        const scalar deltaN_;

        //- Interpolator for alpha.air at the particle position
        autoPtr<interpolation<scalar>> alphaAirInterpPtr_;

        //- Interpolator for grad(alpha.air) at the particle position
        autoPtr<interpolation<vector>> gradAlphaAirInterpPtr_;

        //- Interpolator for alpha.oil at the particle position
        autoPtr<interpolation<scalar>> alphaOilInterpPtr_;

        //- Interpolator for grad(alpha.oil) at the particle position
        autoPtr<interpolation<vector>> gradAlphaOilInterpPtr_;


    // Protected Member Functions

        //- Single-interface contribution:
        //  -c*sigma*dp*L(alpha)*(alpha-0.5)*nHat, see class header
        vector interfaceForce
        (
            const scalar sigma,
            const scalar dp,
            const scalar alpha,
            const vector& gradAlpha
        ) const;


public:

    //- Runtime type information
    TypeName("interfaceAttachment");


    // Constructors

        //- Construct from mesh
        InterfaceAttachmentForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        InterfaceAttachmentForce(const InterfaceAttachmentForce& iaf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new InterfaceAttachmentForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~InterfaceAttachmentForce();


    // Member Functions

        // Access

            //- Return the alpha.air interpolator
            inline const interpolation<scalar>& alphaAirInterp() const;

            //- Return the grad(alpha.air) interpolator
            inline const interpolation<vector>& gradAlphaAirInterp() const;

            //- Return the alpha.oil interpolator
            inline const interpolation<scalar>& alphaOilInterp() const;

            //- Return the grad(alpha.oil) interpolator
            inline const interpolation<vector>& gradAlphaOilInterp() const;


        // Evaluation

            //- Cache carrier fields/interpolators (called once per time
            //  step by KinematicCloud, exactly as for PressureGradientForce)
            virtual void cacheFields(const bool store);

            //- Calculate the coupled force (explicit Su only, Sp = 0)
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InterfaceAttachmentForceI.H"

#ifdef NoRepository
    #include "InterfaceAttachmentForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
