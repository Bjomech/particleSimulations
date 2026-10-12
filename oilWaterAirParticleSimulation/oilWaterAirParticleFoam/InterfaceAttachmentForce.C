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
  InterfaceAttachmentForce.C  (NEW)

  Implementation. Structure (constructor reading its own coeffs subdict,
  cacheFields storing derived fields + building interpolators once per time
  step, calcCoupled interpolating at the particle's barycentric position)
  is modelled directly on the real PressureGradientForce<CloudType>, see
  the citations in InterfaceAttachmentForce.H. Physical model: see the
  header comment in InterfaceAttachmentForce.H.
\*---------------------------------------------------------------------------*/

#include "InterfaceAttachmentForce.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

template<class CloudType>
Foam::InterfaceAttachmentForce<CloudType>::InterfaceAttachmentForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& forceType
)
:
    ParticleForce<CloudType>(owner, mesh, dict, forceType, true),
    alphaAirName_
    (
        this->coeffs().template lookupOrDefault<word>("alphaAir", "alpha.air")
    ),
    alphaOilName_
    (
        this->coeffs().template lookupOrDefault<word>("alphaOil", "alpha.oil")
    ),
    c_(this->coeffs().template lookupOrDefault<scalar>("c", 1.0)),
    sigma13_(this->coeffs().template lookupOrDefault<scalar>("sigma13", 0.07)),
    sigma23_(this->coeffs().template lookupOrDefault<scalar>("sigma23", 0.05)),
    deltaN_(1e-8/Foam::pow(Foam::average(mesh.V()).value(), 1.0/3.0)),
    alphaAirInterpPtr_(nullptr),
    gradAlphaAirInterpPtr_(nullptr),
    alphaOilInterpPtr_(nullptr),
    gradAlphaOilInterpPtr_(nullptr)
{}


template<class CloudType>
Foam::InterfaceAttachmentForce<CloudType>::InterfaceAttachmentForce
(
    const InterfaceAttachmentForce& iaf
)
:
    ParticleForce<CloudType>(iaf),
    alphaAirName_(iaf.alphaAirName_),
    alphaOilName_(iaf.alphaOilName_),
    c_(iaf.c_),
    sigma13_(iaf.sigma13_),
    sigma23_(iaf.sigma23_),
    deltaN_(iaf.deltaN_),
    alphaAirInterpPtr_(nullptr),
    gradAlphaAirInterpPtr_(nullptr),
    alphaOilInterpPtr_(nullptr),
    gradAlphaOilInterpPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InterfaceAttachmentForce<CloudType>::~InterfaceAttachmentForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::InterfaceAttachmentForce<CloudType>::cacheFields(const bool store)
{
    static word gradAirName("interfaceAttachment:gradAlphaAir");
    static word gradOilName("interfaceAttachment:gradAlphaOil");

    const bool gradAirExists =
        this->mesh().template foundObject<volVectorField>(gradAirName);
    const bool gradOilExists =
        this->mesh().template foundObject<volVectorField>(gradOilName);

    if (store)
    {
        const volScalarField& alphaAir =
            this->mesh().template lookupObject<volScalarField>(alphaAirName_);

        const volScalarField& alphaOil =
            this->mesh().template lookupObject<volScalarField>(alphaOilName_);

        if (!gradAirExists)
        {
            volVectorField* gradAlphaAirPtr =
                new volVectorField(gradAirName, fvc::grad(alphaAir));

            gradAlphaAirPtr->store();
        }

        if (!gradOilExists)
        {
            volVectorField* gradAlphaOilPtr =
                new volVectorField(gradOilName, fvc::grad(alphaOil));

            gradAlphaOilPtr->store();
        }

        const volVectorField& gradAlphaAir =
            this->mesh().template lookupObject<volVectorField>(gradAirName);

        const volVectorField& gradAlphaOil =
            this->mesh().template lookupObject<volVectorField>(gradOilName);

        alphaAirInterpPtr_.reset
        (
            interpolation<scalar>::New
            (
                this->owner().solution().interpolationSchemes(),
                alphaAir
            ).ptr()
        );

        gradAlphaAirInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                gradAlphaAir
            ).ptr()
        );

        alphaOilInterpPtr_.reset
        (
            interpolation<scalar>::New
            (
                this->owner().solution().interpolationSchemes(),
                alphaOil
            ).ptr()
        );

        gradAlphaOilInterpPtr_.reset
        (
            interpolation<vector>::New
            (
                this->owner().solution().interpolationSchemes(),
                gradAlphaOil
            ).ptr()
        );
    }
    else
    {
        alphaAirInterpPtr_.clear();
        gradAlphaAirInterpPtr_.clear();
        alphaOilInterpPtr_.clear();
        gradAlphaOilInterpPtr_.clear();

        if (gradAirExists)
        {
            const volVectorField& gradAlphaAir =
                this->mesh().template lookupObject<volVectorField>(gradAirName);

            const_cast<volVectorField&>(gradAlphaAir).checkOut();
        }

        if (gradOilExists)
        {
            const volVectorField& gradAlphaOil =
                this->mesh().template lookupObject<volVectorField>(gradOilName);

            const_cast<volVectorField&>(gradAlphaOil).checkOut();
        }
    }
}


template<class CloudType>
Foam::vector Foam::InterfaceAttachmentForce<CloudType>::interfaceForce
(
    const scalar sigma,
    const scalar dp,
    const scalar alpha,
    const vector& gradAlpha
) const
{
    // Unit interface normal, regularised exactly as
    // threePhaseInterfaceProperties::nHatfv() -> 0 smoothly in bulk phases
    const vector nHat = gradAlpha/(mag(gradAlpha) + deltaN_);

    // Localisation weight: 0 away from the interface, 1 at alpha = 0.5
    const scalar L = max(min(4.0*alpha*(1.0 - alpha), 1.0), 0.0);

    // Two-sided restoring/pinning force about alpha = 0.5, see class header
    return -c_*sigma*dp*L*(alpha - 0.5)*nHat;
}


template<class CloudType>
Foam::forceSuSp Foam::InterfaceAttachmentForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);

    const scalar alphaAir =
        alphaAirInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices()
        );

    const vector gradAlphaAir =
        gradAlphaAirInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices()
        );

    const scalar alphaOil =
        alphaOilInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices()
        );

    const vector gradAlphaOil =
        gradAlphaOilInterp().interpolate
        (
            p.coordinates(),
            p.currentTetIndices()
        );

    const scalar dp = p.d();

    value.Su() =
        interfaceForce(sigma13_, dp, alphaAir, gradAlphaAir)
      + interfaceForce(sigma23_, dp, alphaOil, gradAlphaOil);

    return value;
}

// ************************************************************************* //
