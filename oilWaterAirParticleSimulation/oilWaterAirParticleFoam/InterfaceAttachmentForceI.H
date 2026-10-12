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
  InterfaceAttachmentForceI.H  (NEW)

  Inline access functions for InterfaceAttachmentForce, following exactly
  the same pattern as the real PressureGradientForceI.H::DUcDtInterp()
  (fatal error if the interpolator has not been cached yet, i.e.
  cacheFields(true) was not called before calcCoupled()). Verified against
  https://cpp.openfoam.org/v13/PressureGradientForceI_8H_source.html
\*---------------------------------------------------------------------------*/

template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::InterfaceAttachmentForce<CloudType>::alphaAirInterp() const
{
    if (!alphaAirInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "alpha.air interpolation object not set"
            << abort(FatalError);
    }

    return alphaAirInterpPtr_();
}


template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::InterfaceAttachmentForce<CloudType>::gradAlphaAirInterp() const
{
    if (!gradAlphaAirInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "grad(alpha.air) interpolation object not set"
            << abort(FatalError);
    }

    return gradAlphaAirInterpPtr_();
}


template<class CloudType>
inline const Foam::interpolation<Foam::scalar>&
Foam::InterfaceAttachmentForce<CloudType>::alphaOilInterp() const
{
    if (!alphaOilInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "alpha.oil interpolation object not set"
            << abort(FatalError);
    }

    return alphaOilInterpPtr_();
}


template<class CloudType>
inline const Foam::interpolation<Foam::vector>&
Foam::InterfaceAttachmentForce<CloudType>::gradAlphaOilInterp() const
{
    if (!gradAlphaOilInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "grad(alpha.oil) interpolation object not set"
            << abort(FatalError);
    }

    return gradAlphaOilInterpPtr_();
}

// ************************************************************************* //
