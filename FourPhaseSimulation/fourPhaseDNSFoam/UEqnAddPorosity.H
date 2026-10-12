/*---------------------------------------------------------------------------*\
  This file is part of particleSimulations (Bjomech/particleSimulations).
-------------------------------------------------------------------------------
Origin
    Derived from OpenFOAM (www.openfoam.com), upstream reference:
      applications/solvers/multiphase/interFoam/ and interFoam/interMixingFoam/ (solver family)
    Copyright (C) OpenFOAM Foundation and OpenCFD Ltd.; see the upstream
    file for the exact original copyright years.
    Developed and ported mainly against OpenFOAM v2512; the exact upstream
    version originally copied has not been verified.

Modifications
    Copyright (C) 2026 Bjoern Klepsch
    Last modified: 2026-07-20 (see git history for the detailed change log)
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

// Including porosity effects in UEqn following:
// Jensen, B., Jacobsen, N. G., & Christensen, E. D. (2014).
// Investigations on the porous media equations and resistance
// coefficients for coastal structures. Coastal Engineering, 84, 56-72.

if (porosityEnabled)
{
    const volScalarField& porosity = tporosity.cref();

    const word porosityModel("JensenEtAl2014");
    const dictionary& dict =
        porosityProperties.subDict(porosityModel + "Coeffs");
    const dimensionedScalar alpha(dimless/dimArea, dict.get<scalar>("alpha"));
    const dimensionedScalar beta(dimless/dimLength, dict.get<scalar>("beta"));
    const dimensionedScalar d50(dimless, dict.get<scalar>("d50"));
    const dimensionedScalar KC(dimless, dict.get<scalar>("KC"));

    // Generating Darcy-Forchheimer coefficient: F = rho*U*(a + b*|U|)
    // Shoud it be mu or muEff in the equation below?
    {
        // Darcy term
        volScalarField DarcyForchheimerCoeff
        (
            alpha*sqr(1 - porosity)*mixture.mu()/sqr(porosity)/sqr(d50)
        );

        // Adding Forchheimer term
        DarcyForchheimerCoeff += rho*mag(U)
            *beta*(1 + pos(KC)*7.5/KC)*(1 - porosity)/sqr(porosity)/d50;

        // Adding Darcy-Forchheimer term as implicit source term
        UEqn += fvm::Sp(DarcyForchheimerCoeff, U);
    }

    {
        // Generating added mass force coefficient
        const dimensionedScalar gamma_p(dimless, dict.get<scalar>("gamma_p"));
        const volScalarField Cm(gamma_p*(1 - porosity));

        UEqn += Cm*fvm::ddt(rho, U);
        UEqn += Cm*MRF.DDt(rho, U);
    }

    // Dividing both matrix entries and source term by porosity to compensate
    // for the fact that the FVM cell volume averages use division by cell
    // volume V whereas only the cell pore volume, porosity*V, is accessible.
    UEqn *= scalar(1)/porosity;
}
