#!/bin/bash

# Create 0 directory if it doesn't exist
mkdir -p 0

# ---------- alpha.air ----------
cat > 0/alpha.air << 'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2512                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0";
    object      alpha.air;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    walls
    {
        type            zeroGradient;
    }
    top
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }
    frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //
EOF

# ---------- alpha.oil ----------
cat > 0/alpha.oil << 'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2512                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0";
    object      alpha.oil;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    walls
    {
        type            zeroGradient;
    }
    top
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }
    frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //
EOF

# ---------- alpha.water ----------
cat > 0/alpha.water << 'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2512                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       volScalarField;
    location    "0";
    object      alpha.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 1;

boundaryField
{
    walls
    {
        type            zeroGradient;
    }
    top
    {
        type            inletOutlet;
        inletValue      uniform 1;
        value           uniform 1;
    }
    frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //
EOF

# ---------- U ----------
cat > 0/U << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    walls        { type fixedValue; value uniform (0 0 0); }
    top          { type pressureInletOutletVelocity; value uniform (0 0 0); }
    frontAndBack { type empty; }
}
EOF

# ---------- p_rgh ----------
cat > 0/p_rgh << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    walls
    {
        type            fixedFluxPressure;
        value           uniform 0;
    }
    top
    {
        type            totalPressure;
        p0              uniform 0;
        value           uniform 0;
    }
    frontAndBack
    {
        type            empty;
    }
}
EOF

echo "0 folder successfully recreated with all required files."
