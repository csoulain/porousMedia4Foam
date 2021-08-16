/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "3600";
    object      Y.H;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
50
(
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
110.7
)
;

boundaryField
{
    outlet
    {
        type            zeroGradient;
    }
    inlet
    {
        type            fixedValue;
        value           uniform 110.7;
    }
    frontAndBack
    {
        type            empty;
    }
}


// ************************************************************************* //
