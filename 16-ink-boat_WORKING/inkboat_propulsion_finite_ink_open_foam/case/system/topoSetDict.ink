/*--------------------------------*- C++ -*----------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict.ink;
}

// Ink release region just behind the stern, at the surface
actions
(
    { name inkSource; type cellSet; action new; source boxToCell; box (-0.0215 -0.005 -0.0012) (-0.015 0.005 0.001); }
);
