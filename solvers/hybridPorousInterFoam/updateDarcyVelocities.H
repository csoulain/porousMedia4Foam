// Calculate Porous Media Velocities
/*
Uwetting = Solid*M1*(
         - fvc::grad(p)
         + (1.-alpha1)*fvc::grad(pc)
         -  pc*fvc::grad(alpha1)
         +  rho1*g
         );

UnonWetting = Solid*M2*(
         -fvc::grad(p)
         -alpha1*fvc::grad(pc)
         -pc*fvc::grad(alpha1)
         +rho2*g
        );
*/


//phiw  = M1f/Mf*phi;
//phinw = M2f/Mf*phi;

phiw  = Solidf*M1f/Mf*phi;
phinw = Solidf*M2f/Mf*phi;


Uwetting = fvc::reconstruct(phiw);
Uwetting.correctBoundaryConditions();


UnonWetting = fvc::reconstruct(phinw);
UnonWetting.correctBoundaryConditions();


//Uwetting.correctBoundaryConditions();
//UnonWetting.correctBoundaryConditions();
Solidf = fvc::interpolate(Solid);
