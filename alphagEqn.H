/*              COMPUTE THE FREE SURFACE DETECTION FUNCTION                   */

//Recompute the gradient of alpha1
gradAlpha1 = fvc::grad(alpha1);


/*           COMPUTE THE SOURCE TERM IN THE ALPHAG EQUATION                   */

// Interface-normal unit vectors
//volVectorField nHat(gradAlpha1/(mag(gradAlpha1) + mixture.deltaN()));

//scalar interfaceAngle = 90 - spillwayAngle.value();
//interfaceAngle = interfaceAngle/180*3.1415926535;

//volVectorField nHat(vector(-Foam::cos(interfaceAngle), -Foam::sin(interfaceAngle), 0));

// Interface-normal velocity derivative, keep positive values only
//
//


UnHat = U & interfaceNormal;
dUdn = Foam::posPart(fvc::grad(UnHat) & interfaceNormal);



// Compute the surface roughness
a = 0.5*Foam::pow(0.09, 0.75)
                   *Foam::pow(turbulence->k(), 1.5)
                   /turbulence->epsilon();

// Compute the source term for the alphag equation
alphagSource = a / phiEnt * dUdn * deltaFreeSurface;

scalar dbba = 0.1;

// Bound the source term in space using appropriate criteria
alphagSource *= pos(1.0 - alphag);
alphagSource *= pos(alpha1 - dbba);
alphagSource *= pos(UnHat);


// Info << "kCritical" << kCritical << endl;
// if (kCritical.value() > 0)
// {
//     alphagSource *= pos(turbulence->k() - kCritical);
// }
// else
// {
//     Pt = turbulence->k()*rho;
//     Pd = rho*mag(g & nHat)*a + sigma/a;  
//     alphagSource *= pos(Pt - Pd);
// }

Info<< "Air entrainment source: "
    << "  Min(alphagSource) = " << min(alphagSource).value()
    << "  Max(alphagSource) = " << max(alphagSource).value()
    << endl;


/*              ASSEMBLE AND SOLVE THE ALPHAG EQUATION                        */

scalarField & penaltySourceField = penaltySource.primitiveFieldRef();
penaltySourceField = -neg(alpha1 - dbba)/runTime.deltaTValue();


// *** Solve momentum equation for Ug
if (momentumType == "slipModel")
{
//    #include "slipVelocity.H"
}
else if (momentumType == "off")
{
    Ug = U;
}

// *** Print the velocity values
Info << "Velocities: "
<< "Min(U) = " << min(mag(U)).value()
<< " Max(U) = " << max(mag(U)).value()
<< " Min(Ug) = " << min(mag(Ug)).value()
<< " Max(Ug) = " << max(mag(Ug)).value()
<< endl;

phiUg = linearInterpolate(Ug) & mesh.Sf();

fvScalarMatrix alphagEqn
(
    fvm::ddt(alphag)
  + fvm::div(phi, alphag)
  - diffusionFlag*Ct*fvm::laplacian(mixture.nuf(), alphag)
        ==
    alphagSource
  + bba*fvm::Sp(penaltySource, alphag)
);


// Only introduce air below h90
alphagEqn.solve();

alphag.correctBoundaryConditions();

Info<< "Phase-g volume fraction = "
    << alphag.weightedAverage(mesh.Vsc()).value()
    << "  Min(alphag) = " << min(alphag).value()
    << "  Max(alphag) = " << max(alphag).value()
    << endl;
