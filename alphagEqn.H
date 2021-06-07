/*              COMPUTE THE FREE SURFACE DETECTION FUNCTION                   */

label n = 4;

// Recompute the gradient of alpha1
gradAlpha1 = fvc::grad(alpha1);

// Compute the critical value of grad(alpha1)
dimensionedScalar gradAlpha1Cr("gradAlpha1Cr", 1/(n*averageCellSize));

Info<< "Critical grad(alpha1) value used in free surface detection = " 
    << gradAlpha1Cr.value() << endl; 

// Distance from the focus point of the parabola to its top point
dimensionedScalar aParabola
(
    "aParabola",
     0.25*sqr(gradAlpha1Cr - deltaFunctionCutoff)
);
    
    // delta_fc thickness
    scalar thickness =  5e-3;

    deltaFreeSurface = pos(alpha1 - 0.5);
    deltaFreeSurface *= neg(fvc::average(fvc::interpolate(deltaFreeSurface)) -1);


    // linear function for reducing delta_fc
    scalar factor = (1 - averageCellSize.value()/thickness);
    scalar i = 1;

    while (factor > 0)
    {
        volScalarField temp(pos(fvc::average(fvc::interpolate(deltaFreeSurface))));
        temp = (temp - deltaFreeSurface)*factor;
        deltaFreeSurface += temp;
        i++;
        factor =  1 - i*averageCellSize.value()/thickness;
    }

    // Kill the expansion upwards
    deltaFreeSurface *= pos(alpha1 - 0.5);


/*           COMPUTE THE SOURCE TERM IN THE ALPHAG EQUATION                   */

// Interface-normal unit vectors
volVectorField nHat(gradAlpha1/(mag(gradAlpha1) + mixture.deltaN()));

// Interface-normal velocity derivative, keep positive values only
UnHat = U & nHat;
dUdn = Foam::posPart(fvc::grad(UnHat) & nHat);

forAll(dUdn,i)
{
    if(mesh.C()[i][0] < 0.5){
        // Info << "cellCentres[i]" << mesh.C()[i] << endl;
        deltaFreeSurface[i] = 0;
    }
}


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
