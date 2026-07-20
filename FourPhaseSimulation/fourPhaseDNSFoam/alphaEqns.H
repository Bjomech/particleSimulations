// Alpha transport with MULES subcycling (2512 API)
{
    label nAlphaCorr
    (
        mesh.solution().subDict("PIMPLE").lookupOrDefault<label>("nAlphaCorr", 1)
    );

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(mixture.cAlpha()*phic, max(phic));

    surfaceScalarField phir(phic*mixture.nHatf());

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha1
        (
            fvc::flux
            (
                phi,
                alpha1,
                "div(phi,alpha)"
            )
          + fvc::flux
            (
                phir,
                scalar(1) - alpha1,
                "div(phirb,alpha)"
            )
        );

        MULES::explicitSolve(alpha1, phi, phiAlpha1, 1, 0);

        surfaceScalarField phiAlpha2
        (
            fvc::flux
            (
                phi,
                alpha2,
                "div(phi,alpha)"
            )
          + fvc::flux
            (
               -phir,
                alpha2,
                "div(phirb,alpha)"
            )
        );

        MULES::explicitSolve(alpha2, phi, phiAlpha2, 1, 0);

        alpha3 = scalar(1) - alpha1 - alpha2;
        alpha3.max(0);
        alpha3.min(1);
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
