{
    dimensionedScalar DC("DC", dimArea/dimTime, crystallizationProperties);

    const scalar ScT =
        crystallizationProperties.lookupOrDefault<scalar>("ScT", 0.7);

    const scalar kv =
        crystallizationProperties.lookupOrDefault<scalar>("kv", 1.0);

    const bool limitC =
        crystallizationProperties.lookupOrDefault<bool>("limitC", true);

    const bool limitToSaturation =
        crystallizationProperties.lookupOrDefault<bool>
        (
            "limitToSaturation",
            true
        );

    const bool crystalAdvectionCorrection =
        crystallizationProperties.lookupOrDefault<bool>
        (
            "crystalAdvectionCorrection",
            true
        );

    volScalarField DCEff
    (
        IOobject
        (
            "DCEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        DC + turbulence->nut()/ScT
    );

    // Volumetric flux.  Refreshed every time step.
    phiByRho = phi/fvc::interpolate(rho);

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phiByRho, C)
      - fvm::laplacian(DCEff, C)
    );

    CEqn.relax();
    CEqn.solve();

    C.correctBoundaryConditions();

    const volScalarField& m3 =
        mesh.lookupObject<volScalarField>("moment.3.populationBalance");

    const volScalarField m3Old("m3Old", m3);

    // Advective part of the change of m3 over this step
    tmp<volScalarField> tdivM3 = fvc::div(phiByRho, m3);

    populationBalance->solve();

    volScalarField dm3("dm3", m3 - m3Old);

    if (crystalAdvectionCorrection)
    {
        dm3 += runTime.deltaT()*tdivM3();
    }

    tdivM3.clear();

    // Solute consumed by the growth and nucleation of crystals.
    // The split of dm3 is not exact at steep gradients,
    //  clip the exchange at Csat.
    {
        const scalarField& T0 = thermo.T().primitiveField();
        scalarField& Ci = C.primitiveFieldRef();
        const scalarField& dm3i = dm3.primitiveField();

        forAll(Ci, celli)
        {
            const scalar dC = rhoCrystal.value()*kv*dm3i[celli];

            if (!limitToSaturation)
            {
                Ci[celli] -= dC;
                continue;
            }

            const scalar Csat = saturation(T0[celli]);

            if (Ci[celli] > Csat)
            {
                // crystallizing: cannot consume more than the excess
                Ci[celli] = max(Ci[celli] - dC, Csat);
            }
            else
            {
                // dissolving: cannot release more than the deficit
                Ci[celli] = min(Ci[celli] - dC, Csat);
            }
        }
    }

    if (limitC)
    {
        C.max(0.0);
    }

    C.correctBoundaryConditions();
}
