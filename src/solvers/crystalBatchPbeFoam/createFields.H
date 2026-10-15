Info<< "Reading crystallizationProperties\n" << endl;

IOdictionary crystallizationProperties
(
    IOobject
    (
        "crystallizationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoCrystal
(
    "rhoCrystal",
    dimMass/dimVolume,
    crystallizationProperties
);

const scalar kv =
    crystallizationProperties.lookupOrDefault<scalar>("kv", 1.0);

const bool limitC =
    crystallizationProperties.lookupOrDefault<bool>("limitC", true);

// Solubility curve, Csat(Tc) with Tc in degC
autoPtr<Function1<scalar>> solubilityCurve
(
    Function1<scalar>::New("solubility", crystallizationProperties, &mesh)
);

// Prescribed temperature programme, T [K] as a function of time [s]
autoPtr<Function1<scalar>> temperatureProgram
(
    Function1<scalar>::New("temperatureProgram", crystallizationProperties)
);

auto saturation = [&](const scalar T) -> scalar
{
    const scalar Tc = T - 273.15;
    return solubilityCurve->value(Tc);
};


Info<< "Reading concentration C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading temperature T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Quiescent flow: the population balance interface needs a flux, but the
// moment advection model of a 0-D case is "noAdvection", so it is never used.
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, Zero),
    "zeroGradient"
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);

Info<< "Reading populationBalanceProperties\n" << endl;

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "populationBalance", populationBalanceProperties, phi
    )
);

const volScalarField& m0 =
    mesh.lookupObject<volScalarField>("moment.0.populationBalance");
const volScalarField& m1 =
    mesh.lookupObject<volScalarField>("moment.1.populationBalance");
const volScalarField& m2 =
    mesh.lookupObject<volScalarField>("moment.2.populationBalance");
const volScalarField& m3 =
    mesh.lookupObject<volScalarField>("moment.3.populationBalance");
const volScalarField& m4 =
    mesh.lookupObject<volScalarField>("moment.4.populationBalance");
