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

const bool limitC =
    crystallizationProperties.lookupOrDefault<bool>("limitC", true);

// Prescribed temperature programme, T [K] as a function of time [s]
autoPtr<Function1<scalar>> temperatureProgram
(
    Function1<scalar>::New("temperatureProgram", crystallizationProperties)
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

// Several populations (solid forms) are given as
//
//     populations (cocrystal formA);
//     cocrystal { populationBalanceModel univariate; univariateCoeffs {...} }
//     formA     { ... }
//
// with moments moment.<k>.<name> and quadratureProperties.<name>.  Without
// `populations` the whole dictionary is the single population
// "populationBalance", as before.
const bool multiPopulation = populationBalanceProperties.found("populations");

const wordList populationNames
(
    multiPopulation
  ? populationBalanceProperties.get<wordList>("populations")
  : wordList(1, word("populationBalance"))
);

if (populationNames.empty())
{
    FatalIOErrorInFunction(populationBalanceProperties)
        << "populations is empty" << exit(FatalIOError);
}

// Solutes, solved for by a mass balance each.  Default: the single solute C.
const wordList soluteNames
(
    crystallizationProperties.getOrDefault<wordList>
    (
        "solutes",
        wordList(1, word("C"))
    )
);

if (!multiPopulation && soluteNames.size() != 1)
{
    FatalIOErrorInFunction(crystallizationProperties)
        << "several solutes need several populations: give `populations` in "
        << "populationBalanceProperties and a subdictionary with "
        << "`consumption` for each of them here" << exit(FatalIOError);
}


Info<< "Reading solute concentrations " << soluteNames << nl << endl;

PtrList<volScalarField> solutes(soluteNames.size());

forAll(soluteNames, i)
{
    solutes.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                soluteNames[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
}

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


// --- Crystal properties of each population ---------------------------------
//
// Single population: rhoCrystal, kv and the solubility curve at the top
// level of crystallizationProperties, the solute C consumed one to one.
//
// Several populations: one subdictionary per population,
//
//     cocrystal
//     {
//         rhoCrystal   1.2e-9;
//         kv           0.5236;
//         consumption  { CA 0.55; CB 0.45; }   // solute per unit kv rho m3
//
//         // optional, only for the report of S in batch.csv
//         solutes      { CA 1; CB 1; }
//         solubility   constant 4;         // Ksp(T [degC]), any Function1
//     }
//
// `consumption` is the amount of each solute taken up by a unit amount of
// this solid, in the units of C: mass fractions of the solid if C and
// rhoCrystal are masses, the stoichiometric coefficients if they are moles.
//
// Without `consumption` it is computed as the mass fractions of the solid,
// nu_s M_s/sum_j(nu_j M_j), from the population's `solutes` exponents and
//
//     molarMass { CA 151.1; CB 123.1; }    // top level, any common unit

const label nPopulations = populationNames.size();
const label nSolutes = soluteNames.size();

scalarList rhoCrystal(nPopulations);
scalarList kv(nPopulations);
List<scalarList> consumption(nPopulations, scalarList(nSolutes, Zero));
PtrList<crystalSupersaturation> supersaturation(nPopulations);

forAll(populationNames, p)
{
    const dictionary& dict =
        multiPopulation
      ? crystallizationProperties.subDict(populationNames[p])
      : crystallizationProperties;

    rhoCrystal[p] =
        dimensionedScalar("rhoCrystal", dimMass/dimVolume, dict).value();

    kv[p] = dict.lookupOrDefault<scalar>("kv", 1.0);

    if (multiPopulation)
    {
        // Given directly, or from the stoichiometry and the molar masses
        const bool fromMolarMass = !dict.found("consumption");

        const dictionary& cDict =
            fromMolarMass ? dict.subDict("solutes") : dict.subDict("consumption");

        scalar molarMassSolid = 0;

        for (const word& s : cDict.toc())
        {
            const label si = soluteNames.find(s);

            if (si < 0)
            {
                FatalIOErrorInFunction(cDict)
                    << "population " << populationNames[p] << " consumes "
                    << s << ", which is not one of the solutes "
                    << soluteNames << exit(FatalIOError);
            }

            consumption[p][si] = cDict.get<scalar>(s);

            if (fromMolarMass)
            {
                consumption[p][si] *=
                    crystallizationProperties.subDict("molarMass")
                    .get<scalar>(s);
                molarMassSolid += consumption[p][si];
            }
        }

        // Mass fraction of each solute in the solid, nu_s M_s/sum(nu M)
        if (fromMolarMass)
        {
            forAll(soluteNames, si)
            {
                consumption[p][si] /= molarMassSolid;
            }
        }
    }
    else
    {
        consumption[p][0] = 1;
    }

    if (!multiPopulation || dict.found("solubility"))
    {
        supersaturation.set(p, new crystalSupersaturation(dict, mesh));
    }

    Info<< "Population " << populationNames[p] << ": rhoCrystal "
        << rhoCrystal[p] << ", kv " << kv[p] << ", consumes";
    forAll(soluteNames, s)
    {
        Info<< ' ' << soluteNames[s] << ' ' << consumption[p][s];
    }
    Info<< endl;
}


// --- Population balances ---------------------------------------------------

PtrList<populationBalanceModel> populations(nPopulations);

forAll(populationNames, p)
{
    Info<< "Creating population " << populationNames[p] << endl;

    populations.set
    (
        p,
        populationBalanceModel::New
        (
            populationNames[p],
            multiPopulation
          ? populationBalanceProperties.subDict(populationNames[p])
          : populationBalanceProperties,
            phi
        )
    );
}

// Moments of each population, moment.<k>.<name> for k = 0, 1, ...
List<UPtrList<const volScalarField>> moments(nPopulations);

forAll(populationNames, p)
{
    DynamicList<const volScalarField*> mp;

    for (label k = 0; ; ++k)
    {
        const word name
        (
            "moment." + Foam::name(k) + "." + populationNames[p]
        );

        if (!mesh.foundObject<volScalarField>(name))
        {
            break;
        }

        mp.append(&mesh.lookupObject<volScalarField>(name));
    }

    if (mp.size() < 4 || (!multiPopulation && mp.size() < 5))
    {
        FatalErrorInFunction
            << "population " << populationNames[p] << " has " << mp.size()
            << " moments; at least m0 .. m3 (m0 .. m4 for a single "
            << "population) are needed" << exit(FatalError);
    }

    moments[p].resize(mp.size());
    forAll(mp, k)
    {
        moments[p].set(k, mp[k]);
    }
}
