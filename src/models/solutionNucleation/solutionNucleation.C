/*---------------------------------------------------------------------------*\
  Solution crystallization nucleation model for OpenQBMM.
\*---------------------------------------------------------------------------*/

#include "solutionNucleation.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(solutionNucleation, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        solutionNucleation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::solutionNucleation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    kS_(dimensionedScalar::getOrDefault("kS", dict, inv(dimTime*dimVolume), 0)),
    Es_(dimensionedScalar::getOrDefault("Es", dict, dimEnergy/dimMoles, 0)),
    R_
    (
        dimensionedScalar::getOrDefault
        (
            "R",
            mesh.lookupObject<IOdictionary>("crystallizationProperties"),
            dimEnergy/dimMoles/dimTemperature,
            constant::physicoChemical::R.value()
        )
    ),
    s_(dict.lookupOrDefault<scalar>("s", 1.0)),
    Ms_(dict.lookupOrDefault<scalar>("Ms", 2.0)),
    // rhoCrystal and kv: from this model's dictionary if given there (one
    // population of several), else from the top level of
    // crystallizationProperties
    rhoCrystal_
    (
        dimensionedScalar::getOrDefault
        (
            "rhoCrystal",
            dict.found("rhoCrystal")
          ? dict
          : mesh.lookupObject<IOdictionary>("crystallizationProperties"),
            dimMass/dimVolume,
            1.2e-9
        )
    ),
    kv_
    (
        (
            dict.found("kv")
          ? dict
          : mesh.lookupObject<IOdictionary>("crystallizationProperties")
        ).lookupOrDefault<scalar>("kv", 1.0)
    ),
    nucleationSize_(dict.lookupOrDefault<scalar>("nucleationSize", 0.0)),
    epsilonExp_(dict.lookupOrDefault<scalar>("epsilonExp", 0.0)),
    epsilonRef_(dict.lookupOrDefault<scalar>("epsilonRef", 1.0)),
    epsilonName_(dict.lookupOrDefault<word>("epsilon", "epsilon")),
    supersaturation_(dict, mesh),
    kb1_(dict.lookupOrDefault<scalar>("primaryRate", 0)),
    cntA_(0),
    m3_
    (
        mesh.lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("m3", "moment.3.populationBalance")
        )
    ),
    epsilonPtr_(nullptr),
    CsatName_(dict.lookupOrDefault<word>("Csat", "Csat")),
    BrateName_(dict.lookupOrDefault<word>("Brate", "Brate")),
    CsatPtr_(nullptr),
    BratePtr_(nullptr),
    printCoeffs_(dict.getOrDefault<Switch>("printCoeffs", false))
{
    if (kb1_ > 0)
    {
        // Classical nucleation theory: 16 pi v^2 sigma^3 / (3 k^3) [K^3]
        const scalar sigma = dict.get<scalar>("interfacialEnergy");   // J/m2
        const scalar v = dict.get<scalar>("molecularVolume");         // m3
        const scalar kB = 1.380649e-23;

        cntA_ =
            16*constant::mathematical::pi*sqr(v)*pow3(sigma)/(3*pow3(kB));

        Info<< "solutionNucleation: primary nucleation " << kb1_
            << " exp(-" << cntA_ << "/(T^3 ln(S)^2))" << endl;
    }

    if (mag(epsilonExp_) > SMALL)
    {
        epsilonPtr_ = &mesh.lookupObject<volScalarField>(epsilonName_);

        Info<< "solutionNucleation: secondary nucleation scaled by ("
            << epsilonName_ << '/' << epsilonRef_ << ")^" << epsilonExp_
            << endl;
    }

    if (!mesh.foundObject<volScalarField>(CsatName_))
    {
        Info<< CsatName_ << " not found, creating\n";

        CsatPtr_ =
        new volScalarField
        (
            IOobject
            (
                "Csat",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimDensity, Zero),
            "zeroGradient"
        );
    }
    else
    {
        Info<< "Using existing " << CsatName_ << "\n";
        CsatPtr_ = &mesh.lookupObjectRef<volScalarField>(CsatName_);
    }

    if (!mesh.foundObject<volScalarField>(BrateName_))
    {
        Info<< BrateName_ << " not found, creating\n";

        BratePtr_ =
        new volScalarField
        (
            IOobject
            (
                BrateName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", inv(dimVolume*dimTime), Zero),
            "zeroGradient"
        );
    }
    else
    {
        Info<< "Using existing " << BrateName_ << "\n";
        BratePtr_ = &mesh.lookupObjectRef<volScalarField>(BrateName_);
    }

    printCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::~solutionNucleation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::saturationConcentration(const scalar T) const
{
    return supersaturation_.Ceq(T);
}


Foam::scalar
Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::nucleationSource
(
    const label& momentOrder,
    const label celli,
    const label environment
) const
{
    // Effective concentration and its equilibrium value; for a single
    // solute just C and Csat (see crystalSupersaturation.H)
    const scalar T = supersaturation_.T(celli);
    const scalar C = supersaturation_.Ceff(celli);
    const scalar Csat = saturationConcentration(T);

    (*CsatPtr_)[celli] = Csat;

    // Driving force: C - Csat, or (C - Csat)/Csat
    scalar deltaC = max(C - Csat, scalar(0));

    if (supersaturation_.relative())
    {
        deltaC /= Csat;
    }

    if (deltaC < SMALL)
    {
        (*BratePtr_)[celli] = 0.0;
        return 0.0;
    }

    // Suspension density, kg solid / m3
    const scalar M = max(kv_*rhoCrystal_.value()*m3_[celli], scalar(0));

    scalar B =
        kS_.value()
       *Foam::pow(deltaC, s_)
       *Foam::exp(-Es_.value()/(R_.value()*T))
       *Foam::pow(M, Ms_);

    if (epsilonPtr_)
    {
        B *=
            Foam::pow
            (
                max((*epsilonPtr_)[celli], scalar(0))/epsilonRef_,
                epsilonExp_
            );
    }

    // Optional primary nucleation from classical nucleation theory,
    // B1 = kb1 exp(-A/(T^3 ln(S)^2)), added to the (secondary) term above
    // and not scaled by the dissipation.  The chemical potential difference
    // per formula unit of a solid of n = sum(nu_i) solute molecules is
    // kT ln(prod C_i^nu_i/Ksp) = n kT ln(S), hence the factor n (1 for a
    // single solute).
    if (kb1_ > 0)
    {
        const scalar lnS = supersaturation_.n()*Foam::log(C/Csat);

        if (lnS > SMALL)
        {
            B += kb1_*Foam::exp(-cntA_/(pow3(T)*sqr(lnS)));
        }
    }

    if (momentOrder == 0)
    {
        (*BratePtr_)[celli] = B;
        return B;
    }

    return B*Foam::pow(nucleationSize_, momentOrder);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void
Foam::populationBalanceSubModels::nucleationModels::solutionNucleation
::printCoeffs()
{
    if (!printCoeffs_)
    {
        return;
    }

    Info<< "Coefficients:" << endl;
    Info<< kS_ << endl;
    Info<< Es_ << endl;
    Info<< R_ << endl;
    Info<< "s " << s_ << endl;
    Info<< "Ms " << Ms_ << endl;
    Info<< rhoCrystal_ << endl;
    Info<< "kv " << kv_ << endl;
    Info<< "nucleationSize " << nucleationSize_ << endl;
    Info<< "epsilonExp " << epsilonExp_ << endl;
    Info<< "epsilonRef " << epsilonRef_ << endl;
    Info<< "epsilon " << epsilonName_ << endl;
    Info<< "relativeSupersaturation " << supersaturation_.relative() << endl;
    Info<< "solutes " << supersaturation_.soluteNames() << endl;
    Info<< "solubilityCurve " << supersaturation_.solubilityCurve() << endl;
    Info<< "primaryRate " << kb1_ << endl;
}

// ************************************************************************* //
