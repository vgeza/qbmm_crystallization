/*---------------------------------------------------------------------------*\
  Solution crystallization nucleation model for OpenQBMM.
\*---------------------------------------------------------------------------*/

#include "solutionNucleation.H"
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
    kS_(dict.get<scalar>("kS")),
    Es_("Es", dimEnergy/dimMoles, dict),
    R_("R", dimEnergy/dimMoles/dimTemperature, dict),
    s_(dict.lookupOrDefault<scalar>("s", 1.0)),
    Ms_(dict.lookupOrDefault<scalar>("Ms", 2.0)),
    rhoCrystal_(dict.lookupOrDefault<scalar>("rhoCrystal", 1.2e-9)),
    kv_(dict.lookupOrDefault<scalar>("kv", 1.0)),
    nucleationSize_(dict.lookupOrDefault<scalar>("nucleationSize", 0)),
    epsilonExp_(dict.lookupOrDefault<scalar>("epsilonExp", 0)),
    epsilonRef_(dict.lookupOrDefault<scalar>("epsilonRef", 1.0)),
    epsilonName_(dict.lookupOrDefault<word>("epsilon", "epsilon")),
    relativeSupersaturation_(false),
    printCoeffs_(dict.getOrDefault<Switch>("printCoeffs", false)),
    C_(mesh.lookupObject<volScalarField>(dict.lookupOrDefault<word>("C", "C"))),
    T_(mesh.lookupObject<volScalarField>(dict.lookupOrDefault<word>("T", "T"))),
    m3_
    (
        mesh.lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("m3", "moment.3.populationBalance")
        )
    ),
    epsilonPtr_(nullptr)
{
    solubilityCurve_.reset
    (
        Function1<scalar>::New
        (
            "solubility",
            dict,
            &mesh
        ).ptr()
    );

    const word ss(dict.getOrDefault<word>("supersaturation", "absolute"));

    if (ss == "relative")
    {
        relativeSupersaturation_ = true;
    }
    else if (ss != "absolute")
    {
        FatalIOErrorInFunction(dict)
            << "supersaturation must be absolute or relative, not " << ss
            << exit(FatalIOError);
    }

    if (mag(epsilonExp_) > SMALL)
    {
        epsilonPtr_ = &mesh.lookupObject<volScalarField>(epsilonName_);

        Info<< "solutionNucleation: secondary nucleation scaled by ("
            << epsilonName_ << '/' << epsilonRef_ << ")^" << epsilonExp_
            << endl;
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
    const scalar Tc = T - scalar(273.15);

    return solubilityCurve_->value(Tc);
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
    const scalar Csat = saturationConcentration(T_[celli]);

    // Driving force: C - Csat, or (C - Csat)/Csat
    scalar deltaC = max(C_[celli] - Csat, scalar(0));

    if (relativeSupersaturation_)
    {
        deltaC /= Csat;
    }

    if (deltaC < SMALL)
    {
        return 0.0;
    }

    // Suspension density, mg solid / mL
    const scalar M = max(kv_*rhoCrystal_*m3_[celli], scalar(0));

    scalar B =
        kS_
       *Foam::pow(deltaC, s_)
       *Foam::exp(-Es_.value()/(R_.value()*T_[celli]))
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

    if (momentOrder == 0)
    {
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

    // TODO: do not print the name of a dimensionedScalar twice
    Info<< "Coefficients:" << endl;
    Info<< "kS " << kS_ << endl;
    Info<< "Es " << Es_ << endl;
    Info<< "R " << R_ << endl;
    Info<< "s " << s_ << endl;
    Info<< "Ms " << Ms_ << endl;
    Info<< "rhoCrystal " << rhoCrystal_ << endl;
    Info<< "kv " << kv_ << endl;
    Info<< "nucleationSize " << nucleationSize_ << endl;
    Info<< "epsilonExp " << epsilonExp_ << endl;
    Info<< "epsilonRef " << epsilonRef_ << endl;
    Info<< "epsilon " << epsilonName_ << endl;
    Info<< "relativeSupersaturation " << relativeSupersaturation_ << endl;
    Info<< "solubilityCurve " << *solubilityCurve_ << endl;
}

// ************************************************************************* //
