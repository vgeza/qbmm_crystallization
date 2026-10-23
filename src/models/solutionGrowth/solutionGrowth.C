/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi and 2017 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
-------------------------------------------------------------------------------
2017-03-28 Alberto Passalacqua: Adapted to single scalar calculation.
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solutionGrowth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(solutionGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        solutionGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::solutionGrowth
::solutionGrowth
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    growthModel(dict, mesh),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", 0.0)),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", GREAT)),
    kG_(dimensionedScalar::getOrDefault("kG", dict, dimLength/dimTime, 0)),
    g_(dict.lookupOrDefault("g", 0.0)),
    Eg_(dimensionedScalar::getOrDefault("Eg", dict, dimEnergy/dimMoles, 0)),
    alphaG_(dimensionedScalar::getOrDefault("alphaG", dict, inv(dimLength), 0)),
    Lg_(dict.lookupOrDefault("Lg", 0.0)),
    kD_(dimensionedScalar::getOrDefault("kD", dict, dimLength/dimTime, 0)),
    d_(dict.lookupOrDefault<scalar>("d", 0.0)),
    Ed_(dimensionedScalar::getOrDefault("Ed", dict, dimEnergy/dimMoles, 0)),
    alphaD_(dimensionedScalar::getOrDefault("alphaD", dict, inv(dimLength), 0)),
    Ld_(dict.lookupOrDefault<scalar>("Ld", 0.0)),
    minAbscissaD_(dict.lookupOrDefault<scalar>("minAbscissaD", 0.0)),
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
    relativeSupersaturation_(false),
    CsatName_(dict.lookupOrDefault<word>("Csat", "Csat")),
    CsatPtr_(nullptr),
    printCoeffs_(dict.getOrDefault<Switch>("printCoeffs", false)),
    celli_(0)
{
    solubilityCurve_.reset
    (
        Function1<scalar>::New
        (
            "solubility",
            mesh.lookupObject<IOdictionary>("crystallizationProperties"),
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

    printCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::solutionGrowth
::~solutionGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth
::saturationConcentration
(
    const scalar T
) const
{
    const scalar Tc = T - scalar(273.15);

    return solubilityCurve_->value(Tc);
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth
::growthRate
(
    const scalar L,
    const scalar C,
    const scalar T
) const
{
    const scalar Csat = saturationConcentration(T);

    // Driving force: C - Csat, or (C - Csat)/Csat
    const scalar scale = relativeSupersaturation_ ? 1.0/Csat : 1.0;

    if (C >= Csat)
    {
        return
            kG_.value()
           *Foam::pow((C - Csat)*scale, g_)
           *Foam::exp(-Eg_.value()/(R_.value()*T))
           *Foam::pow(1.0 + alphaG_.value()*L, Lg_);
    }

    if (mag(kD_.value()) < SMALL || L < minAbscissaD_)
    {
        return 0.0;
    }

    // kD is negative, so the returned rate shrinks the crystals
    return
        kD_.value()
       *Foam::pow((Csat - C)*scale, d_)
       *Foam::exp(-Ed_.value()/(R_.value()*T))
       *Foam::pow(1.0 + alphaD_.value()*L, Ld_);
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth::Kg
(
    const scalar& abscissa,
    const bool lengthBased,
    const label environment
) const
{
    const volScalarField& C =
        mesh_.lookupObject<volScalarField>("C");

    const volScalarField& T =
        mesh_.lookupObject<volScalarField>("T");

    (*CsatPtr_)[celli_] = saturationConcentration(T[celli_]);

    return
        growthRate
        (
            abscissa,
            C[celli_],
            T[celli_]
        )
       *pos0(abscissa - minAbscissa_)
       *neg0(abscissa - maxAbscissa_);
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth
::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    celli_ = celli;

    return growthModel::phaseSpaceConvection(momentOrder, celli, quadrature);
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::solutionGrowth
::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    celli_ = celli;

    return growthModel::phaseSpaceConvection(momentOrder, celli, quadrature);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void
Foam::populationBalanceSubModels::growthModels::solutionGrowth
::printCoeffs()
{
    if (!printCoeffs_)
    {
        return;
    }

    Info<< "Coefficients:" << endl;
    Info<< "minAbscissa " << minAbscissa_ << endl;
    Info<< "maxAbscissa " << maxAbscissa_ << endl;
    Info<< kG_ << endl;
    Info<< "g " << g_ << endl;
    Info<< Eg_ << endl;
    Info<< alphaG_ << endl;
    Info<< "Lg " << Lg_ << endl;
    Info<< kD_ << endl;
    Info<< "d " << d_ << endl;
    Info<< Ed_ << endl;
    Info<< alphaD_ << endl;
    Info<< "Ld " << Ld_ << endl;
    Info<< "minAbscissaD " << minAbscissaD_ << endl;
    Info<< R_ << endl;
    Info<< "relativeSupersaturation " << relativeSupersaturation_ << endl;
    Info<< "solubilityCurve " << *solubilityCurve_ << endl;
}

// ************************************************************************* //
