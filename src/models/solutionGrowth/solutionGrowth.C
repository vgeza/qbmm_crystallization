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
    minAbscissa_(dict.lookupOrDefault("minAbscissa", scalar(0))),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", GREAT)),
    kG_("kG", dimLength/dimTime, dict),
    g_(dict.lookupOrDefault("g", scalar(1))),
    Eg_("Eg", dimEnergy/dimMoles, dict),
    alphaG_("alphaG", inv(dimLength), dict),
    Lg_(dict.lookupOrDefault("Lg", scalar(1))),
    kD_(dict.lookupOrDefault<scalar>("kD", 0)),
    d_(dict.lookupOrDefault<scalar>("d", 1)),
    Ed_(dict.lookupOrDefault<scalar>("Ed", 0)),
    alphaD_(dict.lookupOrDefault<scalar>("alphaD", 0)),
    Ld_(dict.lookupOrDefault<scalar>("Ld", 0)),
    minAbscissaD_(dict.lookupOrDefault<scalar>("minAbscissaD", 0)),
    R_("R", dimEnergy/dimMoles/dimTemperature, dict),
    relativeSupersaturation_(false),
    printCoeffs_(dict.getOrDefault<Switch>("printCoeffs", false)),
    celli_(0)
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

    if (mag(kD_) < SMALL || L < minAbscissaD_)
    {
        return 0.0;
    }

    // kD is negative, so the returned rate shrinks the crystals
    return
        kD_
       *Foam::pow((Csat - C)*scale, d_)
       *Foam::exp(-Ed_/(R_.value()*T))
       *Foam::pow(1.0 + alphaD_*L, Ld_);
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

    // TODO: do not print the name of a dimensionedScalar twice
    Info<< "Coefficients:" << endl;
    Info<< "minAbscissa " << minAbscissa_ << endl;
    Info<< "maxAbscissa " << maxAbscissa_ << endl;
    Info<< "kG " << kG_ << endl;
    Info<< "g " << g_ << endl;
    Info<< "Eg " << Eg_ << endl;
    Info<< "alphaG " << alphaG_ << endl;
    Info<< "Lg " << Lg_ << endl;
    Info<< "kD " << kD_ << endl;
    Info<< "d " << d_ << endl;
    Info<< "Ed " << Ed_ << endl;
    Info<< "alphaD " << alphaD_ << endl;
    Info<< "Ld " << Ld_ << endl;
    Info<< "minAbscissaD " << minAbscissaD_ << endl;
    Info<< "R " << R_ << endl;
    Info<< "relativeSupersaturation " << relativeSupersaturation_ << endl;
    Info<< "solubilityCurve " << *solubilityCurve_ << endl;
}

// ************************************************************************* //
