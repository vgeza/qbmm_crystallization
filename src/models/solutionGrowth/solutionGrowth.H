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
2020-03-12 Alberto Passalacqua: Added access to fvMesh.
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

Class
    Foam::populationBalanceSubModels::growthModels::solutionGrowth

Description
    Growth from, and dissolution into, a solution, following

        Y. C. Liu et al., "Population Balance Model Development Verification
        and Validation of Cooling Crystallization of Carbamazepine",
        Cryst. Growth Des. 2020, 20, 5235-5250, eq. (29) / (30):

        supersaturated, C > Csat:

            G = kG (C - Csat)^g exp(-Eg/(R T)) (1 + alphaG L)^Lg     > 0

        undersaturated, C < Csat:

            D = kD (Csat - C)^d exp(-Ed/(R T)) (1 + alphaD L)^Ld     < 0

    The two branches are returned through the same Kg() interface, so a
    negative rate simply convects the number density towards smaller sizes.
    Dissolution is off by default (kD = 0).

    With `supersaturation relative` the driving force is the relative
    supersaturation sigma = (C - Csat)/Csat = S - 1 instead of C - Csat (and
    (Csat - C)/Csat for dissolution), as in power-law kinetics of the form
    G = kG sigma^g.

    All coefficients are in 
        L in um
        C in mg/mL

SourceFiles
    solutionGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef solutionGrowth_H
#define solutionGrowth_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class solutionGrowth Declaration
\*---------------------------------------------------------------------------*/

class solutionGrowth
:
    public growthModel
{
private:

    // Private data

        //- Minimum value of the abscissa to apply the model
        scalar minAbscissa_;

        //- Maximum value of the abscissa to apply the model
        scalar maxAbscissa_;

        // Growth kinetics
        dimensionedScalar kG_;
        scalar g_;
        dimensionedScalar Eg_;
        dimensionedScalar alphaG_;
        scalar Lg_;

        // Dissolution kinetics (kD <= 0; zero disables dissolution)
        scalar kD_;
        scalar d_;
        scalar Ed_;
        scalar alphaD_;
        scalar Ld_;

        //- Crystals smaller than this do not dissolve any further.  A moment
        //  method has no sink for m0, so without a floor the quadrature
        //  abscissae are driven through zero and the moment set stops being
        //  realizable once the fines have dissolved.  Only active in the
        //  dissolution branch, so it does not prevent nuclei from growing.
        scalar minAbscissaD_;

        //- Gas constant
        dimensionedScalar R_;

        //- Solubility curve, Csat(Tc), Tc in degC
        autoPtr<Function1<scalar>> solubilityCurve_;

        //- Driving force (C - Csat)/Csat instead of C - Csat
        bool relativeSupersaturation_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Cell currently being integrated.
        //  The Kg() interface of OpenQBMM does not carry a cell index, only an
        //  environment index, so a growth rate that depends on the local
        //  concentration and temperature would otherwise evaluate every cell
        //  with the values of cell 0.  phaseSpaceConvection() is given the
        //  cell, so it is recorded here on the way in.
        mutable label celli_;


protected:
    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();


public:

        //- Runtime type information
        TypeName("solution");


    // Constructors

        //- Construct from components
        solutionGrowth
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~solutionGrowth();


    // Member Functions

        //- Solubility [mg/mL] at temperature T [K]
        scalar saturationConcentration
        (
            const scalar T
        ) const;

        //- Signed rate of change of the crystal size [um/s]
        scalar growthRate
        (
            const scalar L,
            const scalar C,
            const scalar T
        ) const;

        //- Growth model
        virtual scalar Kg
        (
            const scalar& abscissa,
            const bool lengthBased,
            const label environment = 0
        ) const;

        //- Record the cell, then use the standard quadrature summation
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        );

        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
