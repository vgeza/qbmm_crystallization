/*---------------------------------------------------------------------------*\
  Solution crystallization nucleation model for OpenQBMM.

Class
    Foam::populationBalanceSubModels::nucleationModels::solutionNucleation

Description
    Secondary nucleation from a supersaturated solution, following

        Y. C. Liu et al., "Population Balance Model Development Verification
        and Validation of Cooling Crystallization of Carbamazepine",
        Cryst. Growth Des. 2020, 20, 5235-5250, eq. (13) / (25):

            B = kS (C - Csat(T))^s exp(-Es/(R T)) M^Ms
            M = kv rhoCrystal m3

    The nuclei are born at the size nucleationSize (L0), so the source term
    of the moment of order k is

            S_k = B L0^k

    The paper takes L0 -> 0, which only feeds moment 0; that is the default
    here.  Note that a finite L0 injects solid mass at a rate B kv rhoc L0^3,
    which feeds back into M^Ms and can trigger a runaway.  Values of L0 above
    a few micrometres are therefore not advisable.

    With `supersaturation relative` the driving force is the relative
    supersaturation sigma = (C - Csat)/Csat = S - 1 instead of C - Csat, i.e.
    B = kS sigma^s ..., the form of most power-law primary nucleation fits.

    All coefficients are in 
        L in um
        C and M in mg/mL
        f in #/um/mL
        rates are per second

SourceFiles
    solutionNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef solutionNucleation_H
#define solutionNucleation_H

#include "Function1.H"
#include "nucleationModel.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

class solutionNucleation
:
    public nucleationModel
{
private:

    // Private data

        //- Nucleation rate constant, per second
        scalar kS_;

        //- Activation energy and gas constant
        dimensionedScalar Es_;
        dimensionedScalar R_;

        //- Exponent on the supersaturation
        scalar s_;

        //- Exponent on the suspension density
        scalar Ms_;

        //- Crystal density and volume shape factor used to build the
        //  suspension density M = kv rhoc m3
        scalar rhoCrystal_;
        scalar kv_;

        //- Size at which new crystals enter the population balance
        scalar nucleationSize_;

        //- Optional dependence on the turbulent dissipation rate,
        //  B *= (epsilon/epsilonRef)^epsilonExp.  With epsilonExp = 0 (the
        //  default) the model reduces to the 0-D form of the paper, in which
        //  the effect of the stirring is lumped into kS.
        scalar epsilonExp_;
        scalar epsilonRef_;
        word epsilonName_;

        //- Solubility curve, Csat(Tc) with Tc in degC
        autoPtr<Function1<scalar>> solubilityCurve_;

        //- Driving force (C - Csat)/Csat instead of C - Csat
        bool relativeSupersaturation_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Concentration, temperature and third moment fields
        const volScalarField& C_;
        const volScalarField& T_;
        const volScalarField& m3_;

        //- Turbulent dissipation rate, only looked up when epsilonExp != 0
        const volScalarField* epsilonPtr_;


    // Private Member Functions

        //- Solubility [mg/mL] at temperature T [K]
        scalar saturationConcentration(const scalar T) const;


protected:
    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();


public:

    //- Runtime type information
    TypeName("solution");


    // Constructors

        //- Construct from components
        solutionNucleation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~solutionNucleation();


    // Member Functions

        //- Nucleation rate of the moment of the given order
        virtual scalar nucleationSource
        (
            const label& momentOrder,
            const label celli,
            const label environment = 0
        ) const;
};

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

#endif

// ************************************************************************* //
