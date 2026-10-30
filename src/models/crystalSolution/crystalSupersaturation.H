/*---------------------------------------------------------------------------*\
Class
    Foam::crystalSupersaturation

Description
    Supersaturation of a solution with respect to one solid form, shared by
    the solution growth and nucleation models and by the solvers.

    The solid may be built from several solutes (a cocrystal A_a B_b, a salt,
    a solvate), each with a stoichiometric exponent nu_i:

        solutes   { CA 1; CB 1; }      // field name and exponent nu_i

    and its solubility is then a solubility product Ksp(T).  With
    n = sum(nu_i) the effective and the equilibrium concentrations are

        Ceff = (prod C_i^nu_i)^(1/n),    Ceq = Ksp(T)^(1/n),    S = Ceff/Ceq

    so that S is the usual mean-ionic / cocrystal supersaturation ratio and
    the driving force is

        absolute:   Ceff - Ceq
        relative:   (Ceff - Ceq)/Ceq = S - 1

    With a single solute of exponent 1 (the default, solutes { C 1; }) this
    is exactly C - Csat and (C - Csat)/Csat, evaluated in the same order of
    operations as before, so single-solute results are bit-for-bit unchanged.

    Entries read from the dictionary:

        solutes          { C 1; }           optional, default from `C`
        C                C;                 legacy: name of the single solute
        T                T;                 temperature field name
        solubility       <Function1>;       Csat or Ksp of T [degC], see below
        supersaturation  absolute;          or relative

    `solubility` is read from the given dictionary if it is there, and
    otherwise from the top level of constant/crystallizationProperties, where
    a single-population case keeps it.  Several populations can thus each
    carry their own Csat or Ksp curve in their own dictionary.

    Fields are looked up in the mesh registry on first use, so the object
    may be constructed before the fields exist.

\*---------------------------------------------------------------------------*/

#ifndef crystalSupersaturation_H
#define crystalSupersaturation_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "IOdictionary.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class crystalSupersaturation
{
    // Private data

        const fvMesh& mesh_;

        //- Solute field names and their stoichiometric exponents
        wordList soluteNames_;
        scalarList nu_;

        //- Sum of the exponents
        scalar n_;

        //- Single solute with exponent 1: the legacy formulae, bit-for-bit
        bool single_;

        word TName_;

        //- Csat (one solute) or Ksp as a function of T [degC]
        autoPtr<Function1<scalar>> solubilityCurve_;

        bool relative_;

        //- Resolved on first use
        mutable UPtrList<const volScalarField> C_;
        mutable const volScalarField* T_;


    // Private Member Functions

        void resolve() const
        {
            if (T_)
            {
                return;
            }

            T_ = &mesh_.lookupObject<volScalarField>(TName_);

            C_.resize(soluteNames_.size());
            forAll(soluteNames_, i)
            {
                C_.set(i, &mesh_.lookupObject<volScalarField>(soluteNames_[i]));
            }
        }


public:

    // Constructors

        crystalSupersaturation(const dictionary& dict, const fvMesh& mesh)
        :
            mesh_(mesh),
            soluteNames_(),
            nu_(),
            n_(0),
            single_(false),
            TName_(dict.getOrDefault<word>("T", "T")),
            solubilityCurve_(nullptr),
            relative_(false),
            C_(),
            T_(nullptr)
        {
            solubilityCurve_.reset
            (
                Function1<scalar>::New
                (
                    "solubility",
                    dict.found("solubility")
                  ? dict
                  : mesh.lookupObject<IOdictionary>
                    (
                        "crystallizationProperties"
                    ),
                    &mesh
                ).ptr()
            );

            if (dict.found("solutes"))
            {
                const dictionary& sDict = dict.subDict("solutes");

                soluteNames_ = sDict.toc();
                nu_.resize(soluteNames_.size());

                forAll(soluteNames_, i)
                {
                    nu_[i] = sDict.get<scalar>(soluteNames_[i]);

                    if (nu_[i] <= 0)
                    {
                        FatalIOErrorInFunction(dict)
                            << "solutes: the exponent of " << soluteNames_[i]
                            << " must be positive" << exit(FatalIOError);
                    }

                    n_ += nu_[i];
                }

                if (soluteNames_.empty())
                {
                    FatalIOErrorInFunction(dict)
                        << "solutes is empty" << exit(FatalIOError);
                }
            }
            else
            {
                soluteNames_ = wordList(1, dict.getOrDefault<word>("C", "C"));
                nu_ = scalarList(1, scalar(1));
                n_ = 1;
            }

            single_ = (soluteNames_.size() == 1 && nu_[0] == 1);

            const word ss
            (
                dict.getOrDefault<word>("supersaturation", "absolute")
            );

            if (ss == "relative")
            {
                relative_ = true;
            }
            else if (ss != "absolute")
            {
                FatalIOErrorInFunction(dict)
                    << "supersaturation must be absolute or relative, not "
                    << ss << exit(FatalIOError);
            }
        }


    // Member Functions

        const wordList& soluteNames() const
        {
            return soluteNames_;
        }

        bool relative() const
        {
            return relative_;
        }

        //- Number of solute molecules per formula unit, sum(nu_i)
        scalar n() const
        {
            return n_;
        }

        //- Csat (one solute) or Ksp at temperature T [K]
        scalar solubility(const scalar T) const
        {
            return solubilityCurve_->value(T - scalar(273.15));
        }

        //- The Csat or Ksp curve, of T [degC]
        const Function1<scalar>& solubilityCurve() const
        {
            return *solubilityCurve_;
        }

        //- Equilibrium concentration Ceq = Ksp^(1/n) at temperature T [K]
        scalar Ceq(const scalar T) const
        {
            const scalar K = solubility(T);

            if (single_)
            {
                return K;
            }

            return Foam::pow(max(K, scalar(0)), 1.0/n_);
        }

        //- Effective concentration (prod C_i^nu_i)^(1/n) in cell celli
        scalar Ceff(const label celli) const
        {
            resolve();

            if (single_)
            {
                return C_[0][celli];
            }

            scalar lnC = 0;

            forAll(C_, i)
            {
                const scalar Ci = C_[i][celli];

                if (Ci <= 0)
                {
                    return 0;
                }

                lnC += nu_[i]*Foam::log(Ci);
            }

            return Foam::exp(lnC/n_);
        }

        scalar T(const label celli) const
        {
            resolve();
            return (*T_)[celli];
        }

        //- Supersaturation ratio S = Ceff/Ceq
        scalar S(const label celli) const
        {
            return Ceff(celli)/Ceq(T(celli));
        }

        //- Signed driving force, Ceff - Ceq or (Ceff - Ceq)/Ceq.
        //  For one solute: C - Csat or (C - Csat)/Csat.
        scalar drivingForce(const label celli) const
        {
            const scalar Ce = Ceq(T(celli));
            const scalar dC = Ceff(celli) - Ce;

            return relative_ ? dC/Ce : dC;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
