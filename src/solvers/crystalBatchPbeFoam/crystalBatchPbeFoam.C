/*---------------------------------------------------------------------------*\
Application
    crystalBatchPbeFoam

Description
    Zero-dimensional (well mixed batch) crystallization solver.

    Solves the univariate population balance together with the solute mass
    balance for a perfectly mixed batch crystallizer, with the temperature
    prescribed as a Function1 of time.  No momentum, energy or species
    transport is solved, so this is the direct OpenFOAM counterpart of the
    0-D model of

        Y. C. Liu et al., Cryst. Growth Des. 2020, 20, 5235-5250

    and is meant for verifying the growth/nucleation/dissolution sub-models
    against that reference before they are used in 3-D.

        dm_k/dt = (population balance sources)
        dC/dt   = -kv rhoCrystal dm3/dt

    Several populations (solid forms, e.g. a cocrystal and its pure
    components) can share the solution: list them in `populations` in
    populationBalanceProperties, list the solutes in `solutes` in
    crystallizationProperties, and give each population its rhoCrystal, kv
    and `consumption` of every solute there (see createFields.H).  Then

        dC_s/dt = - sum_p consumption_ps kv_p rhoCrystal_p dm3_p/dt

    The mesh may contain any number of cells; every cell is integrated
    independently, so a single cell mesh is the usual choice.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Function1.H"
#include "OFstream.H"
#include "populationBalanceModel.H"
#include "crystalSupersaturation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Zero-dimensional batch crystallization solver: population balances "
        "coupled to the solute mass balances at a prescribed temperature."
    );

    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;

    // Header of the csv log.  A single population keeps the columns it
    // always had; with several, every solute and every population gets its
    // own columns.
    autoPtr<OFstream> logFile;
    {
        fileName dir(runTime.path()/"postProcessing"/"batch");
        mkDir(dir);
        logFile.reset(new OFstream(dir/"batch.csv"));

        if (!multiPopulation)
        {
            logFile()
                << "t,T,Csat,C,supersaturation,"
                << "m0,m1,m2,m3,m4,d10,d43,M" << endl;
        }
        else
        {
            logFile() << "t,T";

            for (const word& s : soluteNames)
            {
                logFile() << ',' << s;
            }

            forAll(populationNames, p)
            {
                const word& pn = populationNames[p];

                if (supersaturation.set(p))
                {
                    logFile() << ',' << pn << ".Ceq," << pn << ".S";
                }

                forAll(moments[p], k)
                {
                    logFile() << ',' << pn << ".m" << k;
                }

                logFile()
                    << ',' << pn << ".d10," << pn << ".d43," << pn << ".M";
            }

            logFile() << endl;
        }
    }

    // Third moment of every population at the start of the step
    List<scalarField> m3Old(nPopulations);

    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Prescribed temperature -------------------------------------
        const scalar Tnow = temperatureProgram->value(runTime.value());
        T == dimensionedScalar("T", dimTemperature, Tnow);

        // --- Population balances, with C and T frozen over the step ------
        //  Every population sees the same solution state; the solutes are
        //  only updated once all of them have been integrated, so the
        //  result does not depend on the order of the populations.
        forAll(populations, p)
        {
            m3Old[p] = moments[p][3].primitiveField();

            populations[p].solve();
        }

        // --- Solute mass balances ---------------------------------------
        //  dC_s = - sum_p consumption[p][s] kv_p rhoCrystal_p dm3_p
        forAll(solutes, s)
        {
            volScalarField& C = solutes[s];

            forAll(populations, p)
            {
                const scalar coeff = consumption[p][s];

                if (coeff == 0)
                {
                    continue;
                }

                C.primitiveFieldRef() -=
                    coeff*rhoCrystal[p]*kv[p]
                   *(moments[p][3].primitiveField() - m3Old[p]);
            }

            if (limitC)
            {
                C.max(0.0);
            }

            C.correctBoundaryConditions();
        }

        // --- Report ------------------------------------------------------
        if (!multiPopulation)
        {
            const volScalarField& C = solutes[0];
            const UPtrList<const volScalarField>& m = moments[0];

            const scalar Csat = supersaturation[0].solubility(Tnow);
            const scalar d43 = m[3][0] > SMALL ? m[4][0]/m[3][0] : 0.0;
            const scalar d10 = m[0][0] > SMALL ? m[1][0]/m[0][0] : 0.0;
            const scalar M = kv[0]*rhoCrystal[0]*m[3][0];

            Info<< "    T = " << Tnow << " K, Csat = " << Csat
                << ", C = " << C[0] << " kg/m3, sigma = " << C[0] - Csat
                << nl
                << "    m0 = " << m[0][0] << ", m3 = " << m[3][0]
                << ", M = " << M << " kg/m3, d43 = " << d43 << " m"
                << endl;

            logFile()
                << runTime.value() << ',' << Tnow << ',' << Csat << ','
                << C[0] << ',' << C[0] - Csat << ','
                << m[0][0] << ',' << m[1][0] << ',' << m[2][0] << ','
                << m[3][0] << ',' << m[4][0] << ','
                << d10 << ',' << d43 << ',' << M << endl;
        }
        else
        {
            Info<< "    T = " << Tnow << " K";
            logFile() << runTime.value() << ',' << Tnow;

            forAll(solutes, s)
            {
                Info<< ", " << soluteNames[s] << " = " << solutes[s][0];
                logFile() << ',' << solutes[s][0];
            }

            Info<< endl;

            forAll(populations, p)
            {
                const UPtrList<const volScalarField>& m = moments[p];

                const scalar d43 =
                    m.size() > 4 && m[3][0] > SMALL ? m[4][0]/m[3][0] : 0.0;
                const scalar d10 = m[0][0] > SMALL ? m[1][0]/m[0][0] : 0.0;
                const scalar M = kv[p]*rhoCrystal[p]*m[3][0];

                Info<< "    " << populationNames[p] << ':';

                if (supersaturation.set(p))
                {
                    const scalar Ceq = supersaturation[p].Ceq(Tnow);
                    const scalar S = supersaturation[p].S(0);

                    Info<< " S = " << S << ',';
                    logFile() << ',' << Ceq << ',' << S;
                }

                forAll(m, k)
                {
                    logFile() << ',' << m[k][0];
                }

                logFile() << ',' << d10 << ',' << d43 << ',' << M;

                Info<< " m0 = " << m[0][0] << ", M = " << M
                    << ", d43 = " << d43 << endl;
            }

            logFile() << endl;
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
