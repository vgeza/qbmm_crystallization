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

    The mesh may contain any number of cells; every cell is integrated
    independently, so a single cell mesh is the usual choice.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Function1.H"
#include "OFstream.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Zero-dimensional batch crystallization solver: population balance "
        "coupled to the solute mass balance at a prescribed temperature."
    );

    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;

    // Header of the csv log
    autoPtr<OFstream> logFile;
    {
        fileName dir(runTime.path()/"postProcessing"/"batch");
        mkDir(dir);
        logFile.reset(new OFstream(dir/"batch.csv"));
        logFile()
            << "t,T,Csat,C,supersaturation,"
            << "m0,m1,m2,m3,m4,d10,d43,M" << endl;
    }

    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Prescribed temperature -------------------------------------
        const scalar Tnow = temperatureProgram->value(runTime.value());
        T == dimensionedScalar("T", dimTemperature, Tnow);

        // --- Population balance, with C and T frozen over the step -------
        scalarField m3Old(m3.primitiveField());

        populationBalance->solve();

        // --- Solute mass balance ----------------------------------------
        C.primitiveFieldRef() -=
            rhoCrystal.value()*kv*(m3.primitiveField() - m3Old);

        if (limitC)
        {
            C.max(0.0);
        }

        C.correctBoundaryConditions();

        // --- Report ------------------------------------------------------
        const scalar Csat = saturation(Tnow);
        const scalar d43 = m3[0] > SMALL ? m4[0]/m3[0] : 0.0;
        const scalar d10 = m0[0] > SMALL ? m1[0]/m0[0] : 0.0;
        const scalar M = kv*rhoCrystal.value()*m3[0];

        Info<< "    T = " << Tnow << " K, Csat = " << Csat
            << ", C = " << C[0] << " mg/mL, sigma = " << C[0] - Csat << nl
            << "    m0 = " << m0[0] << ", m3 = " << m3[0]
            << ", M = " << M << " mg/mL, d43 = " << d43 << " um" << endl;

        logFile()
            << runTime.value() << ',' << Tnow << ',' << Csat << ','
            << C[0] << ',' << C[0] - Csat << ','
            << m0[0] << ',' << m1[0] << ',' << m2[0] << ','
            << m3[0] << ',' << m4[0] << ','
            << d10 << ',' << d43 << ',' << M << endl;

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
