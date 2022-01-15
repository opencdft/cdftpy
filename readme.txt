Code structure:

    tests/
        collection of tests for automated
        pytest testing
    examples/cdft1d/
        Contains example input files
    data
        Data directory containing
        solvent smdl files
    cdft1d/
        Main directory containing all the python modules
        for 1D CDFT calculations
        cli.py
            wrappers for executable script
        coulomb.py
            Coulomb potential and energy calculations
        diis.py
            DIIS-related routines
        exceptions.py
            Custom exception classes
        globals.py
            Global configuration variables
        io_utils.py
            Input/Output routines
        loggers.py
            Logging related routines
        potential.py
            Short range potential (non-Coulomb)
        rad_fft.py
            Radial FFT functions
        rdf.py
            Density distribution calculations
        rism.py
            RISM solver
        rsdft.py
            RSDFT solver
        simulation.py
            Simulation class
        solvent.py
            Solvent related routine
        units.py
            Units conversion
        viz.py
            Dashboard generation
        workflow.py
            Task wrappers

Installation instructions:

1. Ensure that the version of python is 3.9 or higher
    python3 --version

2. Change into the directory containing code distribution
3. Create virtual environment
    python3 -m venv venv
4. Activate virtual environment
    source venv/bin/activate
5. Install
    pip install -e .
