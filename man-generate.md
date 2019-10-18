# GENERATOR

## Name
    Generator - generate events.

## Synopsis
    Generator OPTIONS [--EventType EVENT_TYPE] [--nEvents EVENTS] [--BlockSize BLOCKSIZE] [--Seed SEED]
              [--Output OUTPUT] [--Type TYPE] [--Library LIBRARY] [--nBins BINS] [--nCores CORES]
              [--GenerateTimeDependent] [--help]

## Description
    Generator generates events according to a model and outputs the result as a ROOT binary. This should contain a tree
    of candidates with full four-vectors, as well as one- and two-dimensional projections.

## Options
    ### Generic Program Information
    --help Output a usage message.

    ### Required Arguments
    OPTIONS
        Options file containing settings and the description of particle decays.

    --EventType=EVENT_TYPE
        EventType to generate, in format *parent daughter1 daughter2 ...*

    ### Optional Arguments
    --nEvents=EVENTS
        Totoal number of events to generate. Defaults to 1.

    --BlockSize=BLOCKSIZE
        Number of events to generate per block. Defaults to 10000.

    --Seed=SEED
        Random integer seed to use in event generation. Used by ROOT TRandom3 to generate a random number.

    --Output=OUTPUT
        Relative filepath to write output data to. Defaults to Generate_Output.root

    --Type=TYPE
        Generator configuration to use, of type AmpGen::generatorType.
            * CoherentSum  : Full phase-space generator with (pseudo)scalar amplitude
            * PolarisedSum : Full phase-space generator with particles carrying spin in the initial/final states
            * FixedLib     : Full phase-space generator with an amplitude from a precompiled library
            * RGenerator   : Recursive phase-space generator for intermediate (quasi)stable states such as the D-mesons

    --Library=LIBRARY
        Name of library to use for a fixed library generation.

    --nBins=BINS
        Number of bins for monitoring plots. Defaults to 100.

    --nCores=CORES
        Number of cores to use; OpenMP only. Defaults to 12.

    --GenerateTimeDependent
        Boolean flag to include possible time dependence of the amplitude. Off by default.

## Generator Configurations
The Generator can be configured to generate amplitudes in a number of different ways, specified via the --Type option.

### CoherentSum
    A coherent sum of resonant contributions, where the probability density is found from summing the products of isobar
    channel couplings and amplitudes. See [https://goofit.github.io/AmpGen/d1/d91/class_amp_gen_1_1_coherent_sum.html](/uri "here").

