# GENERATOR

## Name
Generator - generate events.

## Synopsis
    Generator OPTIONS [--EventType EVENT_TYPE] [--nEvents EVENTS] [--BlockSize BLOCKSIZE] [--Seed SEED]
              [--Output OUTPUT] [--Library LIBRARY] [--nBins BINS] [--nCores CORES]
              [--GenerateTimeDependent] [--help] [--Type TYPE]

## Description
Generator generates events according to a model and outputs the result as a ROOT binary. This should contain a tree
of candidates with full four-vectors, as well as one- and two-dimensional projections.

## Options
### Generic Program Information
**--help**  
Output a usage message.

### Required Arguments
**OPTIONS**  
Options file containing settings and the description of particle decays.

**--EventType=EVENT_TYPE**  
EventType to generate, in format *parent daughter1 daughter2 ...*

### Optional Arguments
**--nEvents=EVENTS**  
Total number of events to generate. Defaults to 1.

**--BlockSize=BLOCKSIZE**  
Number of events to generate per block. Defaults to 10000.

**--Seed=SEED**  
Random integer seed to use in event generation. Used by ROOT [TRandom3](https://root.cern.ch/doc/master/classTRandom3.html)
to generate a random number.

**--Output=OUTPUT**  
Relative filepath to write output data to. Defaults to `Generate_Output.root`

**--Library=LIBRARY**  
Name of library to use for a fixed library generation.

**--nBins=BINS**  
Number of bins for monitoring plots. Defaults to 100.

**--nCores=CORES**  
Number of cores to use; OpenMP only. Defaults to 12.

**--GenerateTimeDependent**  
Boolean flag to include possible time dependence of the amplitude. Off by default.

**--Type=TYPE**
Generator configuration to use, of type AmpGen::generatorType:
* CoherentSum  : Full phase-space generator with (pseudo)scalar amplitude
* PolarisedSum : Full phase-space generator with particles carrying spin in the initial/final states
* FixedLib     : Full phase-space generator with an amplitude from a precompiled library
* RGenerator   : Recursive phase-space generator for intermediate (quasi)stable states such as the D-mesons

## Generator Configurations
The Generator can be configured to generate amplitudes in a number of different ways, specified via the `--Type` option.

### CoherentSum
A coherent sum of resonant contributions, where the probability density is found from summing the products of isobar
channel couplings and amplitudes. Used for (pseudo-)scalar amplitudes; see [here](https://goofit.github.io/AmpGen/d1/d91/class_amp_gen_1_1_coherent_sum.html).

#### PolarisedSum
A sum of contributions taking into account the spins in the initial/final states; see [here](https://goofit.github.io/AmpGen/d4/d31/class_amp_gen_1_1_polarised_sum.html#abec291964694d5a25040f22e5188b464)
for API documentation.

#### FixedLib
Use amplitudes from a prebuilt library; specify this using `--Library`.

#### RGenerator
Recursive phase-space generation, for quasistable states. Each quasistable particle's decay is handled by the
[PhaseSpace](https://goofit.github.io/AmpGen/d6/d53/class_amp_gen_1_1_phase_space.html) generator. See the API
documentation for RecursivePhaseSpace [here](https://goofit.github.io/AmpGen/dd/daf/class_amp_gen_1_1_recursive_phase_space.html).

## Examples
Example options files are provided with AmpGen in `./options/`.

To generate one million D0 -> K- pi+ pi+ pi- events, run:
```
./Generator options/D02Kpipipi.opt --EventType "D0 K- pi+ pi+ pi-" --nEvents 1000000 --Output=my_file.root
```
The output can then be analysed with e.g.
```
root my_file.root  
```
The options file `D02Kpipipi.opt` describes the Cabbibo favoured decay D0 -> K3pi. A minimal options file might look like:
```
EventType D0 K- pi+ pi+ pi-
D0{K*(892)0{K+,pi-},rho(770)0{pi+,pi-}}     2   1   0       2   0   0
```
See the docs on options files for more.

