# OPTION FILES

## Overview
An options file contains descriptions of particle decays and AmpGen settings such as input/output locations, global
flags and particle data if it differs from that given in the PDG. The name of an options file conventionally ends in
`.opt`.

## Syntax
This list is not exhaustive, and covers only some common syntax.
### Comments
- Lines beginning with `*`, `#` or a null character (`\0`) are interpreted as comments and are ignored.  
Empty lines are also ignored.

### Keywords
- `CouplingConstant::Coordinates (polar|cartesian)`: specify the co-ordinate system used.
- `CouplingConstant::AngularUnits (deg|rad)` : specify the angular unit if using polar co-ordinates.
- `EventType` : specify type of event.
- `Import` : import another options file.
- `Type` : specify the generator configuration to use; must be of type `AmpGen::generatorType`.

#### Particle Properties
The particles available and their default properties can be found in `options/mass_width.csv` using the MC format of
the 2008 PDG. Additional pseudoparticles, such as nonresonant states and terms for use in conjunction with K-matrices
are defined by `options/MintDalitzSpecialParticles.csv`. Any additional user defined particles should be added here. For
the default lineshape (the relavistic Breit-Wigner or BW), there are three parameters: The mass, the width and the Blatt
Weisskopf radius. These default to their PDG values, but can be overridden in the options file with parameters:
`particleName_mass`, `particleName_width`, `particleName_radius`. To vary the mass of e.g. the K1(1270)+ meson, the
line:
```
K(1)(1270)+_mass 0 1.270 0.01
```
should be added to the user configuration. Other lineshapes may define other parameters, for example channel couplings
or pole masses in the case of the K-matrices can be set or varied in a similar manner.

#### Spin Formalism
AmpGen implements both the covariant tensor (or Rarita-Schwinger) and canonical helicity formalism for describing the
angular momentum component of decays. Both formalisms refer to states of well-defined orbital angular momentum, as
opposed to the helicity states, as the states with well-defined orbital angular momentum have a straightforward parity
and momentum dependences. The default formalism is the covariant tensor formalism, but this can be switched to the
canonical formalism changing the flag
- `Particle::SpinFormalism` : choose which spin formalism to use. Allowed values: `Canonical`, `Covariant`.
This option can also be set in an individual decay chain's decay descriptor, e.g.
```
D0[SpinFormalim=Canonical]{K*(892)bar0,rho(770)0}
```
When using the canonical formalism, the user can also specify systems of helicty coupling:
```
D0[SpinFormalism=Canonical;helAmp=Long]{K*(892)bar0,rho(770)0}
D0[SpinFormalism=Canonical;helAmp=t1]{K*(892)bar0,rho(770)0}
D0[SpinFormalism=Canonical;helAmp=t2]{K*(892)bar0,rho(770)0}
```
These helicity systems must then be defined by the user:
```
Long {
  1.0 0 0
}

t1 {
  0.707106781 +1 +1
  0.707106781 -1 -1
}

t2 {
   0.707106781 +1 +1
  -0.707106781 -1 -1
}
```
First the coupling, then the two particle helicities are specified.

#### Fit Parameters and Expressions

Parameters can either be specified by three parameters, in the case of a scalar parameter such as a mass or a width, or
with six parameters in the case of a complex parameter such as a coupling. Upper and lower bounds on parameters can
also be set by specifying a parameter with five parameters or ten parameters for a scalar or complex, respectively. For
example, if we wished to restrict the allowed values of K meson mass:
```
K(1)(1270)+_mass 0 1.27 0.01 0.0 2.0
```
Parameters can also be related to each other via expressions, Suppose for example we have K(1270)+ and in the same fit
K(1270)-, The properties of one can be allowed to vary and the other fixed to the same value, using:
```
K(1)(1270)+_mass 0 1.27 0.01 0.0 2.0
K(1)(1270)bar-_mass = K(1)(1270)+_mass
```

Parameter expressions are whitespace delimited due to the abundance of odd glyphs such as brackets and +/- in the names
of parameters. Expressions support the binary operations`+-/*`, as well as common unary functions such as sqrt,
trigonometric functions etc.


## Particle Data
An example of specifying data for a decay is below.
```
D0{K*(892)bar0{K-,pi+},rho(770)0{pi+,pi-}}      0       0.196     0.001     0       -22.4       0.4
```
The decay products of particles are provided within curly braces. These can either be provided inline or split over
multiple steps.
In general a complex number is needed to specify an amplitude, so two parameters must be specified.
These parameters are provided in the form of three numbers, the *fix* flag; the initial value; and the step size.
```
#decay chain                                    ----Real/Amplitude----      ----Imaginary/Phase----
#                                               Fix?    Value     Step      Fix?    Value       Step
D0{K*(892)bar0{K-,pi+},rho(770)0{pi+,pi-}}      0       0.196     0.001     0       -22.4       0.4
```
Whether these parameters are Real/Imaginary parts or Magnitude/Phase is specified with `CouplingConstant::Coordinates`.

#### Fix Flag
The fix flag takes three possible options:
- `0`; Free, and a nonzero step size.
- `2`; Fixed
- `3`; Compile-Time-Constant, which indicates that this parameter should be treated as a (JIT) compile-time constant.
In some cases, this allows for more aggressive optimisation.

#### Initial Value and Step Size
The initial value and step size can either be specified in term of an Argand or Polar complex number `x+iy` or
`rexp(theta)` respectively. This is chosen using the `CouplingConstant::Coordinates` option.

## Examples
Example options files can be found at `./options/`.
