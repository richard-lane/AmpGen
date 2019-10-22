# OPTION FILES

## Overview
Options files will generally contain the description of one or more particle decays,
as well as other settings such as input/output locations, global flags such as
whether complex numbers should be interpreted as cartesian or polar, and other parameters
such as the masses and widths of particles is these differ from those given in the PDG.  
An options file conventially ends in .opt.
### Table of contents
* [Syntax](#syntax)
* [Particle Data](#particle-data)
* [Advanced Options](#advanced-options)
* [Examples](#examples)


## Syntax
This list is not exhaustive, and covers only some common syntax.
### Comments
- Lines beginning with `*`, `#` or a null character (`\0`) are interpreted as comments and are ignored.  
Empty lines are also ignored.

### Keywords
- `EventType` : specify type of event.
- `Import` : import another options file; see the example below.
- `Type` : specify the generator configuration to use; must be of type `AmpGen::generatorType`.
- `CouplingConstant::Coordinates (polar|cartesian)`: specify the co-ordinate system used.
- `CouplingConstant::AngularUnits (deg|rad)` : specify the angular unit if using polar co-ordinates.

### Using Multiple Options Files
Configuration can be split over multiple files by the using _Import_ keyword, for example, to import the parameters for
the isoscalar K-matrix, the line
```
Import $AMPGENROOT/options/kMatrix.opt
```

can be added to options file. Multiple user configuration files can also be specified by including multiple files on the
command line.


## Particle Data
An example minimal options file for the Generator application could contain:
```
EventType D0 K+ pi- pi- pi+
D0{K*(892)bar0{K-,pi+},rho(770)0{pi+,pi-}}      0       0.196     0.001     0       -22.4       0.4
```
The EventType specifies the initial and final states requested by the user.  
This gives the ordering of particles used in input data, in output code, and used in internal computations.
It also defines how the amplitude source code must be interfaced with external packages, i.e. MC generators such as
EvtGen.  
The decay products of particles are provided within curly braces. These can either be provided inline or split over
multiple steps- see [Cascade Decays](#cascade-decays).

In general a complex number is needed to specify an amplitude, so two parameters must be specified.
These parameters are provided in the form of three numbers, the *fix* flag; the initial value; and the step size.
```
# Decay chain                                   ----Real/Amplitude----      ----Imaginary/Phase----
#                                               Fix?    Value     Step      Fix?    Value       Step
D0{K*(892)bar0{K-,pi+},rho(770)0{pi+,pi-}}      0       0.196     0.001     0       -22.4       0.4
```
Whether these parameters are Real/Imaginary parts or Magnitude/Phase is specified with `CouplingConstant::Coordinates`.

#### Fix Flag
The fix flag determines whether a parameter is allowed to vary and takes three possible options:
- `0`; Free, and a nonzero step size.
- `2`; Fixed
- `3`; Compile-Time-Constant, which indicates that this parameter should be treated as a (JIT) compile-time constant.
In some cases, this allows for more aggressive optimisation.

#### Initial Value and Step Size
The initial value and step size can either be specified in term of an Argand or Polar complex number `x+iy` or
`rexp(theta)` respectively. This is chosen using the `CouplingConstant::Coordinates` option.


### Cascade Decays
Decays can either be specified fully inline, as above, or split into multiple steps, which is useful for treating the so called _cascade_ decays, an example of which is shown below.
```
D0{K(1)(1270)+,pi-}                         0     1      0.1       0     0      0.1   

K(1)(1270)+{rho(770)0{pi+,pi-},K+}          2     1      0         2     0      0
K(1)(1270)+{K*(892)0{K+,pi-},pi+}           0     1      0.1       0     0      0.1
```
The production/decay couplings of the  <img src="figs/tex/fig0.png" style="margin-bottom:-5px" />  resonance are now
defined in terms of the coupling to the  <img src="figs/tex/fig1.png" style="margin-bottom:-5px" />  channel,
which can be useful in making comparisons between different production modes of a resonance.  
Additional care must be taken in such a case to not introduce redundant degrees of freedom.


## Advanced Options

#### Particle Properties and Lineshape Parameters
The particles available and their default properties can be found in `$AMPGENROOT/options/mass_width.csv` using the MC
format of the 2008 PDG. Additional pseudoparticles, such as nonresonant states and terms for use in conjunction with
K-matrices are defined by `options/MintDalitzSpecialParticles.csv`. Any additional user defined particles can be
added here. For the default lineshape (the relavistic Breit-Wigner or BW), there are three parameters: the mass, the
width and the Blatt Weisskopf radius. These default to their PDG values, but can be overridden in the options file
with parameters: `particleName_mass`, `particleName_width`, `particleName_radius`. To vary the mass of e.g. the
<img src="figs/tex/fig5.png" style="margin-bottom:-5px" />  meson, the line:
```
K(1)(1270)+_mass 0 1.270 0.01
```
should be added to the user configuration. Other lineshapes may define other parameters, for example channel couplings
or pole masses in the case of the K-matrices can be set or varied in a similar manner.  

For more details about the API for describing particle decays, see [AmpGen::Particle](https://goofit.github.io/AmpGen/de/dd7/class_amp_gen_1_1_particle.html).  
These options can be used in the Generator application.

#### Spin Formalism
AmpGen implements both the covariant tensor (or Rarita-Schwinger) and canonical helicity formalism for describing the
angular momentum component of decays. Both formalisms refer to states of well-defined orbital angular momentum, as
opposed to the helicity states, as the states with well-defined orbital angular momentum have a straightforward parity
and momentum dependences. The default formalism is the covariant tensor formalism, but this can be switched to the
canonical formalism changing the flag:
- `Particle::SpinFormalism : choose which spin formalism to use. Allowed values: Canonical, Covariant`.

This option can also be set in an individual decay chain's decay descriptor, e.g.
```
D0[SpinFormalim=Canonical]{K*(892)bar0,rho(770)0}
```

selects the S-wave of the  <img src="figs/tex/fig13.png" style="margin-bottom:0px" />  system. The user can also
specify systems of helicity couplings in the canonical formalism, using the attribute _helAmp_. For example, suppose the
transversity amplitudes were used rather than the canonical, then the user can specify systems of helicity coupling:  
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
In this example, the longitudinal amplitude is the  <img src="figs/tex/fig14.png" style="margin-bottom:0px" /> 
helicity state, while the two transverse amplitudes and the sum and difference of the two other helicity amplitudes.

#### Fit Parameters and Expressions

Parameters can either be specified by three parameters, in the case of a scalar parameter such as a mass or a width, or
with six parameters in the case of a complex parameter such as a coupling. Upper and lower bounds on parameters can
also be set by specifying a parameter with five parameters or ten parameters for a scalar or complex, respectively.  
For example, if we wished to vary the mass of the  <img src="figs/tex/fig6.png" style="margin-bottom:-5px" />  meson in
the above example, but restricting the allowed values in the range  <img src="figs/tex/fig7.png" style="margin-bottom:-5px" /> :
```
K(1)(1270)+_mass 0 1.27 0.01 0.0 2.0
```
Parameters can also be related to each other via expressions.
Suppose for example we have  <img src="figs/tex/fig8.png" style="margin-bottom:-5px" />  
and  <img src="figs/tex/fig9.png" style="margin-bottom:-5px" />  in the same fit
(for example, for  <img src="figs/tex/fig10.png" style="margin-bottom:0px" /> )
The properties of one can be allowed to vary, for example the  <img src="figs/tex/fig11.png" style="margin-bottom:-5px" /> , 
and the other fixed to the same value, using:
```
K(1)(1270)+_mass 0 1.27 0.01 0.0 2.0
K(1)(1270)bar-_mass = K(1)(1270)+_mass
```

Parameter expressions are whitespace delimited due to the abundance of odd glyphs such as brackets and +/- in the names
of parameters. Expressions support the binary operations <img src="figs/tex/fig12.png" style="margin-bottom:-5px" />,
as well as common unary functions such as sqrt, trigonometric functions etc.

## Examples
Example options files can be found at `.$AMPGENROOT/options/`.

