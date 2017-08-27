# SiteResponse2D

A program to perform a non-linear site response analysis in 2 dimensions.

## Description

Given an input ground motion, a soil profile with properties of each layer a non-linear site response analysis can be carried out. Dispite its name, the program can also be used to perform the more conventional 1D response analysis. In this sense it is very similar to the [DEEPSOIL](http://deepsoil.cee.illinois.edu/) software. The extension to 2D was carried out using a circular yield function for the individual elasto-plastic springs (or layers) that make up the Masing rules. This is equivalent to the behaviour of MAT_HYSTERETIC_SOIL in the general purpose finite element software [LS-DYNA](http://www.lstc.com/products/ls-dyna). In fact, one obtains the *exact* same response for the undamped case, the damped response is (very) slightly different due to implementation differences.

The program calculates the response using the finite element method where the soil is discretised into many elements with depth. Use is made of the central difference explictit time integration scheme to calculate the response numerically. Only the stiffness-proportional part of Rayleigh damping is implemented, since the mass proportional part is deemed inappriate in this case. 

## Input file format

For the program to work the file name of the input file should be passed as an argument. It has several options that can be specified. See the example below for a fully defined input file. All of the other files are in the plain text comma (,) separated value (CSV) format. Units should be consistent, for formatting purposes it is assumed SI units are adopted. Ground motions should be provided as accelerograms, one column time in seconds, the other column acceleration in m/s^2. File names in the input file can also be (relative) paths. Description of the available options in the input file:

Parameter              | Description
:--------------------- | :----------
ground_motion_x        | File name of the ground motion in x-direction. [Required]
ground_motion_y        | File name of the ground motion in y-direction.
output_ground_motion_x | File name of the output ground motion in x-direction. [Required]
output_ground_motion_y | File name of the output ground motion in y-direction.
soil_elements          | File name of the soil elements file. This file contains information about the elements that are used to discretise the soil profile. See below for the specification of this format. [Required]
scale_factor           | Scale factor by which the input ground motion is multiplied. The default is 1.0 (i.e. no scaling).
damping_ratio          | Damping specified as ratio of critical; used to determine the stiffness-proportional part of Rayleigh damping. The default is 0.01 (i.e. 1% of critical). Note that this is supposed to represent the small strain damping; additional hysteretic damping is provided by the yielding of the material. 
damping_period         | Specifies the period at which the damping has the specified ratio. Note that due to the nature of stiffness-proportional damping, smaller periods have more damping and larger periods less. The default is 0.1 s.
damping_tangent        | If this parameter is set to 'true' the stiffness upon which the stiffness-proportional damping is based is adjusted as the soil material goes into a yielding state. This means that the effective damping becomes less when this happens. Note that LS-DYNA uses this behaviour. The default setting is 'true', but it can be set to 'false' if the damping should be based on the initial stiffness.
timestep               | Specifies timestep at which the calculations are performed. The default is 5e-4 which is a good compromise between accuracy and speed in most cases.
output_timestep        | Specifies timestep that is used for writing the output ground motion files. The default is 5e-3 s.
output_extent          | Can be used to add some additional time to the input motion so that if the soil is swaying in the end of the analysis it is allowed to come to a rest. It is specified in seconds and the default is 0 s.
verbosity              | Defines the level of output that is printed. The default is 1, specify 2 for extra information.

## Soil elements file format

The soil element file contains the elements that are used to model the soil. The format is a plain text comma (,) separated value (CSV) file. The first line contains the header, with the names of all parameters. A description is provided below.

Parameter&nbsp;&nbsp;&nbsp;&nbsp;| Description
:--------------------- | :----------
name                   | Each soil element can be given a name, this value is not used by the program.
z                      | The top of the element is given by the depth location (z). Values can be positive or negative but their absolute values should monotonically increasing.
rho                    | Density of soil element in kg/m^3
Vs                     | Shear wave velocity of soil element in m/s
strain value 1         | G_sec/G_0 value 1; The following columns define the secant shear stiffness (G_sec) divided by small-strain stiffness (G_0) curve. Provide as many points as needed, normally log-spaced.
strain value 2         | G_sec/G_0 value 2
strain value n         | G_sec/G_0 value n

The last row of the file defines the elastic half-space boundary that is implemented using the specification of the damper constant as the impedance (C = rho x Vs).

## Example

In the example folder of this repository an example is provided. It contains several files, including a fully specified input and soil elements file. For completeness, the contents of the input file are:

```
ground_motion_x = input_x.csv
ground_motion_y = input_y.csv
output_ground_motion_x = output_x.csv
output_ground_motion_y = output_y.csv
soil_elements = soil_elements.csv
scale_factor = 1.0
damping_ratio = 0.01
damping_period = 0.1
damping_tangent = true
timestep = 5e-4
output_timestep = 5e-3
output_extent = 0.0
verbosity = 1
```

The soil profile was discretised into 40 elements. For the top elements a finer discretisation is used than for the bottom ones. At a depth of 30 m the elastic half-space boundary is defined. The G_sec/G_0 curves have been discretised into 10 points, up to a strain of 1e-7 the elements behave elastic. In the following table the soil elements input file is displayed (note that some rows and columns are omitted).

name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|z|rho|Vs|1e-7|1e-6|1e-5|3e-5|1e-4|3e-4|0.001|0.003
:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--
Shallow clay|0|1650|100|1|0.998|0.979|0.939|0.832|0.623|0.341|0.181
Shallow clay|-0.5|1650|100|1|0.998|0.979|0.939|0.832|0.623|0.341|0.181
Shallow clay|-1|1650|100|1|0.998|0.979|0.939|0.832|0.623|0.341|0.181
Shallow clay|-1.5|1650|100|1|0.998|0.979|0.939|0.832|0.623|0.341|0.181
Deeper clay|-2|1650|101|1|0.998|0.982|0.948|0.85|0.66|0.392|0.212
Deeper clay|-2.5|1650|102.5|1|0.998|0.982|0.948|0.85|0.66|0.392|0.212
...|||||||||||
Deeper sand|-26|1950|276|1|0.997|0.969|0.917|0.782|0.566|0.302|0.137
Deeper sand|-27|1950|282|1|0.997|0.969|0.917|0.782|0.566|0.302|0.137
Deeper sand|-28|1950|288|1|0.997|0.969|0.917|0.782|0.566|0.302|0.137
Deeper sand|-29|1950|294|1|0.997|0.969|0.917|0.782|0.566|0.302|0.137
Boundary|-30|1950|300|||||||

When the example is run, it is seen that the non-linear modification is largest in the y-direction. In the following figures the response is shown as time-history and response spectrum. In the response spectrum plot the response of a 1D analysis is shown as well – to illustrate the effect of a 2D analysis.

![time-history plot](/example/plots/time-history.png)

![response spectrum plot](/example/plots/response_spectrum.png)

## Compiling the code

The source code can be compiled with any C++11 compatible compiler. In the source directory two folders are located (VS2015 and GNU), which contain the project file for Microsoft Visual Studio (Express) 2015 and a makefile for the GNU C++ compiler, respectively. For the latter `cd` into the GNU directory and type `make`, this will build the application for your platform.

Note that in the release directory builds are available for Windows and Linux (Ubuntu platform), both 64-bit.

## License

This software uses the [GNU GPLv3](/LICENSE.md) license. In general terms it means that this license requires anyone who distributes the code, or a derivative work, to make the source available under the same terms (i.e. open-source software).
