# MPT\_Surr #

MPT\_Surr is a MATLAB code used to calculate limits in composition for transportation fuel mixtures based on experimental uncertainty using several tools from optimization theory. It incorporates MW, H/C ratio, TSI (threshold sooting index), ignition delay time and distillation curve error norm for the time being and is very amenable to utilizing other constraints.

## How do I get set up? ##

Once [pyIDT](https://github.com/gpavanb/pyIDT) has been used to create the search space, the current code can be used to do the rest. The details are all described in the accompanying publication [1]. 

There are only 3 files which you will need to work with
* `construct_linear` : For MW and H/C constraints only
* `construct_nonlinear` : For both linear and non-linear constraints
* `construct_mpt` : For multi-parametric optimization to obtain bounds as function of error threshold (linear constraints only)

Just modify the `INPUTS` section in the beginning of each of the file to your prescribed surrogates. You might need to modify `set_target_details.m` in the include directory if you want to try out your own surrogate. Obviously, distillation and IDT data have to be computed in such a situation.

### Prerequisites ###
This code goes along with the [pyIDT](https://github.com/gpavanb/pyIDT) package, which is used to create the search space based on ignition delay time. 

All the required packages are based in MATLAB and can be easily installed by downloading the files and adding them to the search path. Please refer to the links provided alongside each package for further instructions regarding installation.

* [cantera - MATLAB](http://www.cantera.org/docs/sphinx/html/install.html)
* [CVX](http://cvxr.com/cvx/)
* [MPT Toolbox](http://control.ee.ethz.ch/~mpt/3/Main/Installation)
* [Ellipsoidal Toolbox](https://code.google.com/archive/p/ellipsoids/downloads) : This works only with Version 1.1.3!
* [0d\_multicomp](https://bitbucket.org/gpavanb/0d_multicomp)
* [DrawCuboid](https://www.mathworks.com/matlabcentral/fileexchange/25559-draw-cuboid)

## License ##
Please refer to the LICENSE.pdf in the repository. Note that this code requires PRIOR PERMISSION FROM AUTHORS FOR COMMERCIAL PURPOSES.

## Who do I talk to? ##

* Repo owner or admin : [Pavan Bharadwaj](https://github.com/gpavanb)
* Other community or team contact : The code was developed at the Flow Physics and Computational Engineering group at Stanford University. Please direct any official queries to [Prof. Matthias Ihme](mailto:mihme@stanford.edu)

## References ##
[1] P.B. Govindaraju, M. Ihme, "Formulation of Optimal Surrogate Descriptions of Fuels Considering Sensitivities to Experimental Uncertainties", Combustion and Flame, Vol. 188, Pages 337-356

Several software packages have also been utilized for this code and their contribution is acknowledged in the 'References' section of [1]. 
