Fitting a straight line with non-symmetric uncertainties
========================================================

This repository is an exercise I had with `Iskren Georgiev`_ when he was
looking at linear correlations between the stellar mass of galaxies and the
properties of their central nuclear star cluster.

The approach here is neither frequentist nor fully Bayesian but more like in
between. In this approach, we sample the likelihood space instead of the prior
one to construct the posterior distribution. 
In addition, we explored how to quickly get fits through his data even when the
uncertainties were non-symmetric or upper limits.

The full mcmc approach is also included in the repository.

.. _Iskren Georgiev: http://www.mpia.de/homes/georgiev/

Example usage
-------------

The example below analyses the sample data from  `Iskren Georgiev`_'s paper. 
This run samples the likelihood using a split-normal distribution to account for
non-symmetric uncertainties, and bootstrap the data to include that the dataset
is relatively sparse.

.. code:: python

        ./fit.py reff_NSC_Mass_late.dat -o reff_NSC_Mass_late.theta.dat \
                -b -n 5000 --xfloor 10 --yfloor 10 \
                --log_xnorm 6 --log_ynorm 1 \
                --x12label '${\cal M}_{\rm NSC}$' \
                --y12label '$r_{\rm eff,NSC}$' \
                -f

.. image:: reff_NSC_Mass_late.png


Command line details
--------------------

Options are all available from the command line and from the (optional)
configuration file.


Options
~~~~~~~

.. code:: bash

        ./fit --help for options

**fitting options**

+-------------------------+-------------------------------------------------------------+
|  -n, --nsamp            |  number of samples to represent per data point uncertainties|
+-------------------------+-------------------------------------------------------------+
|  -b, --bootstrap        |  include bootstrapping from the dataset                     |
+-------------------------+-------------------------------------------------------------+
|  --log_xnorm            |  x-data normalization value                                 |
+-------------------------+-------------------------------------------------------------+
|  --log_ynorm            |  y-data normalization value                                 |
+-------------------------+-------------------------------------------------------------+
|  --xfloor               |  floor of x-value uncertainty (in %)                        |
+-------------------------+-------------------------------------------------------------+
|  --yfloor               |  floor of y-value uncertainty (in %)                        |
+-------------------------+-------------------------------------------------------------+

**output options**

+-------------------------+-------------------------------------------------------------+
|  -o OUTPUT, --output    |  export the samples into a file                             |
+-------------------------+-------------------------------------------------------------+
|  -f, --savefig          |  Generate figures with the desired format (pdf, png...)     |
+-------------------------+-------------------------------------------------------------+

**plotting options**

+------------------------+--------------------------------------------------------------+
|  --xlabel=XLABEL       |   X-label of the top-right plot (it can be in latex form)    |
+------------------------+--------------------------------------------------------------+
|  --ylabel=YLABEL       |   Y-label of the top-right plot (it can be in latex form)    |
+------------------------+--------------------------------------------------------------+
                                                                                                 
**special options**                                                                             

+------------------------+------------------------------------------------------------------------+
|  -c, --config          |   Configuration file to use for default values (see below for details) |
+------------------------+------------------------------------------------------------------------+
                                                                                                 

Configuration file
~~~~~~~~~~~~~~~~~~
All options from the command line have default values. Theses values are
arbitrary but can be set by using a configuration file (option: `-c`).

In this file, any option can be set but the **command line has priority**. To
set one option, you need to use the *longname* option as reference.
