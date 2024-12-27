Torus fixed coils in simsopt
(cost function approach)
Elias G.

-------

In this repository you find prepared scripts for a simsopt coil optimizations calculation 
which takes a plasma surface vmec file as an input. 

In the current configuration of parameters you will perform a minimization forcing the coils onto a torus surface
using the plasma surface from https://quasr.flatironinstitute.org/model/0021326. 

It will also automatically create a poincare plot with the magnetic field of the coils.


-------- RUN -------

To run the code just run the #optimize.py# file or use the terminal:

    python optimize.py output_folder

By changing the parameters in the #config.py# file you can adjust your optimization parameters, cost function weights etc.

The #optimize.py# file  runs the files #main_optimization.py# and #run_poincare_simsopt.py# which can also be run separately.


------- terminal execution and flags -------

Using the terminal you can add different flags:

-o   ... just the optimization without producing the Poincare plots
-p   ... just the Poincare plots from the already optimized data in the result folder
-n   ... neither optimization nor Poincare plots will run
-orp ... use original parameters saved in the specified result folder from the prior optimization instead
    of the parameters in config.py.

-c [parameter] [value] ... change a parameter before optimization/poincare claculation


For example: You have already done an optimization and saved it in the folder ’output’. You have done
a few other optimizations in other folders since then and your config.py file has changed. Now you want to redo the optimization
from ’output’ with the original parameters (in the parameters.json file) but just change the current from 100 to 1000 ampere.

    python optimize.py output -o -orp -c CURRENT 1000

You look at the coils and want to do a Poincare plot with the plot parameters from parameters.json
in the folder:

    python optimize.py output -p

Using this framework you can conduct parameter searches with bash scripts in search of suitable coil
configurations.


Plasma:
https://quasr.flatironinstitute.org/model/0021326

https://github.com/itpplasma/ALPES/tree/main

https://github.com/hiddenSymmetries/simsopt

https://simsopt.readthedocs.io/en/latest/overview.html


