NOTE: 
'$' symbol represents the text following the symbol has to be typed/ pasted on to the terminal. WITHOUT THE '$'. 

Case study: Calcite dissolution considering kinetic equilibrium - with no porosity feedback

Phreeqc test setup:
The phreeqc run file is available in 'constant/phreeqcRun.phr' of folder - dt0_75

To run the phreeqc file:
$phreeqc phreeqcRun.phr phreeqcRun.log

The output of the phreeqc run is written to a file named 'CalciteDiss_KinEq.op' 

The output data used in the plots is available in 'plots/phreeqcDataExtraction' which has been extracted from 'CalciteDiss_KinEq.op'

If user changes the phreeqc input file (i.e. constant/phreeqcRun.phr) the user need to run the following to update data in the plots
$./plot/phreeqcDataExtraction/runFile

Corresponding setup case investigated using pM4F is available in:
(a) dt = 0.75s
