
Case study: Calcite dissolution considering thermodynamic equilibrium - with no porosity feedback

test case described in Pavuluri et al. (2022?)

Phreeqc test setup:
The phreeqc run file is available in 'constant/phreeqcRun.phr'

To run the phreeqc file:
$phreeqc phreeqcRun.phr phreeqcRun.log

The output of the phreeqc run is written to a file named 'CalciteDiss_ThermoEq.op' 

The output data used in the plots is available in 'plots/phreeqcDataExtraction' which has been extracted from 'CalciteDiss_ThermoEq.op'

If user changes the phreeqc input file (i.e. constant/phreeqcRun.phr) the user need to run the following to update data in the plots
$./plot/phreeqcDataExtraction/runFile
