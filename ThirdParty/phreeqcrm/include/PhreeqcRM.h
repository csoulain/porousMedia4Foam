/*! @file PhreeqcRM.h
	@brief C++ Documentation
*/
#if !defined(PHREEQCRM_H_INCLUDED)
#define PHREEQCRM_H_INCLUDED
#ifdef USE_MPI
#include "mpi.h"
#define MP_TYPE MPI_Comm
#else
#define MP_TYPE int
#endif
class IPhreeqcPhast;

class cxxStorageBin;
//class cxxNameDouble;
#include "NameDouble.h"
class cxxSolution;
class PHRQ_io;
#include <vector>
#include <list>
#include <set>
#include <map>
#include <string>

#if defined(_WINDLL)
#define IRM_DLL_EXPORT __declspec(dllexport)
#else
#define IRM_DLL_EXPORT
#endif

class PHRQ_io;
class IPhreeqc;
/**
 * @class PhreeqcRMStop
 *
 * @brief This class is derived from std::exception and is thrown
 * when an unrecoverable error has occurred.
 */
class IRM_DLL_EXPORT PhreeqcRMStop : public std::exception
{
public:
  const char *what() const throw () {return "Failure in PhreeqcRM\n";}
};

/*! @brief Enumeration used to return error codes.
*/
#include "IrmResult.h"
enum {
	METHOD_CREATEMAPPING,
	METHOD_DUMPMODULE,
	METHOD_FINDCOMPONENTS,
	METHOD_GETCONCENTRATIONS,
	METHOD_GETDENSITY,
	METHOD_GETERRORSTRING,
	METHOD_GETPRESSURE,
	METHOD_GETSATURATION,
	METHOD_GETSELECTEDOUTPUT,
	METHOD_GETSOLUTIONVOLUME,
	METHOD_GETSPECIESCONCENTRATIONS,
	METHOD_GETSPECIESLOG10GAMMAS,
	METHOD_GETTEMPERATURE,
	METHOD_INITIALPHREEQC2MODULE,
	METHOD_INITIALPHREEQCCELL2MODULE,
	METHOD_LOADDATABASE,
	METHOD_MPIWORKERBREAK,
	METHOD_RUNCELLS,
	METHOD_RUNFILE,
	METHOD_RUNSTRING,
	METHOD_SETCOMPONENTH2O,
	METHOD_SETCONCENTRATIONS,
	METHOD_SETDENSITY,
	METHOD_SETERRORHANDLERMODE,
	METHOD_SETFILEPREFIX,
	METHOD_SETPARTITIONUZSOLIDS,
	METHOD_SETPOROSITY,
	METHOD_SETPRESSURE,
	METHOD_SETPRINTCHEMISTRYON,
	METHOD_SETPRINTCHEMISTRYMASK,
	METHOD_SETREBALANCEBYCELL,
	METHOD_SETREBALANCEFRACTION,
	METHOD_SETREPRESENTATIVEVOLUME,
	METHOD_SETSATURATION,
	METHOD_SETSELECTEDOUTPUTON,
	METHOD_SETSPECIESSAVEON,
	METHOD_SETTEMPERATURE,
	METHOD_SETTIME,
	METHOD_SETTIMECONVERSION,
	METHOD_SETTIMESTEP,
	METHOD_SETUNITSEXCHANGE,
	METHOD_SETUNITSGASPHASE,
	METHOD_SETUNITSKINETICS,
	METHOD_SETUNITSPPASSEMBLAGE,
	METHOD_SETUNITSSOLUTION,
	METHOD_SETUNITSSSASSEMBLAGE,
	METHOD_SETUNITSSURFACE,
	METHOD_SPECIESCONCENTRATIONS2MODULE,
	METHOD_USESOLUTIONDENSITYVOLUME
} /* MPI_METHOD */;

/**
 * @mainpage PhreeqcRM Library Documentation (3.6.2-15100)
 *
 *  @htmlonly
 *  <table>
 *   <tr><td class="indexkey"><a class="el" href="class_phreeqc_r_m.html">PhreeqRM.h</a> </td><td class="indexvalue">C++ Documentation</td></tr>
 *   <tr><td class="indexkey"><a class="el" href="_r_m__interface___c_8h.html">RM_interface_C.h</a> </td><td class="indexvalue">C Documentation </td></tr>
 *   <tr><td class="indexkey"><a class="el" href="classphreeqcrm.html">RM_interface.F90</a></td><td class="indexvalue">Fortran Documentation </td></tr>
 *   <tr><td class="indexkey"><a class="el" href="_irm_result_8h.html">IrmResult.h</a></td><td class="indexvalue">Return codes </td></tr>
 *  </table>
 *  @endhtmlonly
 */

/**
 * @class PhreeqcRM
 *
 * @brief Geochemical reaction module
 */


class IRM_DLL_EXPORT PhreeqcRM
{
public:
	static void             CleanupReactionModuleInstances(void);
	static int              CreateReactionModule(int nxyz, MP_TYPE nthreads);
	static IRM_RESULT       DestroyReactionModule(int n);
	static PhreeqcRM      * GetInstance(int n);

/**
Constructor for the PhreeqcRM reaction module. If the code is compiled with
the preprocessor directive USE_OPENMP, the reaction module use OPENMP and multiple threads.
If the code is compiled with the preprocessor directive USE_MPI, the reaction
module will use MPI and multiple processes. If neither preprocessor directive is used,
the reaction module will be serial (unparallelized).
@param nxyz        The number of grid cells in the users model.
@param thread_count_or_communicator If multithreaded, the number of threads to use
in parallel segments of the code.
If @a thread_count_or_communicator is <= 0, the number of threads is set equal to the number of processors in the computer.
If multiprocessor, the MPI communicator to use within the reaction module.
@param io        Optionally, a PHRQ_io input/output object can be provided to the constructor. By default
a PHRQ_io object is constructed to handle reading and writing files.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int nxyz = 40;
#ifdef USE_MPI
  PhreeqcRM phreeqc_rm(nxyz, MPI_COMM_WORLD);
  int mpi_myself;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
  {
    exit(4);
  }
  if (mpi_myself > 0)
  {
    phreeqc_rm.MpiWorker();
    return EXIT_SUCCESS;
  }
#else
  int nthreads = 3;
  PhreeqcRM phreeqc_rm(nxyz, nthreads);
#endif
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and all workers.
 */
	PhreeqcRM(int nxyz, MP_TYPE thread_count_or_communicator, PHRQ_io * io=NULL);
	~PhreeqcRM(void);
/**
Close the output and log files.
@retval IRM_RESULT   0 is success, negative is failure (See @ref DecodeError).
@see                 @ref OpenFiles, @ref SetFilePrefix
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.CloseFiles();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called only by root.
 */
	IRM_RESULT                                CloseFiles(void);
/**
@a N sets of component concentrations are converted to SOLUTIONs numbered 1-@a n in the Utility IPhreeqc.
The solutions can be reacted and manipulated with the methods of IPhreeqc. If solution concentration units
(@ref SetUnitsSolution) are per liter, one liter of solution is created in the Utility instance; if solution
concentration units are mass fraction, one kilogram of solution is created in the Utility instance.
The motivation for this
method is the mixing of solutions in wells, where it may be necessary to calculate solution properties
(pH for example) or react the mixture to form scale minerals.
The code fragment below makes a mixture of
concentrations and then calculates the pH of the mixture.
@param c             Vector of concentrations to be made SOLUTIONs in Utility IPhreeqc.
Vector contains @a n values for each component (@ref GetComponentCount) in sequence.
@param tc            Vector of temperatures to apply to the SOLUTIONs, in degrees C. Vector of size @a n.
@param p_atm         Vector of pressures to apply to the SOLUTIONs, in atm. Vector of size @a n.
@retval IRM_RESULT   0 is success, negative is failure (See @ref DecodeError).
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector <double> c_well;
c_well.resize(1*ncomps, 0.0);
for (int i = 0; i < ncomps; i++)
{
  c_well[i] = 0.5 * c[0 + nxyz*i] + 0.5 * c[9 + nxyz*i];
}
std::vector<double> tc, p_atm;
tc.resize(1, 15.0);
p_atm.resize(1, 3.0);
IPhreeqc * util_ptr = phreeqc_rm.Concentrations2Utility(c_well, tc, p_atm);
input = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1";
int iphreeqc_result;
util_ptr->SetOutputFileName("utility_cpp.txt");
util_ptr->SetOutputFileOn(true);
iphreeqc_result = util_ptr->RunString(input.c_str());
phreeqc_rm.ErrorHandler(iphreeqc_result, "IPhreeqc RunString failed");
int vtype;
double pH;
char svalue[100];
util_ptr->SetCurrentSelectedOutputUserNumber(5);
iphreeqc_result = util_ptr->GetSelectedOutputValue2(1, 0, &vtype, &pH, svalue, 100);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called only by root.
 */
	IPhreeqc * Concentrations2Utility(std::vector<double> &c,
		   std::vector<double> tc, std::vector<double> p_atm);
/**
Provides a mapping from grid cells in the user's model to reaction cells for which chemistry needs to be run.
The mapping is used to eliminate inactive cells and to use symmetry to decrease the number of cells
for which chemistry must be run. 
The array @a grid2chem of size @a nxyz (the number of grid cells, @ref GetGridCellCount) 
must contain the set of all integers 0 <= @a i < @a count_chemistry, 
where @a count_chemistry is a number less than or equal to @a nxyz.
Inactive cells are assigned a negative integer. 
The mapping may be many-to-one to account for symmetry.
Default is a one-to-one mapping--all user grid cells are reaction cells
(equivalent to @a grid2chem values of 0,1,2,3,...,nxyz-1).
@param grid2chem        A vector of integers: Nonnegative is a reaction-cell number (0 based),
negative is an inactive cell. Vector is of size @a nxyz (number of grid cells, @ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// For demonstation, two equivalent rows by symmetry
std::vector<int> grid2chem;
grid2chem.resize(nxyz, -1);
for (int i = 0; i < nxyz/2; i++)
{
  grid2chem[i] = i;
  grid2chem[i + nxyz/2] = i;
}
status = phreeqc_rm.CreateMapping(grid2chem);
if (status < 0) phreeqc_rm.DecodeError(status);
int nchem = phreeqc_rm.GetChemistryCellCount();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                CreateMapping(std::vector<int> &grid2chem);
/**
If @a result is negative, this method prints an error message corresponding to IRM_RESULT @a result.
If @a result is non-negative, no action is taken.
@param result               An IRM_RESULT value returned by one of the reaction-module methods.
@par IRM_RESULT definition:
@htmlonly
<CODE>
<PRE>
typedef enum {
  IRM_OK            =  0,  //Success
  IRM_OUTOFMEMORY   = -1,  //Failure, Out of memory
  IRM_BADVARTYPE    = -2,  //Failure, Invalid VAR type
  IRM_INVALIDARG    = -3,  //Failure, Invalid argument
  IRM_INVALIDROW    = -4,  //Failure, Invalid row
  IRM_INVALIDCOL    = -5,  //Failure, Invalid column
  IRM_BADINSTANCE   = -6,  //Failure, Invalid rm instance id
  IRM_FAIL          = -7,  //Failure, Unspecified
} IRM_RESULT;
</PRE>
</CODE>
@endhtmlonly
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.CreateMapping(grid2chem);
phreeqc_rm.DecodeError(status);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Can be called by root and (or) workers.
 */
	void                                      DecodeError(int result);
/**
Writes the contents of all workers to file in _RAW formats (see appendix of PHREEQC version 3 manual),
including SOLUTIONs and all reactants.
@param dump_on          Signal for writing the dump file, true or false.
@param append           Signal to append to the contents of the dump file, true or false.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetDumpFileName.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
bool dump_on = true;
bool append = false;
status = phreeqc_rm.SetDumpFileName("advection_cpp.dmp");
status = phreeqc_rm.DumpModule(dump_on, append);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root; workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                DumpModule(bool dump_on, bool append = false);
/**
Checks @a result for an error code. If result is negative, the result is decoded (@ref DecodeError),
and printed as an error message along with the @a e_string, and an exception is thrown. If the result
is nonnegative, no action is taken.
@param result           IRM_RESULT to be checked for an error.
@param e_string         String to be printed if an error is found.
@see                    @ref DecodeError, @ref ErrorMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
iphreeqc_result = util_ptr->RunString(input.c_str());
if (iphreeqc_result != 0)
{
  phreeqc_rm.ErrorHandler(IRM_FAIL, "IPhreeqc RunString failed");
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	void                                      ErrorHandler(int result, const std::string &e_string);
/**
Send an error message to the screen, the output file, and the log file.
@param error_string      String to be printed.
@param prepend           True, prepends @a error_string with "Error: "; false, @a error_string is used with no prepended text.
@see                    @ref OpenFiles, @ref LogMessage, @ref ScreenMessage, @ref WarningMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
phreeqc_rm.ErrorMessage("Goodbye world");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers; root writes to output and log files.
 */
	void                                      ErrorMessage(const std::string &error_string, bool prepend = true);
/**
This method accumulates a list of elements. Elements are those that have been
defined in a solution or any other reactant
(EQUILIBRIUM_PHASE, KINETICS, and others), including charge imbalance.
This method can be called multiple times and the list that is created is cummulative.
The list is the set of components that needs to be transported. By default the list
includes water, excess H and excess O (the H and O not contained in water);
alternatively, the list may be set to contain total H and total O (@ref SetComponentH2O),
which requires transport results to be accurate to eight or nine significant digits.
If multicomponent diffusion (MCD) is to be modeled,
there is a capability to retrieve aqueous species concentrations
(@ref GetSpeciesConcentrations) and to set new solution concentrations after
MCD by using individual species concentrations
(@ref SpeciesConcentrations2Module).
To use these methods the save-species property needs to be turned on (@ref SetSpeciesSaveOn).
If the save-species property is on, FindComponents will generate
a list of aqueous species (@ref GetSpeciesCount, @ref GetSpeciesNames), 
their diffusion coefficients at 25 C (@ref GetSpeciesD25),
and their charge (@ref GetSpeciesZ).
@retval              Number of components currently in the list, or IRM_RESULT error code 
(negative value, see @ref DecodeError).
@see                   @ref GetComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount, 
@ref GetSpeciesD25, 
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SetComponentH2O,
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module. 
@par The @ref FindComponents method also generates lists of reactants--equilibrium phases,
exchangers, gas components, kinetic reactants, solid solution components, and surfaces. 
The lists are cumulative, including all reactants that were
defined in the initial phreeqc instance at any time FindComponents was called.
In addition, a list of phases is generated for which saturation indices may be calculated from the
cumulative list of components.
@see also
@ref GetEquilibriumPhases,
@ref GetEquilibriumPhasesCount,
@ref GetExchangeNames,
@ref GetExchangeSpecies,
@ref GetExchangeSpeciesCount,
@ref GetGasComponents,
@ref GetGasComponentsCount,
@ref GetKineticReactions,
@ref GetKineticReactionsCount,
@ref GetSICount,
@ref GetSINames,
@ref GetSolidSolutionComponents,
@ref GetSolidSolutionComponentsCount,
@ref GetSolidSolutionNames,
@ref GetSurfaceNames,
@ref GetSurfaceSpecies,
@ref GetSurfaceSpeciesCount,
@ref GetSurfaceTypes.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int ncomps = phreeqc_rm.FindComponents();
const std::vector<std::string> &components = phreeqc_rm.GetComponents();
const std::vector < double > & gfw = phreeqc_rm.GetGfw();
for (int i = 0; i < ncomps; i++)
{
  std::ostringstream strm;
  strm.width(10);
  strm << components[i] << "    " << gfw[i] << "\n";
  phreeqc_rm.OutputMessage(strm.str());
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	int                                       FindComponents();
/**
Returns a vector of vectors,
where the @a nth vector is a vector of grid-cell numbers
that are mapped to reaction-cell number @a n.
Each reaction-cell number has a vector of one or more grid-cell numbers.
@retval              Vector of vectors of ints. For each reaction cell @a n,
the @a nth vector in the vector of vectors contains
the grid-cell numbers that map to the reaction cell.
@see                 @ref CreateMapping, @ref GetForwardMapping, @ref GetChemistryCellCount, 
@ref GetGridCellCount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector < std::vector <int> > & back = phreeqcrm_ptr->GetBackwardMapping();
if (option == "HYDRAULIC_K")
{
  return (*data_ptr->hydraulic_K)[back[rm_cell_number][0]];
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root or workers.
 */
	const std::vector < std::vector <int> > & GetBackwardMapping(void) {return this->backward_mapping;}
/**
Returns the number of reaction cells in the reaction module. The number of reaction cells is defined by
the set of non-negative integers in the mapping from grid cells (@ref CreateMapping), or, by default,
the number of grid cells (@ref GetGridCellCount).
The number of reaction cells is less than or equal to the number of grid cells in the user's model.
@retval              Number of reaction cells.
@see                 @ref CreateMapping, @ref GetGridCellCount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.CreateMapping(grid2chem);
std::ostringstream oss;
oss << "Number of reaction cells in the reaction module: "
    << phreeqc_rm.GetChemistryCellCount() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetChemistryCellCount(void) const {return this->count_chemistry;}
/**
Returns the number of components in the reaction-module component list.
@retval                 The number of components in the reaction-module component list. The component list is
generated by calls to @ref FindComponents.
The return value from the last call to @ref FindComponents is equal to the return value from GetComponentCount.
@see                    @ref FindComponents, @ref GetComponents.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Number of components for transport: " << phreeqc_rm.GetComponentCount() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	int                                       GetComponentCount(void) const {return (int) this->components.size();}
/**
Returns a reference to the reaction-module component list that was generated by calls to @ref FindComponents.
@retval const std::vector<std::string>&       A vector of strings; each string is a component name.
@see                    @ref FindComponents, @ref GetComponentCount
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<std::string> &components = phreeqc_rm.GetComponents();
const std::vector < double > & gfw = phreeqc_rm.GetGfw();
int ncomps = phreeqc_rm.GetComponentCount();
for (int i = 0; i < ncomps; i++)
{
  std::ostringstream strm;
  strm.width(10);
  strm << components[i] << "    " << gfw[i] << "\n";
  phreeqc_rm.OutputMessage(strm.str());
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
const std::vector<std::string> &          GetComponents(void) const {return this->components;}
	
/**
Transfer solution concentrations from each reaction cell
to the concentration vector given in the argument list (@a c).
Units of concentration for @a c are defined by @ref SetUnitsSolution.
For per liter concentration units,
solution volume is used to calculate the concentrations for @a c.
For mass-fraction concentration units, the solution mass is used to calculate concentrations for @a c.
Two options are available for the volume and mass of solution
that are used in converting to transport concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or (2) the volume of solution is the product of saturation (@ref SetSaturation),
porosity (@ref SetPorosity), and representative volume (@ref SetRepresentativeVolume),
and the mass of solution is volume times density as defined by @ref SetDensity.
@ref UseSolutionDensityVolume determines which option is used.
For option 1, the databases that have partial molar volume definitions needed
to accurately calculate solution volume are
phreeqc.dat, Amm.dat, and pitzer.dat.

@param c                Vector to receive the concentrations.
Dimension of the vector is set to @a ncomps times @a nxyz,
where,  ncomps is the result of @ref FindComponents or @ref GetComponentCount,
and @a nxyz is the number of user grid cells (@ref GetGridCellCount).
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref FindComponents, @ref GetComponentCount, @ref GetSaturation, @ref SetConcentrations,
@ref SetDensity, @ref SetRepresentativeVolume, @ref SetSaturation, @ref SetUnitsSolution, @ref UseSolutionDensityVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> c;
status = phreeqc_rm.RunCells();
status = phreeqc_rm.GetConcentrations(c);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                GetConcentrations(std::vector<double> &c);
/**
Returns the file name of the database. Should be called after @ref LoadDatabase.
@retval std::string      The file name defined in @ref LoadDatabase.
@see                    @ref LoadDatabase.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Database: " << phreeqc_rm.GetDatabaseFileName() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	std::string                               GetDatabaseFileName(void) {return this->database_file_name;}
/**
Transfer solution densities from the reaction-module workers to the vector given in the argument list (@a density).
@param density              Vector to receive the densities. Dimension of the array is set to @a nxyz,
where @a nxyz is the number of user grid cells (@ref GetGridCellCount).
Values for inactive cells are set to 1e30.
Densities are those calculated by the reaction module.
Only the following databases distributed with PhreeqcRM have molar volume information needed
to accurately calculate density: phreeqc.dat, Amm.dat, and pitzer.dat.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetSolutionVolume, @ref SetDensity.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.RunCells();
status = phreeqc_rm.GetConcentrations(c);
std::vector<double> density;
status = phreeqc_rm.GetDensity(density);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                GetDensity(std::vector<double> & density);
/**
Returns a vector of integers that contains the largest reaction-cell number assigned to each worker.
Each worker is assigned a range of reaction-cell numbers that are run during a call to @ref RunCells.
The range of reaction cells for a worker may vary as load rebalancing occurs.
At any point in the calculations, the first cell and last cell to be run by a worker can be found
in the vectors returned by @ref GetStartCell and @ref GetEndCell.
Each method returns a vector of integers that has length of the number of threads (@ref GetThreadCount),
if using OPENMP, or the number of processes (@ref GetMpiTasks), if using MPI.

@retval IRM_RESULT      Vector of integers, one for each worker, that gives the last reaction cell
to be run by each worker.
@see                    @ref GetStartCell, @ref GetThreadCount, @ref GetMpiTasks.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Current distribution of cells for workers\n";
oss << "Worker First Cell   Last Cell\n";
int n;
n = phreeqc_rm.GetThreadCount() * phreeqc_rm.GetMpiTasks();
for (int i = 0; i < n; i++)
{
	oss << i << "      "
	    << phreeqc_rm.GetStartCell()[i]
	    << "            "
		<< phreeqc_rm.GetEndCell()[i] << "\n";
}
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector < int> &                GetEndCell(void) const {return this->end_cell;}
/**
Returns a reference to the vector of all equilibrium phases.
The list includes all phases included in any EQUILIBRIUM_PHASES definitions in
the initial-phreeqc module.
@ref FindComponents must be called before @ref GetEquilibriumPhases.
This method may be useful when generating selected output definitions related to equilibrium phases.

@retval const std::vector<std::string>&       A vector of strings; each string is a unique
equilibrium phases name.

@see                    @ref FindComponents,
@ref GetEquilibriumPhasesCount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -equilibrium_phases " << "\n";
// equilibrium phases
const std::vector<std::string> &eq_phases = phreeqc_rm.GetEquilibriumPhases();
for (size_t i = 0; i < phreeqc_rm.GetEquilibriumPhasesCount(); i++)
{
oss << "    " << eq_phases[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetEquilibriumPhases(void) const { return this->EquilibriumPhasesList; }
/**
Returns the number of equilibrium phases in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetEquilibriumPhasesCount.
This method may be useful when generating selected output definitions related to
equilibrium phases.

@retval                 The number of equilibrium phases in the initial-phreeqc module.

@see                    @ref FindComponents,
@ref GetEquilibriumPhases.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -equilibrium_phases " << "\n";
// equilibrium phases
const std::vector<std::string> &eq_phases = phreeqc_rm.GetEquilibriumPhases();
for (size_t i = 0; i < phreeqc_rm.GetEquilibriumPhasesCount(); i++)
{
oss << "    " << eq_phases[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
int                                       GetEquilibriumPhasesCount(void) const { return (int) this->EquilibriumPhasesList.size(); }


/**
Get the setting for the action to be taken when the reaction module encounters an error.
Options are 0, return to calling program with an error return code (default);
1, throw an exception, which can be caught in C++ (for C and Fortran, the program will exit);
2, attempt to exit gracefully.
@retval IRM_RESULT      Current setting for the error handling mode: 0, 1, or 2.
@see                    @ref SetErrorHandlerMode.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Error handler mode: " << phreeqc_rm.GetErrorHandlerMode() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetErrorHandlerMode(void) {return this->error_handler_mode;}
/**
Returns a standard string containing error messages related to the last call to a PhreeqcRM method.
@retval                 Error messages related to the last call to a PhreeqcRM method.
@see                    @ref GetErrorHandlerMode.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
if (status != IRM_OK)
{
  std::cerr << phreeqc_rm.GetErrorString();
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	std::string                                  GetErrorString(void);

/**
Returns a reference to the vector of exchange names (such as "X") that correspond with
the exchange species names.
@ref FindComponents must be called before @ref GetExchangeNames.
The exchange names vector is the same length as the exchange species names vector
and provides the corresponding exchange site.
This method may be useful when generating selected output definitions related to exchangers.

@retval const std::vector<std::string>&       A vector of strings; each string is an
exchange name corresponding to the exchange species vector; an exchange name may occur
multiple times.

@see                    @ref FindComponents,
@ref GetExchangeSpeciesCount, @ref GetExchangeSpecies.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// molalities of exchange species
const std::vector<std::string> &ex_species = phreeqc_rm.GetExchangeSpecies();
const std::vector<std::string> &ex_names = phreeqc_rm.GetExchangeNames();
for (size_t i = 0; i < phreeqc_rm.GetExchangeSpeciesCount(); i++)
{

oss << "    ";
oss.width(15);
oss << std::left << ex_species[i];
oss << " # " << ex_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetExchangeNames(void) const { return this->ExchangeNamesList; }
/**
Returns a reference to the vector of exchange species names (such as "NaX").
The list of exchange species (such as "NaX") is derived from the list of components
(@ref FindComponents) and the list of all exchange names (such as "X")
that are included in EXCHANGE definitions in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetExchangeSpecies.
This method may be useful when generating selected output definitions related to exchangers.

@retval const std::vector<std::string>&       A vector of strings; each string is a
unique exchange species name.

@see                    @ref FindComponents,
@ref GetExchangeSpeciesCount, @ref GetExchangeNames.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// molalities of exchange species
const std::vector<std::string> &ex_species = phreeqc_rm.GetExchangeSpecies();
const std::vector<std::string> &ex_names = phreeqc_rm.GetExchangeNames();
for (size_t i = 0; i < phreeqc_rm.GetExchangeSpeciesCount(); i++)
{

oss << "    ";
oss.width(15);
oss << std::left << ex_species[i];
oss << " # " << ex_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetExchangeSpecies(void) const { return this->ExchangeSpeciesNamesList; }
/**
Returns the number of exchange species in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetExchangeSpeciesCount.
This method may be useful when generating selected output definitions related to exchangers.

@retval                 The number of exchange species in the initial-phreeqc module.

@see                    @ref FindComponents,
@ref GetExchangeSpecies, @ref GetExchangeNames.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
// molalities of exchange species
const std::vector<std::string> &ex_species = phreeqc_rm.GetExchangeSpecies();
const std::vector<std::string> &ex_names = phreeqc_rm.GetExchangeNames();
for (size_t i = 0; i < phreeqc_rm.GetExchangeSpeciesCount(); i++)
{

oss << "    ";
oss.width(15);
oss << std::left << ex_species[i];
oss << " # " << ex_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
int                                       GetExchangeSpeciesCount(void) const { return (int) this->ExchangeSpeciesNamesList.size(); }


/**
Returns the file prefix for the output (.chem.txt) and log files (.log.txt).
@retval std::string     The file prefix as set by @ref SetFilePrefix, or "myrun", by default.
@see                    @ref SetFilePrefix.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "File prefix: " << phreeqc_rm.GetFilePrefix() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	std::string                               GetFilePrefix(void) {return this->file_prefix;}
/**
Returns a reference to a vector of ints that is a mapping from grid cells to
reaction cells.
The mapping is used to eliminate cells that are inactive and cells that are unnecessary because of symmetry
from the list of cells for which reactions must be run. The mapping may be many-to-one to account for symmetry.
The mapping is set by @ref CreateMapping, or, by default, is a one-to-one mapping--all grid cells are
reaction cells (vector contains 0,1,2,3,...,@a nxyz-1).
@retval const std::vector < int >&      A vector of integers of size @a nxyz (number of grid cells, @ref GetGridCellCount).
Nonnegative is a reaction-cell number (0 based), negative is an inactive cell.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<int> &f_map = phreeqc_rm.GetForwardMapping();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
const std::vector < int > &               GetForwardMapping(void) {return this->forward_mapping_root;}

/**
Returns a reference to the vector of all gas components in the initial-phreeqc module.
The list includes all gas components included in any GAS_PHASE definitions in
the initial-phreeqc module.
@ref FindComponents must be called before @ref GetGasComponents.
This method may be useful when generating selected output definitions related to gas phases.

@retval const std::vector<std::string>&       A vector of strings; each string is a unique
gas component name.

@see                    @ref FindComponents,
@ref GetGasComponentsCount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -gases " << "\n";
// gas components
const std::vector<std::string> &gas_phases = phreeqc_rm.GetGasComponents();
for (size_t i = 0; i < phreeqc_rm.GetGasComponentsCount(); i++)
{
oss << "    " << gas_phases[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetGasComponents(void) const { return this->GasComponentsList; }
/**
Returns the number of gas phase components in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetGasComponentsCount.
This method may be useful when generating selected output definitions related to
gas phases.

@retval                 The number of gas phase components in the initial-phreeqc module.

@see                    @ref FindComponents,
@ref GetGasComponents.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -gases " << "\n";
// gas components
const std::vector<std::string> &gas_phases = phreeqc_rm.GetGasComponents();
for (size_t i = 0; i < phreeqc_rm.GetGasComponentsCount(); i++)
{
oss << "    " << gas_phases[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
int                                       GetGasComponentsCount(void) const { return (int) this->GasComponentsList.size(); }

/**
Returns a reference to a vector of doubles that contains the gram-formula weight of
each component. Called after @ref FindComponents. Order of weights corresponds to the list of components from
@ref GetComponents.
@retval const std::vector<double>&       A vector of doubles; each value is a component gram-formula weight, g/mol.
@see                    @ref FindComponents, @ref GetComponents.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<std::string> &components = phreeqc_rm.GetComponents();
const std::vector < double > & gfw = phreeqc_rm.GetGfw();
int ncomps = phreeqc_rm.GetComponentCount();
for (int i = 0; i < ncomps; i++)
{
  std::ostringstream strm;
  strm.width(10);
  strm << components[i] << "    " << gfw[i] << "\n";
  phreeqc_rm.OutputMessage(strm.str());
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector < double > &            GetGfw(void) {return this->gfw;}
/**
Returns the number of grid cells in the user's model, which is defined in
the call to the constructor for the reaction module.
The mapping from grid cells to reaction cells is defined by @ref CreateMapping.
The number of reaction cells may be less than the number of grid cells if
there are inactive regions or there is symmetry in the model definition.
@retval                 Number of grid cells in the user's model.
@see                    @ref PhreeqcRM::PhreeqcRM ,  @ref CreateMapping.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Number of grid cells in the user's model: " << phreeqc_rm.GetGridCellCount() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetGridCellCount(void) {return this->nxyz;}
/**
Returns an IPhreeqc pointer to the @a ith IPhreeqc instance in the reaction module.
For the threaded version, there are @a nthreads + 2 IPhreeqc instances, where
@a nthreads is defined in the constructor (@ref PhreeqcRM::PhreeqcRM).
The number of threads can be determined by @ref GetThreadCount.
The first @a nthreads (0 based) instances will be the workers, the
next (@a nthreads) is the InitialPhreeqc instance, and the next (@a nthreads + 1) is the Utility instance.
Getting the IPhreeqc pointer for one of these instances allows the user to use any of the IPhreeqc methods
on that instance.
For MPI, each process has exactly three IPhreeqc instances, one worker (number 0),
one InitialPhreeqc instance (number 1), and one Utility instance (number 2).
@param i                The number of the IPhreeqc instance (0 based) to be retrieved.
@retval                 IPhreeqc pointer to the @a ith IPhreeqc instance (0 based) in the reaction module.
@see                    @ref PhreeqcRM::PhreeqcRM, @ref GetThreadCount. See IPhreeqc documentation for descriptions of IPhreeqc methods.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// Utility pointer is worker nthreads + 1
IPhreeqc * util_ptr = phreeqc_rm.GetIPhreeqcPointer(phreeqc_rm.GetThreadCount() + 1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
IPhreeqc *                                GetIPhreeqcPointer(int i);

/**
Returns a reference to the vector of all kinetic reactions in the initial-phreeqc module.
The list includes all kinetic reactions included in any KINETICS definitions in
the reaction model.
@ref FindComponents must be called before @ref GetKineticReactions.
This method may be useful when generating selected output definitions related to kinetic reactions.

@retval const std::vector<std::string>&       A vector of strings; each string is a unique
kinetic reaction name.

@see                    @ref FindComponents,
@ref GetKineticReactionsCount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -kinetics " << "\n";
// kinetic reactions
const std::vector<std::string> &kin_reactions = phreeqc_rm.GetKineticReactions();
for (size_t i = 0; i < phreeqc_rm.GetKineticReactionsCount(); i++)
{
oss << "    " << kin_reactions[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetKineticReactions(void) const { return this->KineticReactionsList; }
/**
Returns the number of kinetic reactions in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetKineticReactionsCount.
This method may be useful when generating selected output definitions related to
kinetic reactions.

@retval                 The number of kinetic reactions in the initial-phreeqc module.

@see                    @ref FindComponents,
@ref GetKineticReactions.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -kinetics " << "\n";
// kinetic reactions
const std::vector<std::string> &kin_reactions = phreeqc_rm.GetKineticReactions();
for (size_t i = 0; i < phreeqc_rm.GetKineticReactionsCount(); i++)
{
oss << "    " << kin_reactions[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
int                                       GetKineticReactionsCount(void) const { return (int) this->KineticReactionsList.size(); }

/**
Returns the MPI process (task) number. For the MPI version,
the root task number is zero, and all MPI tasks have unique task numbers greater than zero.
The number of tasks can be obtained with @ref GetMpiTasks. The number of
tasks and computer hosts is determined at run time by the mpiexec command, and the
number of reaction-module processes is defined by the communicator used in
constructing the reaction modules (@ref PhreeqcRM::PhreeqcRM).
For the OPENMP version, the task number is always
zero, and the result of @ref GetMpiTasks is one.
@retval                 The MPI task number for a process.
@see                    @ref GetMpiTasks, @ref PhreeqcRM::PhreeqcRM.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "MPI task number: " << phreeqc_rm.GetMpiMyself() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                 GetMpiMyself(void) const {return this->mpi_myself;}
/**
Returns the number of MPI processes (tasks) assigned to the reaction module.
For the MPI version, the number of
tasks and computer hosts is specified at run time by the mpiexec command. The number of MPI processes
used for reaction calculations is determined by the MPI communicator
used in constructing the reaction modules. The communicator may define a subset of the
total number of MPI processes.
The root task number is zero, and all other MPI tasks have unique task numbers greater than zero.
For the OPENMP version, the number of tasks is
one, and the task number returned by @ref GetMpiMyself is zero.
@retval                 The number of MPI processes assigned to the reaction module.
@see                    @ref GetMpiMyself, @ref PhreeqcRM::PhreeqcRM.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Number of MPI processes: " << phreeqc_rm.GetMpiTasks() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                 GetMpiTasks(void) const {return this->mpi_tasks;}
/**
Returns the user number for the @a nth selected-output definition.
Definitions are sorted by user number. Phreeqc allows multiple selected-output
definitions, each of which is assigned a nonnegative integer identifier by the
user. The number of definitions can be obtained by @ref GetSelectedOutputCount.
To cycle through all of the definitions, GetNthSelectedOutputUserNumber
can be used to identify the user number for each selected-output definition
in sequence. @ref SetCurrentSelectedOutputUserNumber is then used to select
that user number for selected-output processing.
@param n                The sequence number of the selected-output definition for which the user number will be returned.
Fortran, 1 based; C, 0 based.
@retval                 The user number of the @a nth selected-output definition, negative is failure (See @ref DecodeError).
@see                    @ref GetSelectedOutput,
@ref GetSelectedOutputColumnCount, @ref GetSelectedOutputCount,
@ref GetSelectedOutputHeading,
@ref GetSelectedOutputRowCount, @ref SetCurrentSelectedOutputUserNumber, @ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
{
  int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
  status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
  std::cerr << "Selected output sequence number: " << isel << "\n";
  std::cerr << "Selected output user number:     " << n_user << "\n";
  std::vector<double> so;
  int col = phreeqc_rm.GetSelectedOutputColumnCount();
  status = phreeqc_rm.GetSelectedOutput(so);
  // Process results here
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	int                                       GetNthSelectedOutputUserNumber(int n);
/**
Returns the setting for partitioning solids between the saturated and unsaturated
parts of a partially saturated cell. The option is intended to be used by saturated-only
flow codes that allow a variable water table.
The value has meaning only when saturations
less than 1.0 are encountered. The partially saturated cells
may have a small water-to-rock ratio that causes
reactions to proceed slowly relative to fully saturated cells.
By setting  @ref SetPartitionUZSolids to true, the
amounts of solids and gases are partioned according to the saturation.
If a cell has a saturation of 0.5, then
the water interacts with only half of the solids and gases; the other half is unreactive
until the water table rises. As the saturation in a cell varies,
solids and gases are transferred between the
saturated and unsaturated (unreactive) reservoirs of the cell.
Unsaturated-zone flow and transport codes will probably use the default (false),
which assumes all gases and solids are reactive regardless of saturation.
@retval bool       @a True, the fraction of solids and gases available for
reaction is equal to the saturation;
@a False (default), all solids and gases are reactive regardless of saturation.
@see                @ref SetPartitionUZSolids.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "Partioning of UZ solids: " << phreeqc_rm.GetPartitionUZSolids();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	bool                                GetPartitionUZSolids(void) const {return this->partition_uz_solids;}
#ifdef USE_RV
/**
Returns the current set of pore volumes as
defined by the last use of @ref SetPoreVolume or the default (0.1 L).
Pore volume is used with cell volume (@ref SetCellVolume) in calculating porosity.
Pore volumes may change as a function of pressure, in which case they can be updated
with @ref SetPoreVolume.
@retval const std::vector<double>&       A vector reference to the pore volumes.
Size of vector is @a nxyz, the number of grid cells in the user's model (@ref GetGridCellCount).
@see                 @ref GetCellVolume, @ref SetCellVolume, @ref SetPoreVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<double> & vol = phreeqc_rm.GetPoreVolume();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	std::vector<double> &                     GetPoreVolume(void) {return this->pore_volume;}
#endif
/**
Returns the pressure for each cell.
By default, the pressure vector is initialized with 1 atm;
if @ref SetPressure has not been called, worker solutions will have pressures as defined in
input files (@ref RunFile) or input strings (@ref RunString); if @ref SetPressure has been called,
worker solutions will have the pressures as defined by @ref SetPressure.
Pressure effects are considered by three PHREEQC databases: phreeqc.dat, Amm.dat, and pitzer.dat.
@retval const std::vector<double>&       A vector reference to the pressures in each cell, in atm.
Size of vector is @a nxyz, the number of grid cells in the user's model (@ref GetGridCellCount).
@see                 @ref SetPressure, @ref GetTemperature, @ref SetTemperature.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<double> & p_atm = phreeqc_rm.GetPressure();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector<double> &                     GetPressure(void);
/**
Return a reference to the vector of print flags that enable or disable detailed output for each cell.
Printing for a cell will occur only when the
printing is enabled with @ref SetPrintChemistryOn, and the value in the vector for the cell is 1.
@retval std::vector<int> &      Vector of integers. Size of vector is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref GetGridCellCount). A value of 0 for a cell indicates
printing is disabled;
a value of 1 for a cell indicates printing is enabled.
@see                    @ref SetPrintChemistryOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<int> & print_chemistry_mask1 = phreeqc_rm.GetPrintChemistryMask();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	const std::vector<int> &                  GetPrintChemistryMask (void) {return this->print_chem_mask_root;}
/**
Returns a vector reference to the current print flags for detailed output for the three sets of IPhreeqc instances:
the workers, the InitialPhreeqc instance, and the Utility instance. Dimension of the vector is 3.
Printing of detailed output from reaction calculations to the output file
is enabled when the vector value is true, disabled when false.
The detailed output prints all of the output
typical of a PHREEQC reaction calculation, which includes solution descriptions and the compositions of
all other reactants. The output can be several hundred lines per cell, which can lead to a very
large output file (prefix.chem.txt, @ref OpenFiles). For the worker instances,
the output can be limited to a set of cells
(@ref SetPrintChemistryMask) and, in general, the
amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
(applied by using @ref RunFile or @ref RunString).
Printing the detailed output for the workers is generally used only for debugging,
and PhreeqcRM will run faster when printing detailed output for the workers is disabled (@ref SetPrintChemistryOn).
@retval const std::vector<bool> & Print flag for the workers, InitialPhreeqc, and Utility IPhreeqc instances, respectively.
@see                     @ref SetPrintChemistryOn, @ref SetPrintChemistryMask.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<bool> & print_on = phreeqc_rm.GetPrintChemistryOn();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector <bool> &                GetPrintChemistryOn(void) const {return this->print_chemistry_on;}
/**
Get the load-rebalancing method used for parallel processing.
PhreeqcRM attempts to rebalance the load of each thread or
process such that each
thread or process takes the same amount of time to run its part of a @ref RunCells
calculation. Two algorithms are available: one accounts for cells that were not run
because saturation was zero (true), and the other uses the average time
to run all of the cells assigned to a process or thread (false), .
The methods are similar, but preliminary results indicate the default is better in most cases.
@retval bool           @a True indicates individual
cell run times are used in rebalancing (default); @a False, indicates average run times are used in rebalancing.
@see                    @ref GetRebalanceFraction, @ref SetRebalanceByCell, @ref SetRebalanceFraction.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
bool rebalance = phreeqc_rm.GetRebalanceByCell();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	bool                                      GetRebalanceByCell(void) const {return this->rebalance_by_cell;}
/**
Get the fraction used to determine the number of cells to transfer among threads or processes.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a @ref RunCells
calculation. The rebalancing transfers cell calculations among threads or processes to
try to achieve an optimum balance. @ref SetRebalanceFraction
adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
determine the number of cell transfers that actually are made. A value of zero eliminates
load rebalancing. A value less than 1.0 is suggested to avoid possible oscillations,
where too many cells are transferred at one iteration, requiring reverse transfers at the next iteration.
Default is 0.5.
@retval int       Fraction used in rebalance, 0.0 to 1.0.
@see                    @ref GetRebalanceByCell, @ref SetRebalanceByCell, @ref SetRebalanceFraction.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
double f_rebalance = phreeqc_rm.GetRebalanceFraction();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	double                                    GetRebalanceFraction(void) const {return this->rebalance_fraction;}
/**
Returns a vector of saturations (@a sat) as calculated by the reaction module.
Reactions will change the volume of solution in a cell.
The transport code must decide whether to ignore or account for this change in solution volume due to reactions.
Following reactions, the cell saturation is calculated as solution volume (@ref GetSolutionVolume)
divided by the product of representative volume (@ref SetRepresentativeVolume) and the porosity (@ref SetPorosity).
The cell saturation returned by @a GetSaturation may be less than or greater than the saturation set by the transport code
(@ref SetSaturation), and may be greater than or less than 1.0, even in fully saturated simulations.
Only the following databases distributed with PhreeqcRM have molar volume information needed
to accurately calculate solution volume and saturation: phreeqc.dat, Amm.dat, and pitzer.dat.

@param sat              Vector to receive the saturations. Dimension of the array is set to @a nxyz,
where @a nxyz is the number of user grid cells (@ref GetGridCellCount).
Values for inactive cells are set to 1e30.

@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).

@see                    @ref GetSolutionVolume, @ref SetPorosity, @ref SetRepresentativeVolume, @ref SetSaturation.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> sat;
status = phreeqc_rm.GetSaturation(sat);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
IRM_RESULT               GetSaturation(std::vector<double> & sat);
/**
Returns the array of selected-output values for the current selected-output definition.
@ref SetCurrentSelectedOutputUserNumber
specifies which of the selected-output definitions is returned to the vector (@a so).
@param so               A vector to contain the selected-output values.
Size of the vector is set to @a col times @a nxyz, where @a col is the number of
columns in the selected-output definition (@ref GetSelectedOutputColumnCount),
and @a nxyz is the number of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetNthSelectedOutputUserNumber,
@ref GetSelectedOutputColumnCount, @ref GetSelectedOutputCount, @ref GetSelectedOutputHeading,
@ref GetSelectedOutputRowCount, @ref SetCurrentSelectedOutputUserNumber, @ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
{
  int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
  status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
  std::vector<double> so;
  int col = phreeqc_rm.GetSelectedOutputColumnCount();
  status = phreeqc_rm.GetSelectedOutput(so);
  // Process results here
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                GetSelectedOutput(std::vector<double> &so);
/**
Returns the number of columns in the current selected-output definition.
@ref SetCurrentSelectedOutputUserNumber specifies which of the selected-output definitions is used.
@retval                 Number of columns in the current selected-output definition,
negative is failure (See @ref DecodeError).
@see                    @ref GetNthSelectedOutputUserNumber, @ref GetSelectedOutput,
@ref GetSelectedOutputCount, @ref GetSelectedOutputHeading,
@ref GetSelectedOutputRowCount, @ref SetCurrentSelectedOutputUserNumber, @ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
{
  int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
  status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
  std::vector<double> so;
  int col = phreeqc_rm.GetSelectedOutputColumnCount();
  status = phreeqc_rm.GetSelectedOutput(so);
  // Print results
  for (int i = 0; i < phreeqc_rm.GetSelectedOutputRowCount()/2; i++)
  {
    std::vector<std::string> headings;
    headings.resize(col);
    std::cerr << "     Selected output: " << "\n";
    for (int j = 0; j < col; j++)
    {
      status = phreeqc_rm.GetSelectedOutputHeading(j, headings[j]);
      std::cerr << "          " << j << " " << headings[j] << ": " << so[j*nxyz + i] << "\n";
    }
  }
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	int                                       GetSelectedOutputColumnCount(void);
/**
Returns the number of selected-output definitions.
@ref SetCurrentSelectedOutputUserNumber specifies which of the selected-output definitions is used.
@retval                 Number of selected-output definitions, negative is failure (See @ref DecodeError).
@see                    @ref GetNthSelectedOutputUserNumber, @ref GetSelectedOutput,
@ref GetSelectedOutputColumnCount, @ref GetSelectedOutputHeading,
@ref GetSelectedOutputRowCount, @ref SetCurrentSelectedOutputUserNumber, @ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
{
  int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
  status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
  std::vector<double> so;
  int col = phreeqc_rm.GetSelectedOutputColumnCount();
  status = phreeqc_rm.GetSelectedOutput(so);
  // Process results here
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	int                                       GetSelectedOutputCount(void);
/**
Returns a selected-output heading.
The number of headings is determined by @ref GetSelectedOutputColumnCount.
@ref SetCurrentSelectedOutputUserNumber specifies which of the selected-output definitions is used.
@param icol             The sequence number of the heading to be retrieved, 0 based.
@param heading          A string to receive the heading.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetNthSelectedOutputUserNumber, @ref GetSelectedOutput,
@ref GetSelectedOutputColumnCount, @ref GetSelectedOutputCount,
@ref GetSelectedOutputRowCount, @ref SetCurrentSelectedOutputUserNumber, @ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
{
  int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
  status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
  std::vector<double> so;
  int col = phreeqc_rm.GetSelectedOutputColumnCount();
  status = phreeqc_rm.GetSelectedOutput(so);
  std::vector<std::string> headings;
  headings.resize(col);
  std::cerr << "     Selected output: " << "\n";
  for (int j = 0; j < col; j++)
  {
    status = phreeqc_rm.GetSelectedOutputHeading(j, headings[j]);
    std::cerr << "          " << j << " " << headings[j] << "\n";
  }
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                GetSelectedOutputHeading(int icol, std::string &heading);
/**
Returns the current value of the selected-output property.
A value of true for this property indicates that selected output data will be requested this time step.
A value of false indicates that selected output will not be retrieved for this time step;
processing the selected output is avoided with some time savings.
@retval bool      @a True, selected output will be requested; @a false, selected output will not be retrieved.
@see              @ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
bool so_on = phreeqc_rm.GetSelectedOutputOn();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	bool                                      GetSelectedOutputOn(void) {return this->selected_output_on;}
/**
Returns the number of rows in the current selected-output definition. However, the method
is included only for convenience; the number of rows is always equal to the number of
grid cells in the user's model (@ref GetGridCellCount).
@retval                 Number of rows in the current selected-output definition, negative is failure
(See @ref DecodeError).
@see                    @ref GetNthSelectedOutputUserNumber, @ref GetSelectedOutput, @ref GetSelectedOutputColumnCount,
@ref GetSelectedOutputCount, @ref GetSelectedOutputHeading,
@ref SetCurrentSelectedOutputUserNumber, @ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
{
  int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
  status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
  std::vector<double> so;
  int col = phreeqc_rm.GetSelectedOutputColumnCount();
  status = phreeqc_rm.GetSelectedOutput(so);
  // Print results
  for (int i = 0; i < phreeqc_rm.GetSelectedOutputRowCount()/2; i++)
  {
    std::vector<std::string> headings;
    headings.resize(col);
    std::cerr << "     Selected output: " << "\n";
    for (int j = 0; j < col; j++)
    {
      status = phreeqc_rm.GetSelectedOutputHeading(j, headings[j]);
      std::cerr << "          " << j << " " << headings[j] << ": " << so[j*nxyz + i] << "\n";
    }
  }
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
int                                       GetSelectedOutputRowCount(void);

/**
Returns the number of phases in the initial-phreeqc module for which saturation indices could be calculated.
@ref FindComponents must be called before @ref GetSICount.
This method may be useful when generating selected output definitions related to
saturation indices.

@retval                 The number of phases in the initial-phreeqc module for which saturation indices
could be calculated.

@see                    @ref FindComponents,
@ref GetSINames.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<std::string> &si = phreeqc_rm.GetSINames();
for (size_t i = 0; i < phreeqc_rm.GetSICount(); i++)
{
oss << "    " << si[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
int                                       GetSICount(void) const { return (int) this->SINamesList.size(); }
/**
Returns a reference to the vector of the names of all phases for which
saturation indices (SIs) could be calculated.
The list includes all phases that contain only elements included in the components in
the initial-phreeqc module.
The list assumes that all components are present to be able to calculate the entire list of SIs;
it may be that one or more components are missing in any specific cell.
@ref FindComponents must be called before @ref GetSINames.
This method may be useful when generating selected output definitions related to
saturation indices.

@retval const std::vector<std::string>&       A vector of strings; each string is a unique
phase name.

@see                    @ref FindComponents,
@ref GetSICount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -saturation_indices " << "\n";
// molalities of aqueous species
const std::vector<std::string> &si = phreeqc_rm.GetSINames();
for (size_t i = 0; i < phreeqc_rm.GetSICount(); i++)
{
oss << "    " << si[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetSINames(void) const { return this->SINamesList; }

/**
Returns a reference to the vector of solid solution components.
The list of solid solution components includes all components in any SOLID_SOLUTION
definitions in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetSolidSolutionComponents.
This method may be useful when generating selected output definitions related to
solid solutions.

@retval const std::vector<std::string>&       A vector of strings; each string is a
unique solid solution component.

@see                    @ref FindComponents,
@ref GetSolidSolutionComponentsCount, @ref GetSolidSolutionNames.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -solid_solutions " << "\n";
// solid solutions
const std::vector<std::string> &ss_comps = phreeqc_rm.GetSolidSolutionComponents();
const std::vector<std::string> &ss_names = phreeqc_rm.GetSolidSolutionNames();
for (size_t i = 0; i < phreeqc_rm.GetSolidSolutionComponentsCount(); i++)
{

oss << "    ";
oss.width(15);
oss  << std::left << ss_comps[i];
oss << " # " << ss_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetSolidSolutionComponents(void) const { return this->SolidSolutionComponentsList; }
/**
Returns the number of solid solution components in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetSolidSolutionComponentsCount.
This method may be useful when generating selected output definitions related to solid solutions.

@retval                 The number of solid solution components in the initial-phreeqc module.

@see                    @ref FindComponents,
@ref GetSolidSolutionComponents, @ref GetSolidSolutionNames.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -solid_solutions " << "\n";
// solid solutions
const std::vector<std::string> &ss_comps = phreeqc_rm.GetSolidSolutionComponents();
const std::vector<std::string> &ss_names = phreeqc_rm.GetSolidSolutionNames();
for (size_t i = 0; i < phreeqc_rm.GetSolidSolutionComponentsCount(); i++)
{

oss << "    ";
oss.width(15);
oss  << std::left << ss_comps[i];
oss << " # " << ss_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
int                                       GetSolidSolutionComponentsCount(void) const { return (int) this->SolidSolutionComponentsList.size(); }

/**
Returns a reference to the vector of solid solution names that correspond with
the solid solution components.
@ref FindComponents must be called before @ref GetSolidSolutionNames.
The solid solution names vector is the same length as the solid solution components vector
and provides the corresponding name of solid solution containing the component.
This method may be useful when generating selected output definitions related to solid solutions.

@retval const std::vector<std::string>&       A vector of strings; each string is a
solid solution name corresponding to the solid solution components vector; a solid solution name may occur
multiple times.

@see                    @ref FindComponents,
@ref GetSolidSolutionComponentsCount, @ref GetSolidSolutionComponents.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
oss << "  -solid_solutions " << "\n";
// solid solutions
const std::vector<std::string> &ss_comps = phreeqc_rm.GetSolidSolutionComponents();
const std::vector<std::string> &ss_names = phreeqc_rm.GetSolidSolutionNames();
for (size_t i = 0; i < phreeqc_rm.GetSolidSolutionComponentsCount(); i++)
{

oss << "    ";
oss.width(15);
oss  << std::left << ss_comps[i];
oss << " # " << ss_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetSolidSolutionNames(void) const { return this->SolidSolutionNamesList; }

/**
Return a vector reference to the current solution volumes as calculated by the reaction module.
Dimension of the vector will be @a nxyz, where @a nxyz is the number of user grid cells.
Values for inactive cells are set to 1e30.
Only the following databases distributed with PhreeqcRM have molar volume information
needed to accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
@retval Vector reference to current solution volumes.
@see                    @ref GetSaturation.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.RunCells();
const std::vector<double> &volume = phreeqc_rm.GetSolutionVolume();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	const std::vector<double> &               GetSolutionVolume(void);
/**
Returns a vector reference to aqueous species concentrations (@a species_conc).
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
The list of aqueous species is determined by @ref FindComponents and includes all
aqueous species that can be made from the set of components.
Solution volumes used to calculate mol/L are calculated by the reaction module.
Only the following databases distributed with PhreeqcRM have molar volume information
needed to accurately calculate solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.

@param species_conc     Vector to receive the aqueous species concentrations.
Dimension of the vector is set to @a nspecies times @a nxyz,
where @a nspecies is the number of aqueous species (@ref GetSpeciesCount),
and @a nxyz is the number of grid cells (@ref GetGridCellCount).
Concentrations are moles per liter.
Values for inactive cells are set to 1e30.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref FindComponents, 
@ref GetSpeciesCount, 
@ref GetSpeciesD25, 
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module. 

@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
int ncomps = phreeqc_rm.FindComponents();
int npecies = phreeqc_rm.GetSpeciesCount();
status = phreeqc_rm.RunCells();
std::vector<double> c;
status = phreeqc_rm.GetSpeciesConcentrations(c);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                GetSpeciesConcentrations(std::vector<double> & species_conc);
/**
Returns the number of aqueous species used in the reaction module.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
The list of aqueous species is determined by @ref FindComponents and includes all
aqueous species that can be made from the set of components.
@retval int      The number of aqueous species.
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesD25, 
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module. 
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
int ncomps = phreeqc_rm.FindComponents();
int npecies = phreeqc_rm.GetSpeciesCount();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetSpeciesCount(void) {return (int) this->species_names.size();}
/**
Returns a vector reference to diffusion coefficients at 25C for the set of aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
Diffusion coefficients are defined in SOLUTION_SPECIES data blocks, normally in the database file.
Databases distributed with the reaction module that have diffusion coefficients defined are
phreeqc.dat, Amm.dat, and pitzer.dat.
@retval Vector reference to the diffusion coefficients at 25 C, m^2/s. Dimension of the vector is @a nspecies,
where @a nspecies is the number of aqueous species (@ref GetSpeciesCount).
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount,
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module.  
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
int ncomps = phreeqc_rm.FindComponents();
int npecies = phreeqc_rm.GetSpeciesCount();
const std::vector < double > & species_d = phreeqc_rm.GetSpeciesD25();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector<double> &               GetSpeciesD25(void) {return this->species_d_25;}
	/**
	Returns a vector reference to log10 aqueous species activity coefficients (@a species_log10gammas).
	This method is intended for use with multicomponent-diffusion transport calculations,
	and @ref SetSpeciesSaveOn must be set to @a true.
	The list of aqueous species is determined by @ref FindComponents and includes all
	aqueous species that can be made from the set of components.

	@param species_log10gammas     Vector to receive the log10 aqueous species activity coefficients.
	Dimension of the vector is set to @a nspecies times @a nxyz,
	where @a nspecies is the number of aqueous species (@ref GetSpeciesCount),
	and @a nxyz is the number of grid cells (@ref GetGridCellCount).
	Values for inactive cells are set to 1e30.
	@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
	@see                    @ref FindComponents
    @ref GetSpeciesConcentrations, 
	@ref GetSpeciesCount,
	@ref GetSpeciesD25,
	@ref GetSpeciesNames,
	@ref GetSpeciesSaveOn,
	@ref GetSpeciesStoichiometry,
	@ref GetSpeciesZ,
	@ref SetSpeciesSaveOn,
	@ref SpeciesConcentrations2Module.

	@par C++ Example:
	@htmlonly
	<CODE>
	<PRE>
	status = phreeqc_rm.SetSpeciesSaveOn(true);
	int ncomps = phreeqc_rm.FindComponents();
	int npecies = phreeqc_rm.GetSpeciesCount();
	status = phreeqc_rm.RunCells();
	std::vector<double> species_gammas;
	status = phreeqc_rm.GetSpeciesLog10Gammas(species_gammas);
	</PRE>
	</CODE>
	@endhtmlonly
	@par MPI:
	Called by root, workers must be in the loop of @ref MpiWorker.
	*/
	IRM_RESULT                                GetSpeciesLog10Gammas(std::vector<double> & species_log10gammas);
	
	/**
Returns a vector reference to the names of the aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
The list of aqueous species is determined by @ref FindComponents and includes all
aqueous species that can be made from the set of components.
@retval names      Vector of strings containing the names of the aqueous species. Dimension of the vector is @a nspecies,
where @a nspecies is the number of aqueous species (@ref GetSpeciesCount).
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount,
@ref GetSpeciesD25, 
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module.  
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
int ncomps = phreeqc_rm.FindComponents();
int npecies = phreeqc_rm.GetSpeciesCount();
const std::vector<std::string> &species = phreeqc_rm.GetSpeciesNames();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector<std::string> &          GetSpeciesNames(void) {return this->species_names;}
/**
Returns the value of the species-save property.
By default, concentrations of aqueous species are not saved. Setting the species-save property to true allows
aqueous species concentrations to be retrieved
with @ref GetSpeciesConcentrations, and solution compositions to be set with
@ref SpeciesConcentrations2Module.

@retval True indicates solution species concentrations are saved and can be used for multicomponent-diffusion calculations;
@a False indicates that solution species concentrations are not saved.
@see                    @ref FindComponents, @ref GetSpeciesConcentrations, @ref GetSpeciesCount,
@ref GetSpeciesD25, @ref GetSpeciesSaveOn, @ref GetSpeciesZ,
@ref GetSpeciesNames, @ref SpeciesConcentrations2Module.
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount,
@ref GetSpeciesD25,
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
int ncomps = phreeqc_rm.FindComponents();
int npecies = phreeqc_rm.GetSpeciesCount();
bool species_on = phreeqc_rm.GetSpeciesSaveOn();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	bool                                      GetSpeciesSaveOn(void) {return this->species_save_on;}

/**
Returns a vector reference to the stoichiometry of each aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.

@retval Vector of cxxNameDouble instances (maps) that contain the component names and
associated stoichiometric coefficients for each aqueous species.  Dimension of the vector is @a nspecies,
where @a nspecies is the number of aqueous species (@ref GetSpeciesCount).
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount,
@ref GetSpeciesD25,
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesZ,
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<std::string> &species = phreeqc_rm.GetSpeciesNames();
const std::vector < double > & species_z = phreeqc_rm.GetSpeciesZ();
const std::vector < double > & species_d = phreeqc_rm.GetSpeciesD25();
bool species_on = phreeqc_rm.GetSpeciesSaveOn();
int nspecies = phreeqc_rm.GetSpeciesCount();
for (int i = 0; i < nspecies; i++)
{
  std::ostringstream strm;
  strm << species[i] << "\n";
  strm << "    Charge: " << species_z[i] << std::endl;
  strm << "    Dw:     " << species_d[i] << std::endl;
  cxxNameDouble::const_iterator it = phreeqc_rm.GetSpeciesStoichiometry()[i].begin();
  for (; it != phreeqc_rm.GetSpeciesStoichiometry()[i].end(); it++)
  {
    strm << "          " << it->first << "   " << it->second << "\n";
  }
  phreeqc_rm.OutputMessage(strm.str());
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */

	const std::vector<cxxNameDouble> &        GetSpeciesStoichiometry(void) {return this->species_stoichiometry;}
/**
Returns a vector reference to the charge on each aqueous species.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
@retval Vector containing the charge on each aqueous species. Dimension of the vector is @a nspecies,
where @a nspecies is the number of aqueous species (@ref GetSpeciesCount).
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount,
@ref GetSpeciesD25,
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref SetSpeciesSaveOn,
@ref SpeciesConcentrations2Module.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
int ncomps = phreeqc_rm.FindComponents();
int npecies = phreeqc_rm.GetSpeciesCount();
const std::vector < double > & species_z = phreeqc_rm.GetSpeciesZ();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector<double> &               GetSpeciesZ(void) {return this->species_z;}
/**
Returns a vector of integers that contains the smallest reaction-cell number assigned to each worker.
Each worker is assigned a range of reaction-cell numbers that are run during a call to @ref RunCells.
The range of reaction cell numbers for a worker may vary as load rebalancing occurs.
At any point in the calculations, the first cell and last cell to be run by a worker can be found
in the vectors returned by @a GetStartCell and @ref GetEndCell.
Each method returns a vector of integers that has size of the number of threads (@ref GetThreadCount),
if using OPENMP, or the number of processes (@ref GetMpiTasks), if using MPI.
@retval IRM_RESULT      Vector of integers, one for each worker, that gives the first reaction cell
to be run by each worker.
@see                    @ref GetEndCell, @ref GetThreadCount, @ref GetMpiTasks, @ref RunCells.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Current distribution of cells for workers\n";
oss << "Worker First Cell   Last Cell\n";
int n;
n = phreeqc_rm.GetThreadCount() * phreeqc_rm.GetMpiTasks();
for (int i = 0; i < n; i++)
{
	oss << i << "      "
	    << phreeqc_rm.GetStartCell()[i]
	    << "            "
		<< phreeqc_rm.GetEndCell()[i] << "\n";
}
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
const std::vector < int> &                GetStartCell(void) const {return this->start_cell;}

/**
Returns a reference to the vector of surface names (such as "Hfo") that correspond with
the surface species names. The vectors referenced by @ref GetSurfaceSpecies
and @ref GetSurfaceNames are the same length.
@ref FindComponents must be called before @ref GetSurfaceNames.
This method may be useful when generating selected output definitions related to surfaces.

@retval const std::vector<std::string>&       A vector of strings; each string is a
surface name corresponding to the surface species vector;
a surface name may occur multiple times.

@see                    @ref FindComponents,
@ref GetSurfaceSpeciesCount, @ref GetSurfaceSpecies, @ref GetSurfaceTypes.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// molalities of surface species
const std::vector<std::string> &surf_species = phreeqc_rm.GetSurfaceSpecies();
const std::vector<std::string> &surf_types = phreeqc_rm.GetSurfaceTypes();
const std::vector<std::string> &surf_names = phreeqc_rm.GetSurfaceNames();
for (size_t i = 0; i < phreeqc_rm.GetSurfaceSpeciesCount(); i++)
{
oss << "    ";
oss.width(15);
oss << std::left << surf_species[i];
oss << " # ";
oss.width(15);
oss << surf_types[i] << "   " << surf_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetSurfaceNames(void) const { return this->SurfaceNamesList; }

/**
Returns a reference to the vector of surface species names (such as "Hfo_wOH").
The list of surface species is derived from the list of components
(@ref FindComponents) and the list of all surface site types (such as "Hfo_w")
that are included in SURFACE definitions in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetSurfaceSpecies.
This method may be useful when generating selected output definitions related to surfaces.

@retval const std::vector<std::string>&       A vector of strings; each string is a
unique surface species name.
@see                    @ref FindComponents,
@ref GetSurfaceSpeciesCount, @ref GetSurfaceTypes, @ref GetSurfaceNames.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// molalities of surface species
const std::vector<std::string> &surf_species = phreeqc_rm.GetSurfaceSpecies();
const std::vector<std::string> &surf_types = phreeqc_rm.GetSurfaceTypes();
const std::vector<std::string> &surf_names = phreeqc_rm.GetSurfaceNames();
for (size_t i = 0; i < phreeqc_rm.GetSurfaceSpeciesCount(); i++)
{
oss << "    ";
oss.width(15);
oss << std::left << surf_species[i];
oss << " # ";
oss.width(15);
oss << surf_types[i] << "   " << surf_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetSurfaceSpecies(void) const { return this->SurfaceSpeciesNamesList; }

/**
Returns the number of surface species (such as "Hfo_wOH") in the initial-phreeqc module.
@ref FindComponents must be called before @ref GetSurfaceSpeciesCount.
This method may be useful when generating selected output definitions related to surfaces.

@retval                 The number of surface species in the initial-phreeqc module.

@see                    @ref FindComponents,
@ref GetSurfaceSpecies, @ref GetSurfaceTypes, @ref GetSurfaceNames.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// molalities of surface species
const std::vector<std::string> &surf_species = phreeqc_rm.GetSurfaceSpecies();
const std::vector<std::string> &surf_types = phreeqc_rm.GetSurfaceTypes();
const std::vector<std::string> &surf_names = phreeqc_rm.GetSurfaceNames();
for (size_t i = 0; i < phreeqc_rm.GetSurfaceSpeciesCount(); i++)
{
oss << "    ";
oss.width(15);
oss << std::left << surf_species[i];
oss << " # ";
oss.width(15);
oss << surf_types[i] << "   " << surf_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
int                                       GetSurfaceSpeciesCount(void) const { return (int) this->SurfaceSpeciesNamesList.size(); }


/**
Returns a reference to the vector of surface site types (such as "Hfo_w") that correspond with
the surface species names.
The vectors referenced by @ref GetSurfaceSpecies and
@ref GetSurfaceTypes are the same length.
@ref FindComponents must be called before @ref GetSurfaceTypes.
This method may be useful when generating selected output definitions related to surfaces.

@retval const std::vector<std::string>&       A vector of strings; each string is a
surface site type for the corresponding species in the surface species vector;
a surface site type may occur multiple times.

@see                    @ref FindComponents,
@ref GetSurfaceSpeciesCount, @ref GetSurfaceSpecies, @ref GetSurfaceNames.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
// molalities of surface species
const std::vector<std::string> &surf_species = phreeqc_rm.GetSurfaceSpecies();
const std::vector<std::string> &surf_types = phreeqc_rm.GetSurfaceTypes();
const std::vector<std::string> &surf_names = phreeqc_rm.GetSurfaceNames();
for (size_t i = 0; i < phreeqc_rm.GetSurfaceSpeciesCount(); i++)
{
oss << "    ";
oss.width(15);
oss << std::left << surf_species[i];
oss << " # ";
oss.width(15);
oss << surf_types[i] << "   " << surf_names[i] << "\n";
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
*/
const std::vector<std::string> &          GetSurfaceTypes(void) const { return this->SurfaceTypesList; }

/**
Vector reference to the current temperatures of the cells.
By default, the temperature vector is initialized to 25 C;
if @ref SetTemperature has not been called, worker solutions will have temperatures as defined in
input files (@ref RunFile) or input strings (@ref RunString); if @ref SetTemperature has been called,
worker solutions will have the temperatures as defined by @ref SetTemperature.
@retval Vector of temperatures, in degrees C. Size of vector is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref GetGridCellCount).
@see                    @ref SetTemperature, @ref GetPressure, @ref SetPressure.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector<double> &  tempc = phreeqc_rm.GetTemperature();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	//const std::vector<double> &               GetTemperature(void) {return this->tempc;}
	const std::vector<double> &               GetTemperature(void);
/**
Returns the number of threads, which is equal to the number of workers used to run in parallel with OPENMP.
For the OPENMP version, the number of threads is set implicitly or explicitly
with the constructor (@ref PhreeqcRM::PhreeqcRM).
For the MPI version, the number of threads is always one for each process.
@retval                 The number of threads used for OPENMP parallel processing.
@see                    @ref GetMpiTasks.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Number of threads: " << phreeqc_rm.GetThreadCount() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers; result is always 1.
 */
	int                                       GetThreadCount() {return this->nthreads;}
/**
Returns the current simulation time in seconds.
The reaction module does not change the time value, so the
returned value is equal to the default (0.0) or the last time set by @ref SetTime.
@retval                 The current simulation time, in seconds.
@see                    @ref GetTimeConversion, @ref GetTimeStep, @ref SetTime,
@ref SetTimeConversion, @ref SetTimeStep.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream strm;
strm << "Beginning transport calculation "
     << phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion()
	 << " days\n";
phreeqc_rm.LogMessage(strm.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	double                                    GetTime(void) const {return this->time;}
/**
Returns a multiplier to convert time from seconds to another unit, as specified by the user.
The reaction module uses seconds as the time unit. The user can set a conversion
factor (@ref SetTimeConversion) and retrieve it with GetTimeConversion.
The reaction module only uses the conversion factor when printing the long version
of cell chemistry (@ref SetPrintChemistryOn), which is rare.
Default conversion factor is 1.0.
@retval                 Multiplier to convert seconds to another time unit.
@see                    @ref GetTime, @ref GetTimeStep, @ref SetTime, @ref SetTimeConversion, @ref SetTimeStep.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream strm;
strm << "Beginning transport calculation "
     <<   phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion()
	 << " days\n";
phreeqc_rm.LogMessage(strm.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	double                                    GetTimeConversion(void) {return this->time_conversion;}
/**
Returns the current simulation time step in seconds.
This is the time over which kinetic reactions are integrated in a call to @ref RunCells.
The reaction module does not change the time-step value, so the
returned value is equal to the default (0.0) or the last time step set by @ref SetTimeStep.
@retval                 The current simulation time step, in seconds.
@see                    @ref GetTime, @ref GetTimeConversion, @ref SetTime,
@ref SetTimeConversion, @ref SetTimeStep.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream strm;
strm << "Time step "
     << phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion()
	 << " days\n";
phreeqc_rm.LogMessage(strm.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	double                                    GetTimeStep(void) {return this->time_step;}
/**
Returns the input units for exchangers.
In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
@ref SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

@retval                 Input units for exchangers.
@see                    @ref SetUnitsExchange, @ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int units_exchange = phreeqc_rm.GetUnitsExchange();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetUnitsExchange(void) {return this->units_Exchange;}
/**
Returns the input units for gas phases.
In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
@ref SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

@retval                 Input units for gas phases (0, 1, or 2).
@see                    @ref SetUnitsGasPhase, @ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int units_gas_phase = phreeqc_rm.GetUnitsGasPhase();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetUnitsGasPhase(void) {return this->units_GasPhase;}
/**
Returns the input units for kinetic reactants.
In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
@ref SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

@retval                 Input units for kinetic reactants (0, 1, or 2).
@see                    @ref SetUnitsKinetics, @ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int units_kinetics = phreeqc_rm.GetUnitsKinetics();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetUnitsKinetics(void) {return this->units_Kinetics;}
/**
Returns the input units for pure phase assemblages (equilibrium phases).
In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
@ref SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

@retval                 Input units for equilibrium phases (0, 1, or 2).
@see                    @ref SetUnitsPPassemblage, @ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int units_pp_assemblage = phreeqc_rm.GetUnitsPPassemblage();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetUnitsPPassemblage(void) {return this->units_PPassemblage;}
/**
Returns the units of concentration used by the transport model.
Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
In PHREEQC, solutions are defined by the number of moles of each
element in the solution. The units of transport concentration are used when
transport concentrations are converted to
solution moles by @ref SetConcentrations and @ref Concentrations2Utility.
The units of solution concentration also are used when solution moles are converted to
transport concentrations by
@ref GetConcentrations.
@n@n
To convert from mg/L to moles
of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
multiplied by the solution volume,
which is the product of porosity (@ref SetPorosity), saturation (@ref SetSaturation), and
representative volume (@ref SetRepresentativeVolume).
To convert from mol/L to moles
of element in a cell, mol/L is
multiplied by the solution volume.
To convert from mass fraction to moles
of element in a cell, kg/kgs is converted to mol/kgs, multiplied by density
(@ref SetDensity) and
multiplied by the solution volume.
@n@n
To convert from moles
of element in the representative volume of a reaction cell to mg/L, the number of moles of an element is divided by the
solution volume resulting in mol/L, and then converted to
mg/L.
To convert from moles
of element in the representative volume of a reaction cell to mol/L,  the number of moles of an element is divided by the
solution volume resulting in mol/L.
To convert from moles
of element in the representative volume of a reaction cell to mass fraction,
the number of moles of an element is converted to kg and divided by the total mass of the solution.
Two options are available for the volume and mass of solution
that are used in converting to transport concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or (2) the volume of solution is the product of porosity, saturation, and representative volume,
and the mass of solution is volume times density as defined by @ref SetDensity.
Which option is used is determined by @ref UseSolutionDensityVolume.
@retval                 Units for concentrations in transport.
@see                    @ref Concentrations2Utility, @ref GetConcentrations, @ref SetConcentrations, @ref SetDensity,
@ref SetPorosity, @ref SetRepresentativeVolume, @ref SetUnitsSolution, @ref UseSolutionDensityVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int units_solution = phreeqc_rm.GetUnitsSolution();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetUnitsSolution(void) {return this->units_Solution;}
/**
Returns the input units for solid-solution assemblages.
In PHREEQC input, solid solutions are defined by moles of each component (@a Mp).
@ref SetUnitsSSassemblage specifies how the number of moles in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

@retval                 Input units for solid solutions (0, 1, or 2).
@see                    @ref SetUnitsSSassemblage, @ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int units_ss_exchange = phreeqc_rm.GetUnitsSSassemblage();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetUnitsSSassemblage(void) {return this->units_SSassemblage;}
/**
Returns the input units for surfaces.
In PHREEQC input, surfaces are defined by moles of surface sites  (@a Mp).
@ref SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

@retval                 Input units for solid surfaces (0, 1, or 2).
@see                    @ref SetUnitsSurface, @ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int units_surface = phreeqc_rm.GetUnitsSurface();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	int                                       GetUnitsSurface(void) {return this->units_Surface;}
/**
Returns a reference to the vector of IPhreeqcPhast instances. IPhreeqcPhast
inherits from IPhreeqc, and the vector can be interpreted as a vector of pointers to the worker,
InitialPhreeqc, and Utility IPhreeqc instances. For OPENMP, there are @a nthreads
workers, where @a nthreads is defined in the constructor (@ref PhreeqcRM::PhreeqcRM).
For MPI, there is a single worker. For OPENMP and MPI, there is one InitialPhreeqc and
one Utility instance.
@retval                 Vector of IPhreeqcPhast instances.
@see                    @ref PhreeqcRM::PhreeqcRM.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
const std::vector < IPhreeqcPhast *> & w = phreeqc_rm.GetWorkers();
w[0]->AccumulateLine("Delete; -all");
int iphreeqc_result = w[0]->RunAccumulated();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	const std::vector<IPhreeqcPhast *> &      GetWorkers() {return this->workers;}
/**
Fills a vector (@a destination_c) with concentrations from solutions in the InitialPhreeqc instance.
The method is used to obtain concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell.
@param destination_c       Vector to receive the concentrations.The dimension of @a destination_c is set to @a ncomps times @a n_boundary,
where @a ncomps is the number of components returned from @ref FindComponents or @ref GetComponentCount, and @a n_boundary
is the dimension of the vector @a boundary_solution1.
@param boundary_solution1  Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
Size is @a n_boundary.
@retval IRM_RESULT         0 is success, negative is failure (See @ref DecodeError).
@see                       @ref FindComponents, @ref GetComponentCount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> bc_conc;
std::vector<int> bc1;
int nbound = 1;
bc1.resize(nbound, 0);                      // solution 0 from InitialIPhreeqc instance
status = phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc, bc1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                InitialPhreeqc2Concentrations(
													std::vector < double > & destination_c,
													std::vector < int >    & boundary_solution1);
/**
Fills a vector (@a destination_c) with concentrations from solutions in the InitialPhreeqc instance.
The method is used to obtain concentrations for boundary conditions that are mixtures of solutions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell. Concentrations may be a mixture of two
solutions, @a boundary_solution1 and @a boundary_solution2, with a mixing fraction for @a boundary_solution1 of
@a fraction1 and mixing fraction for @a boundary_solution2 of (1 - @a fraction1).
A negative value for @a boundary_solution2 implies no mixing, and the associated value for @a fraction1 is ignored.
@param destination_c       Vector of concentrations extracted from the InitialPhreeqc instance.
The dimension of @a destination_c is set to @a ncomps times @a n_boundary,
where @a ncomps is the number of components returned from @ref FindComponents or @ref GetComponentCount, and @a n_boundary
is the dimension of the vectors @a boundary_solution1, @a boundary_solution2, and @a fraction1.
@param boundary_solution1  Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
Size is @a n_boundary.
@param boundary_solution2  Vector of solution index numbers that that refer to solutions in the InitialPhreeqc instance
and are defined to mix with @a boundary_solution1.
Size is @a n_boundary.
@param fraction1           Fraction of boundary_solution1 that mixes with (1 - @a fraction1) of @a boundary_solution2.
Size is @a n_boundary.
@retval IRM_RESULT         0 is success, negative is failure (See @ref DecodeError).
@see                  @ref FindComponents, @ref GetComponentCount.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> bc_conc, bc_f1;
std::vector<int> bc1, bc2;
int nbound = 1;
bc1.resize(nbound, 0);                      // solution 0 from InitialIPhreeqc instance
bc2.resize(nbound, 1);                      // solution 1 from InitialIPhreeqc instance
bc_f1.resize(nbound, 0.4);                  // mixing fraction for bc1, result is 0.4/0.6 mix
status = phreeqc_rm.InitialPhreeqc2Concentrations(bc_conc, bc1, bc2, bc_f1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT								  InitialPhreeqc2Concentrations(
													std::vector < double > & destination_c,
													std::vector < int >    & boundary_solution1,
													std::vector < int >    & boundary_solution2,
													std::vector < double > & fraction1);
/**
Transfer solutions and reactants from the InitialPhreeqc instance to the reaction-module workers.
@a Initial_conditions1 is used to select initial conditions, including solutions and reactants,
for each cell of the model, without mixing.
@a Initial_conditions1 is dimensioned 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model
(@ref GetGridCellCount). The dimension of 7 refers to solutions and reactants in the following order:
(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
(5) SOLID_SOLUTIONS, and (6) KINETICS.
The definition initial_solution1[3*nxyz + 99] = 2, indicates that
cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance
(created in the InitialPhreeqc instance either by @ref RunFile or @ref RunString).
@param initial_conditions1 Vector of solution and reactant index numbers that refer to
definitions in the InitialPhreeqc instance.
Size is 7 times @a nxyz. The order of definitions is given above.
Negative values are ignored, resulting in no definition of that entity for that cell.
@retval IRM_RESULT          0 is success, negative is failure (See @ref DecodeError).
@see                        @ref InitialPhreeqcCell2Module.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<int> ic1;
ic1.resize(nxyz*7, -1);
for (int i = 0; i < nxyz; i++)
{
  ic1[i] = 1;              // Solution 1
  ic1[nxyz + i] = -1;      // Equilibrium phases none
  ic1[2*nxyz + i] = 1;     // Exchange 1
  ic1[3*nxyz + i] = -1;    // Surface none
  ic1[4*nxyz + i] = -1;    // Gas phase none
  ic1[5*nxyz + i] = -1;    // Solid solutions none
  ic1[6*nxyz + i] = -1;    // Kinetics none
}
status = phreeqc_rm.InitialPhreeqc2Module(ic1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT  InitialPhreeqc2Module(std::vector < int >    & initial_conditions1);
/**
Transfer solutions and reactants from the InitialPhreeqc instance to the reaction-module workers, possibly with mixing.
In its simplest form, @a  initial_conditions1 is used to select initial conditions, including solutions and reactants,
for each cell of the model, without mixing.
@a Initial_conditions1 is dimensioned 7 times @a  nxyz, where @a  nxyz is the number of grid cells in the user's model
(@ref GetGridCellCount). The dimension of 7 refers to solutions and reactants in the following order:
(0) SOLUTIONS, (1) EQUILIBRIUM_PHASES, (2) EXCHANGE, (3) SURFACE, (4) GAS_PHASE,
(5) SOLID_SOLUTIONS, and (6) KINETICS.
The definition initial_solution1[3*nxyz + 99] = 2, indicates that
cell 99 (0 based) contains the SURFACE definition (index 3) defined by SURFACE 2 in the InitialPhreeqc instance
(either by @ref RunFile or @ref RunString).
@n@n
It is also possible to mix solutions and reactants to obtain the initial conditions for cells. For mixing,
@a initials_conditions2 contains numbers for a second entity that mixes with the entity defined in @a initial_conditions1.
@a Fraction1 contains the mixing fraction for @a initial_conditions1,
whereas (1 - @a fraction1) is the mixing fraction for @a initial_conditions2.
The definitions initial_solution1[3*nxyz + 99] = 2, initial_solution2[3*nxyz + 99] = 3,
fraction1[3*nxyz + 99] = 0.25 indicates that
cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3,
where the surface compositions have been defined in the InitialPhreeqc instance.
If the user number in @a initial_conditions2 is negative, no mixing occurs.
@param initial_conditions1 Vector of solution and reactant index numbers that refer to
definitions in the InitialPhreeqc instance.
Size is 7 times @a nxyz, where @a nxyz is the number of grid cells in the user's model (@ref GetGridCellCount).
The order of definitions is given above.
Negative values are ignored, resulting in no definition of that entity for that cell.
@param initial_conditions2  Vector of solution and reactant index numbers that refer to
definitions in the InitialPhreeqc instance.
Nonnegative values of @a initial_conditions2 result in mixing with the entities defined in @a initial_conditions1.
Negative values result in no mixing.
Size is 7 times @a nxyz. The order of definitions is given above.
@param fraction1           Fraction of @a initial_conditions1 that mixes with (1 - @a fraction1)
of @a initial_conditions2.
Size is 7 times @a nxyz. The order of definitions is given above.
@retval IRM_RESULT          0 is success, negative is failure (See @ref DecodeError).
@see                        @ref InitialPhreeqcCell2Module.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<int> ic1, ic2;
ic1.resize(nxyz*7, -1);
ic2.resize(nxyz*7, -1);
std::vector<double> f1;
f1.resize(nxyz*7, 1.0);
for (int i = 0; i < nxyz; i++)
{
  ic1[i] = 1;              // Solution 1
  ic1[nxyz + i] = -1;      // Equilibrium phases none
  ic1[2*nxyz + i] = 1;     // Exchange 1
  ic1[3*nxyz + i] = -1;    // Surface none
  ic1[4*nxyz + i] = -1;    // Gas phase none
  ic1[5*nxyz + i] = -1;    // Solid solutions none
  ic1[6*nxyz + i] = -1;    // Kinetics none
}
status = phreeqc_rm.InitialPhreeqc2Module(ic1, ic2, f1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT InitialPhreeqc2Module(
		std::vector < int >    & initial_conditions1,
		std::vector < int >    & initial_conditions2,
		std::vector < double > & fraction1);
/**
Fills a vector @a destination_c with aqueous species concentrations from solutions in the InitialPhreeqc instance.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
The method is used to obtain aqueous species concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell.
@param destination_c           Vector of aqueous concentrations extracted from the InitialPhreeqc instance.
The dimension of @a species_c is @a nspecies times @a n_boundary,
where @a nspecies is the number of aqueous species returned from @ref GetSpeciesCount,
and @a n_boundary is the dimension of @a boundary_solution1.
@param boundary_solution1  Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
@retval IRM_RESULT         0 is success, negative is failure (See @ref DecodeError).
@see                  @ref FindComponents, @ref GetSpeciesCount, @ref SetSpeciesSaveOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> bc_conc, bc_f1;
std::vector<int> bc1, bc2;
int nbound = 1;
bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance
status = phreeqc_rm.InitialPhreeqc2SpeciesConcentrations(bc_conc, bc1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                InitialPhreeqc2SpeciesConcentrations(
													std::vector < double > & destination_c,
													std::vector < int >    & boundary_solution1);
/**
Fills a vector @a destination_c with aqueous species concentrations from solutions in the InitialPhreeqc instance.
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
The method is used to obtain aqueous species concentrations for boundary conditions. If a negative value
is used for a cell in @a boundary_solution1, then the highest numbered solution in the InitialPhreeqc instance
will be used for that cell.
Concentrations may be a mixture of two
solutions, @a boundary_solution1 and @a boundary_solution2, with a mixing fraction for @a boundary_solution1 of
@a fraction1 and mixing fraction for @a boundary_solution2 of (1 - @a fraction1).
A negative value for @a boundary_solution2 implies no mixing, and the associated value for @a fraction1 is ignored.
@param destination_c           Vector of aqueous concentrations extracted from the InitialPhreeqc instance.
The dimension of @a species_c is @a nspecies times @a n_boundary,
where @a nspecies is the number of aqueous species returned from @ref GetSpeciesCount, and @a n_boundary is the dimension
of @a boundary_solution1.
@param boundary_solution1  Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance.
@param boundary_solution2  Vector of solution index numbers that refer to solutions in the InitialPhreeqc instance
and are defined to mix with @a boundary_solution1. Size is same as @a boundary_solution1.
@param fraction1           Vector of fractions of @a boundary_solution1 that mix with (1 - @a fraction1) of @a boundary_solution2.
Size is same as @a boundary_solution1.
@retval IRM_RESULT         0 is success, negative is failure (See @ref DecodeError).
@see                  @ref FindComponents, @ref GetSpeciesCount, @ref SetSpeciesSaveOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> bc_conc, bc_f1;
std::vector<int> bc1, bc2;
int nbound = 1;
bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance
bc2.resize(nbound, -1);                     // no bc2 solution for mixing
bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
status = phreeqc_rm.InitialPhreeqc2SpeciesConcentrations(bc_conc, bc1, bc2, bc_f1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT								  InitialPhreeqc2SpeciesConcentrations(
													std::vector < double > & destination_c,
													std::vector < int >    & boundary_solution1,
													std::vector < int >    & boundary_solution2,
													std::vector < double > & fraction1);
/**
A cell numbered @a n in the InitialPhreeqc instance is selected to populate a series of transport cells.
All reactants with the number @a n are transferred along with the solution.
If MIX @a n exists, it is used for the definition of the solution.
If @a n is negative, @a n is redefined to be the largest solution or MIX number in the InitialPhreeqc instance.
All reactants for each cell in the list @a cell_numbers are removed before the cell
definition is copied from the InitialPhreeqc instance to the workers.
@param n                  Number that refers to a solution or MIX and associated reactants in the InitialPhreeqc instance.
@param cell_numbers       A vector of grid-cell numbers (user's grid-cell numbering system) that
will be populated with cell @a n from the InitialPhreeqc instance.
@retval IRM_RESULT        0 is success, negative is failure (See @ref DecodeError).
@see                      @ref InitialPhreeqc2Module.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<int> module_cells;
module_cells.push_back(18);
module_cells.push_back(19);
status = phreeqc_rm.InitialPhreeqcCell2Module(-1, module_cells);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                InitialPhreeqcCell2Module(int n, const std::vector<int> &cell_numbers);
/**
Load a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. All definitions
of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
@param database         String containing the database name.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.LoadDatabase("phreeqc.dat");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                LoadDatabase(const std::string &database);
/**
Print a message to the log file.
@param str              String to be printed.
@see                    @ref OpenFiles, @ref ErrorMessage, @ref OutputMessage, @ref ScreenMessage, @ref WarningMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream strm;
strm << "Beginning transport calculation " <<   phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion() << " days\n";
strm << "          Time step             " <<   phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion() << " days\n";
phreeqc_rm.LogMessage(strm.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	void                                      LogMessage(const std::string &str);
/**
MPI only. Calls MPI_Abort, which aborts MPI, and makes the reaction module unusable.
Should be used only on encountering an unrecoverable error.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
int status = phreeqc_rm.MPI_Abort();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root or workers.
 */
	int                                       MpiAbort();
/**
MPI only. Nonroot processes (processes with @ref GetMpiMyself > 0) must call MpiWorker to be able to
respond to messages from the root to accept data, perform calculations, and
(or) return data within the reaction module.
MpiWorker contains a loop that reads a message from root, performs a
task, and waits for another message from root.
@ref SetConcentrations, @ref RunCells, and @ref GetConcentrations
are examples of methods that send a message from root to get the workers to perform a task.
The workers will respond to all methods that are designated "workers must be in the loop of MpiWorker"
in the MPI section of the method documentation.
The workers will continue to respond to messages from root until root calls
@ref MpiWorkerBreak.
@n@n
(Advanced) The list of tasks that the workers perform can be extended by using @ref SetMpiWorkerCallbackC.
It is then possible to use the MPI processes to perform other developer-defined tasks,
such as transport calculations, without exiting from the MpiWorker loop.
Alternatively, root calls @ref MpiWorkerBreak to allow the workers to continue past a call to MpiWorker.
The workers perform developer-defined calculations, and then MpiWorker is called again to respond to
requests from root to perform reaction-module tasks.

@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
MpiWorker returns a value only when @ref MpiWorkerBreak is called by root.
@see                    @ref MpiWorkerBreak, @ref SetMpiWorkerCallbackC, @ref SetMpiWorkerCallbackCookie.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
PhreeqcRM phreeqc_rm(nxyz, MPI_COMM_WORLD);
int mpi_myself;
if (MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself) != MPI_SUCCESS)
{
  exit(4);
}
if (mpi_myself > 0)
{
  phreeqc_rm.MpiWorker();
  return EXIT_SUCCESS;
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by all workers.
 */
	IRM_RESULT                                MpiWorker();
/**
MPI only. This method is called by root to force nonroot processes (processes with @ref GetMpiMyself > 0)
to return from a call to @ref MpiWorker.
@ref MpiWorker contains a loop that reads a message from root, performs a
task, and waits for another message from root. The workers respond to all methods that are designated
"workers must be in the loop of MpiWorker" in the
MPI section of the method documentation.
The workers will continue to respond to messages from root until root calls MpiWorkerBreak.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref MpiWorker, @ref SetMpiWorkerCallbackC, @ref SetMpiWorkerCallbackCookie (C only).
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.MpiWorkerBreak();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                MpiWorkerBreak();
/**
Opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt
based on the prefix defined by @ref SetFilePrefix.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetFilePrefix, @ref GetFilePrefix, @ref CloseFiles,
@ref ErrorMessage, @ref LogMessage, @ref OutputMessage, @ref WarningMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetFilePrefix("Advect_cpp");
phreeqc_rm.OpenFiles();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                OpenFiles(void);
/**
Print a message to the output file.
@param str              String to be printed.
@see                    @ref ErrorMessage, @ref LogMessage, @ref ScreenMessage, @ref WarningMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream oss;
oss << "Database:                                         " << phreeqc_rm.GetDatabaseFileName().c_str() << "\n";
oss << "Number of threads:                                " << phreeqc_rm.GetThreadCount() << "\n";
oss << "Number of MPI processes:                          " << phreeqc_rm.GetMpiTasks() << "\n";
oss << "MPI task number:                                  " << phreeqc_rm.GetMpiMyself() << "\n";
oss << "File prefix:                                      " << phreeqc_rm.GetFilePrefix() << "\n";
oss << "Number of grid cells in the user's model:         " << phreeqc_rm.GetGridCellCount() << "\n";
oss << "Number of chemistry cells in the reaction module: " << phreeqc_rm.GetChemistryCellCount() << "\n";
oss << "Number of components for transport:               " << phreeqc_rm.GetComponentCount() << "\n";
oss << "Error handler mode:                               " << phreeqc_rm.GetErrorHandlerMode() << "\n";
phreeqc_rm.OutputMessage(oss.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	void                                      OutputMessage(const std::string &str);
/**
Runs a reaction step for all reaction cells in the reaction module.
Normally, tranport concentrations are transferred to the reaction cells (@ref SetConcentrations) before
reaction calculations are run. The length of time over which kinetic reactions are integrated is set
by @ref SetTimeStep. Other properties that may need to be updated as a result of the transport
calculations include porosity (@ref SetPorosity), saturation (@ref SetSaturation),
temperature (@ref SetTemperature), and pressure (@ref SetPressure).

@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetConcentrations,  @ref SetPorosity,
@ref SetTemperature, @ref SetPressure, @ref SetSaturation, @ref SetTimeStep.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSelectedOutputOn(print_selected_output_on);
status = phreeqc_rm.SetPrintChemistryOn(print_chemistry_on, false, false);
status = phreeqc_rm.SetPorosity(por);             // If porosity changes
status = phreeqc_rm.SetSaturation(sat);           // If saturation changes
status = phreeqc_rm.SetTemperature(temperature);  // If temperature changes
status = phreeqc_rm.SetPressure(pressure);        // If pressure changes
status = phreeqc_rm.SetConcentrations(c);         // Transported concentrations
status = phreeqc_rm.SetTimeStep(time_step);		  // Time step for kinetic reactions
time = time + time_step;
status = phreeqc_rm.SetTime(time);
status = phreeqc_rm.RunCells();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                RunCells(void);
/**
Process an IRM_RESULT return code. If the return code is nonnegative, no action is taken. If the return code is negative,
the return code is decoded and printed as an error message along with the second argument (std::string). On an error,
the method will return the same return code, throw an exception, or exit the program depending on the setting for
@ref SetErrorHandlerMode.
@param result          Return code to be processed.
@param e_string        Error message to be printed in case of an error.
@retval IRM_RESULT     The first argument to the method is returned.
@see                    @ref SetErrorHandlerMode.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.ReturnHandler(irm_result, "Previous method failed.");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root or workers.
 */
	IRM_RESULT                                ReturnHandler(IRM_RESULT result, const std::string &e_string);
/**
Run a PHREEQC input file. The first three arguments determine which IPhreeqc instances will run
the file--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
files that modify the thermodynamic database should be run by all three sets of instances.
Files with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
be run by the workers. Files that contain initial conditions or boundary conditions should
be run by the InitialPhreeqc instance.
@param workers          @a True, the workers will run the file; @a False, the workers will not run the file.
@param initial_phreeqc  @a True, the InitialPhreeqc instance will run the file; @a False, the InitialPhreeqc will not run the file.
@param utility          @a True, the Utility instance will run the file; @a False, the Utility instance will not run the file.
@param chemistry_name   Name of the file to run.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref RunString.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.RunFile(true, true, true, "advect.pqi");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                RunFile(bool workers, bool initial_phreeqc, bool utility,  const std::string & chemistry_name);
/**
Run a PHREEQC input string. The first three arguments determine which
IPhreeqc instances will run
the string--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
strings that modify the thermodynamic database should be run by all three sets of instances.
Strings with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
be run by the workers. Strings that contain initial conditions or boundary conditions should
be run by the InitialPhreeqc instance.
@param workers          @a True, the workers will run the string; @a False, the workers will not run the string.
@param initial_phreeqc  @a True, the InitialPhreeqc instance will run the string; @a False, the InitialPhreeqc will not run the string.
@param utility          @a True, the Utility instance will run the string; @a False, the Utility instance will not run the string.
@param input_string     String containing PHREEQC input.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref RunFile.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::string input = "DELETE; -all";
status = phreeqc_rm.RunString(true, false, true, input.c_str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                RunString(bool workers, bool initial_phreeqc, bool utility, const std::string & input_string);
/**
Print message to the screen.
@param str              String to be printed.
@see                    @ref ErrorMessage, @ref LogMessage, @ref OutputMessage, @ref WarningMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::ostringstream strm;
strm << "Beginning transport calculation "
     <<   phreeqc_rm.GetTime() * phreeqc_rm.GetTimeConversion()
	 << " days\n";
strm << "          Time step             "
     <<   phreeqc_rm.GetTimeStep() * phreeqc_rm.GetTimeConversion()
	 << " days\n";
phreeqc_rm.ScreenMessage(strm.str());
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	void                                      ScreenMessage(const std::string &str);
/**
Select whether to include H2O in the component list.
The concentrations of H and O must be known
accurately (8 to 10 significant digits) for the numerical method of
PHREEQC to produce accurate pH and pe values.
Because most of the H and O are in the water species,
it may be more robust (require less accuracy in transport) to
transport the excess H and O (the H and O not in water) and water.
The default setting (@a true) is to include water, excess H, and excess O as components.
A setting of @a false will include total H and total O as components.
SetComponentH2O must be called before @ref FindComponents.

@param tf               @a True (default), excess H, excess O, and water are included in the component list;
@a False, total H and O are included in the component list.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref FindComponents.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetComponentH2O(true);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetComponentH2O(bool tf);
/**
Use the vector of concentrations (@a c) to set the moles of components in each reaction cell.
The volume of water in a cell is the product of porosity (@ref SetPorosity), saturation (@ref SetSaturation),
and reference volume (@ref SetRepresentativeVolume).
The moles of each component are determined by the volume of water and per liter concentrations.
If concentration units (@ref SetUnitsSolution) are mass fraction, the
density (as specified by @ref SetDensity) is used to convert from mass fraction to per mass per liter.
@param c               Vector of component concentrations. Size of vector is @a ncomps times @a nxyz,
where @a ncomps is the number of components as determined
by @ref FindComponents or @ref GetComponentCount and
@a nxyz is the number of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetDensity, @ref SetPorosity, @ref SetRepresentativeVolume,
@ref SetSaturation, @ref SetUnitsSolution.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> c;
c.resize(nxyz * components.size());
...
AdvectCpp(c, bc_conc, ncomps, nxyz, nbound);
status = phreeqc_rm.SetPorosity(por);             // If porosity changes
status = phreeqc_rm.SetSaturation(sat);           // If saturation changes
status = phreeqc_rm.SetTemperature(temperature);  // If temperature changes
status = phreeqc_rm.SetPressure(pressure);        // If pressure changes
status = phreeqc_rm.SetConcentrations(c);         // Transported concentrations
status = phreeqc_rm.SetTimeStep(time_step);		  // Time step for kinetic reactions
time = time + time_step;
status = phreeqc_rm.SetTime(time);
status = phreeqc_rm.RunCells();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetConcentrations(const std::vector<double> &c);
/**
Select the current selected output by user number. The user may define multiple SELECTED_OUTPUT
data blocks for the workers. A user number is specified for each data block. The value of
the argument @a n_user selects which of the SELECTED_OUTPUT definitions will be used
for selected-output operations.
@param n_user           User number of the SELECTED_OUTPUT data block that is to be used.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetNthSelectedOutputUserNumber, @ref GetSelectedOutput,
@ref GetSelectedOutputColumnCount, @ref GetSelectedOutputCount,
@ref GetSelectedOutputRowCount, @ref GetSelectedOutputHeading,
@ref SetSelectedOutputOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
for (int isel = 0; isel < phreeqc_rm.GetSelectedOutputCount(); isel++)
{
  int n_user = phreeqc_rm.GetNthSelectedOutputUserNumber(isel);
  status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(n_user);
  std::cerr << "Selected output sequence number: " << isel << "\n";
  std::cerr << "Selected output user number:     " << n_user << "\n";
  std::vector<double> so;
  status = phreeqc_rm.GetSelectedOutput(so);
  // Process results here
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT								  SetCurrentSelectedOutputUserNumber(int n_user);
/**
Set the density for each reaction cell. These density values are used
when converting from transported mass-fraction concentrations (@ref SetUnitsSolution) to
produce per liter concentrations during a call to @ref SetConcentrations.
They are also used when converting from reaction-cell concentrations to transport concentrations
(@ref GetConcentrations), if @ref UseSolutionDensityVolume is set to @a false.
@param density          Vector of densities. Size of vector is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetConcentrations, @ref GetGridCellCount, @ref SetConcentrations,
@ref SetUnitsSolution, @ref UseSolutionDensityVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> initial_density;
initial_density.resize(nxyz, 1.0);
phreeqc_rm.SetDensity(initial_density);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetDensity(const std::vector<double> &density);
/**
Set the name of the dump file. It is the name used by @ref DumpModule.
@param dump_name        Name of dump file.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref DumpModule.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetDumpFileName("advection_cpp.dmp");
bool dump_on = true;
bool append = false;
status = phreeqc_rm.DumpModule(dump_on, append);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                SetDumpFileName(const std::string & dump_name);
/**
Set the action to be taken when the reaction module encounters an error.
Options are 0, return to calling program with an error return code (default);
1, throw an exception, in C++, the exception can be caught, for C and Fortran, the program will exit; or
2, attempt to exit gracefully.
@param mode             Error handling mode: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@par C++ Example:
@htmlonly
<CODE>
<PRE>
PhreeqcRM phreeqc_rm(nxyz, nthreads);
IRM_RESULT status;
status = phreeqc_rm.SetErrorHandlerMode(1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetErrorHandlerMode(int mode);
/**
Set the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files.
These files are opened by @ref OpenFiles.
@param prefix           Prefix used when opening the output and log files.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref OpenFiles, @ref CloseFiles.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetFilePrefix("Advect_cpp");
phreeqc_rm.OpenFiles();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                SetFilePrefix(const std::string & prefix);
/**
MPI and C/C++ only. Defines a callback function that allows additional tasks to be done
by the workers. The method @ref MpiWorker contains a loop,
where the workers receive a message (an integer),
run a function corresponding to that integer,
and then wait for another message.
SetMpiWorkerCallbackC allows C or C++ developers to add another function
that responds to additional integer messages by calling developer-defined functions
corresponding to those integers.
@ref MpiWorker calls the callback function when the message number
is not one of the PhreeqcRM message numbers.
Messages are unique integer numbers. PhreeqcRM uses integers in a range
beginning at 0. It is suggested that developers use message numbers starting
at 1000 or higher for their tasks.
The callback function calls a developer-defined function specified
by the message number and then returns to @ref MpiWorker to wait for
another message.
@n@n
In C and C++, an additional pointer can be supplied to find the data necessary to do the task.
A void pointer may be set with @ref SetMpiWorkerCallbackCookie. This pointer
is passed to the callback function through a void pointer argument in addition
to the integer message argument. The pointer may be to a struct or class instance
that provides a number of additional pointers to data. @ref SetMpiWorkerCallbackCookie
must be called by each worker before @ref MpiWorker is called.
@n@n
The motivation for this method is to allow the workers to perform other
tasks, for instance, parallel transport calculations, within the structure
of @ref MpiWorker. The callback function
can be used to allow the workers to receive data, perform transport calculations,
and (or) send results, without leaving the loop of @ref MpiWorker. Alternatively,
it is possible for the workers to return from @ref MpiWorker
by a call to @ref MpiWorkerBreak by root. The workers could then call
subroutines to receive data, calculate transport, and send data,
and then resume processing PhreeqcRM messages from root with another
call to @ref MpiWorker.
@param fcn              A function that returns an integer and has an integer argument
and a void * argument.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref MpiWorker, @ref MpiWorkerBreak,
@ref SetMpiWorkerCallbackCookie.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
Code executed by root:
// root calls a function that will involve the workers
int istatus = do_something(&comm);

Code executed by workers:
phreeqc_rm.SetMpiWorkerCallbackC(worker_tasks_cc);
phreeqc_rm.SetMpiWorkerCallbackCookie(&comm);
phreeqc_rm.MpiWorker();

Code executed by root and workers:
int do_something(void *cookie)
{
	int method_number = 1000;
	MP_TYPE *comm = (MP_TYPE *) cookie;
	int mpi_tasks, mpi_myself, worker_number;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	std::stringstream msg;
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, *comm);
		fprintf(stderr, "I am root.\n");
		for (int i = 1; i < mpi_tasks; i++)
		{
			MPI_Status status;
			MPI_Recv(&worker_number, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			fprintf(stderr, "Recieved data from worker number %d.\n", worker_number);
		}
	}
	else
	{
		MPI_Send(&mpi_myself, 1, MPI_INT, 0, 0, *comm);
	}
	return 0;
}

Code called by workers from method MpiWorker:
int worker_tasks_cc(int *task_number, void * cookie)
{
	if (*task_number == 1000)
	{
		do_something(cookie);
	}
	return 0;
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by workers, before call to @ref MpiWorker.
 */
	IRM_RESULT								  SetMpiWorkerCallbackC(int (*fcn)(int *method, void * cookie));
/**
MPI and C/C++ only. Defines a void pointer that can be used by
C and C++ functions called from the callback function (@ref SetMpiWorkerCallbackC)
to locate data for a task. The C callback function
that is registered with @ref SetMpiWorkerCallbackC has
two arguments, an integer message to identify a task, and a void
pointer. SetMpiWorkerCallbackCookie sets the value of the
void pointer that is passed to the callback function.
The void pointer may be a pointer to a struct of class instance that
contains additonal pointers to data.
@param cookie           Void pointer that can be used by subroutines called from the callback function
to locate data needed to perform a task.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref MpiWorker, @ref MpiWorkerBreak,
@ref SetMpiWorkerCallbackC.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
Code executed by root:
// root calls a function that will involve the workers
int istatus = do_something(&comm);

Code executed by workers:
phreeqc_rm.SetMpiWorkerCallbackC(worker_tasks_cc);
phreeqc_rm.SetMpiWorkerCallbackCookie(&comm);
phreeqc_rm.MpiWorker();

Code executed by root and workers:
int do_something(void *cookie)
{
	int method_number = 1000;
	MP_TYPE *comm = (MP_TYPE *) cookie;
	int mpi_tasks, mpi_myself, worker_number;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myself);
	std::stringstream msg;
	if (mpi_myself == 0)
	{
		MPI_Bcast(&method_number, 1, MPI_INT, 0, *comm);
		fprintf(stderr, "I am root.\n");
		for (int i = 1; i < mpi_tasks; i++)
		{
			MPI_Status status;
			MPI_Recv(&worker_number, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			fprintf(stderr, "Recieved data from worker number %d.\n", worker_number);
		}
	}
	else
	{
		MPI_Send(&mpi_myself, 1, MPI_INT, 0, 0, *comm);
	}
	return 0;
}

Code called by workers from method MpiWorker:
int worker_tasks_cc(int *task_number, void * cookie)
{
	if (*task_number == 1000)
	{
		do_something(cookie);
	}
	return 0;
}
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by workers, before call to @ref MpiWorker.
 */
	IRM_RESULT								  SetMpiWorkerCallbackCookie(void * cookie);
/**
MPI and Fortran only. Defines a callback function that allows additional tasks to be done
by the workers. See documentation of PhreeqcRM for Fortran, method SetMpiWorkerCallback.
 */
	IRM_RESULT								  SetMpiWorkerCallbackFortran(int (*fcn)(int *method));
/**
Sets the property for partitioning solids between the saturated and unsaturated
parts of a partially saturated cell.

The option is intended to be used by saturated-only
flow codes that allow a variable water table.
The value has meaning only when saturations
less than 1.0 are encountered. The partially saturated cells
may have a small water-to-rock ratio that causes
reactions to proceed differently relative to fully saturated cells.
By setting  @a SetPartitionUZSolids to true, the
amounts of solids and gases are partioned according to the saturation.
If a cell has a saturation of 0.5, then
the water interacts with only half of the solids and gases; the other half is unreactive
until the water table rises. As the saturation in a cell varies,
solids and gases are transferred between the
saturated and unsaturated (unreactive) reservoirs of the cell.
Unsaturated-zone flow and transport codes will probably use the default (false),
which assumes all gases and solids are reactive regardless of saturation.
@param tf       @a True, the fraction of solids and gases available for
reaction is equal to the saturation;
@a False (default), all solids and gases are reactive regardless of saturation.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                @ref GetPartitionUZSolids.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
phreeqc_rm.SetPartitionUZSolids(false);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetPartitionUZSolids(bool tf);
/**
Set the porosity for each reaction cell.
The volume of water in a reaction cell is the product of porosity, saturation
(@ref SetSaturation), and representative volume (@ref SetRepresentativeVolume).
@param por              Vector of porosities, unitless. Default is 0.1.
Size of vector is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetSaturation, @ref SetRepresentativeVolume, @ref SetSaturation.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> por;
por.resize(nxyz, 0.2);
status = phreeqc_rm.SetPorosity(por);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetPorosity(const std::vector<double> &por);

/**
Set the pressure for each reaction cell. Pressure effects are considered only in three of the
databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.
@param p                Vector of pressures, in atm. Size of vector is @a nxyz,
where @a nxyz is the number of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetTemperature.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> pressure;
pressure.resize(nxyz, 2.0);
phreeqc_rm.SetPressure(pressure);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetPressure(const std::vector<double> &p);
/**
Enable or disable detailed output for each reaction cell.
Printing for a reaction cell will occur only when the
printing is enabled with @ref SetPrintChemistryOn and the @a cell_mask value is 1.
@param cell_mask        Vector of integers. Size of vector is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref GetGridCellCount). A value of 0 will
disable printing detailed output for the cell; a value of 1 will enable printing detailed output for a cell.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetPrintChemistryOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<int> print_chemistry_mask;
print_chemistry_mask.resize(nxyz, 0);
for (int i = 0; i < nxyz/2; i++)
{
  print_chemistry_mask[i] = 1;
}
status = phreeqc_rm.SetPrintChemistryMask(print_chemistry_mask);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetPrintChemistryMask(std::vector<int> & cell_mask);
/**
Set property that enables or disables printing detailed output from reaction calculations
to the output file for a set of cells defined by @ref SetPrintChemistryMask.
The detailed output prints all of the output typical of a PHREEQC reaction calculation,
which includes solution descriptions and the compositions of all other reactants.
The output can be several hundred lines per cell, which can lead to a very
large output file (prefix.chem.txt, @ref OpenFiles).
For the worker instances, the output can be limited to a set of cells
(@ref SetPrintChemistryMask) and, in general, the
amount of information printed can be limited by use of options in the PRINT data block of PHREEQC
(applied by using @ref RunFile or @ref RunString).
Printing the detailed output for the workers is generally used only for debugging,
and PhreeqcRM will run significantly faster
when printing detailed output for the workers is disabled.

@param workers          @a True, enable detailed printing in the worker instances;
@a False, disable detailed printing in the worker instances.
@param initial_phreeqc  @a True, enable detailed printing in the InitialPhreeqc instance;
@a False, disable detailed printing in the InitialPhreeqc instance.
@param utility          @a True, enable detailed printing in the Utility instance;
@a False, disable detailed printing in the Utility instance.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref OpenFiles, @ref RunFile, @ref RunString, @ref SetPrintChemistryMask.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetPrintChemistryOn(false, true, false);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetPrintChemistryOn(bool workers, bool initial_phreeqc, bool utility);
/**
Set the load-balancing algorithm.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a @ref RunCells
calculation. Two algorithms are available; one uses individual times for each cell and
accounts for cells that were not run because
saturation was zero (default), and
the other assigns an average time to all cells.
The methods are similar, but limited testing indicates the default method performs better.
@param tf           @a True, indicates individual cell times are used in rebalancing (default);
@a False, indicates average times are used in rebalancing.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetRebalanceFraction.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetRebalanceByCell(true);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetRebalanceByCell(bool tf);
/**
Sets the fraction of cells that are transferred among threads or processes when rebalancing.
PhreeqcRM attempts to rebalance the load of each thread or process such that each
thread or process takes the same amount of time to run its part of a @ref RunCells
calculation. The rebalancing transfers cell calculations among threads or processes to
try to achieve an optimum balance. @a SetRebalanceFraction
adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
determine the actual number of cell transfers. A value of zero eliminates
load rebalancing. A value less than 1.0 is suggested to slow the approach to the optimum cell
distribution and avoid possible oscillations
when too many cells are transferred at one iteration, requiring reverse transfers at the next iteration.
Default is 0.5.

@param f                Fraction from 0.0 to 1.0.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetRebalanceByCell.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetRebalanceFraction(0.5);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetRebalanceFraction(double f);
/**
Set the representative volume of each reaction cell.
By default the representative volume of each reaction cell is 1 liter.
The volume of water in a reaction cell is determined by the product of the representative volume,
the porosity (@ref SetPorosity), and the saturation (@ref SetSaturation).
The numerical method of PHREEQC is more robust if the water volume for a reaction cell is
within a couple orders of magnitude of 1.0.
Small water volumes caused by small porosities and (or) small saturations (and (or) small representative volumes)
may cause non-convergence of the numerical method.
In these cases, a larger representative volume may help. Note
that increasing the representative volume also increases
the number of moles of the reactants in the reaction cell (minerals, surfaces, exchangers,
and others), which are defined as moles per representative volume.
@param rv              Vector of representative volumes, in liters. Default is 1.0 liter.
Size of array is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetPorosity, @ref SetSaturation.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> rv;
rv.resize(nxyz, 2.0);
status = phreeqc_rm.SetRepresentativeVolume(rv);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetRepresentativeVolume(const std::vector<double> &rv);
/**
Set the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1.
The volume of water in a cell is the product of porosity (@ref SetPorosity), saturation (@a SetSaturation),
and representative volume (@ref SetRepresentativeVolume). As a result of a reaction calculation,
solution properties (density and volume) will change;
the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar volume data to calculate these changes. The methods @ref GetDensity,
@ref GetSolutionVolume, and @ref GetSaturation can be used to account
for these changes in the succeeding transport calculation.
@a SetRepresentativeVolume should be called before initial conditions are defined for the reaction cells.

@param sat              Vector of saturations, unitless. Default 1.0. Size of vector is @a nxyz,
where @a nxyz is the number of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetDensity, @ref GetSaturation, @ref GetSolutionVolume,
@ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> sat;
sat.resize(nxyz, 1.0);
status = phreeqc_rm.SetSaturation(sat);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetSaturation(const std::vector<double> &sat);
/**
Set the property that controls whether messages are written to the screen.
Messages include information about rebalancing during @ref RunCells, and
any messages written with @ref ScreenMessage.

@param tf  @a True, enable screen messages; @a False, disable screen messages. Default is true.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref RunCells, @ref ScreenMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetScreenOn(true);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root.
 */
	IRM_RESULT                                SetScreenOn(bool tf);
/**
Set the property that controls whether selected-output results are available to be retrieved
with @ref GetSelectedOutput. @a True indicates that selected-output results
will be accumulated during @ref RunCells and can be retrieved with @ref GetSelectedOutput;
@a False indicates that selected-output results will not
be accumulated during @ref RunCells.

@param tf  @a True, enable selected output; @a False, disable selected output.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetSelectedOutput, @ref SetPrintChemistryOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSelectedOutputOn(true);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetSelectedOutputOn(bool tf);
/**
Sets the value of the species-save property.
This method enables or disables use of PhreeqcRM with multicomponent-diffusion transport calculations.
By default, concentrations of aqueous species are not saved.
Setting the species-save property to @a true allows
aqueous species concentrations to be retrieved
with @ref GetSpeciesConcentrations, and solution compositions to be set with
@ref SpeciesConcentrations2Module.
@a SetSpeciesSaveOn must be called before calls to @ref FindComponents.

@param save_on          @a True indicates species concentrations are saved;
@a False indicates species concentrations are not saved.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount,
@ref GetSpeciesD25,
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SpeciesConcentrations2Module.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	IRM_RESULT                                SetSpeciesSaveOn(bool save_on);
/**
Set the temperature for each reaction cell. If @a SetTemperature is not called,
worker solutions will have temperatures as defined by initial conditions
(@ref InitialPhreeqc2Module and @ref InitialPhreeqcCell2Module).

@param t                Vector of temperatures, in degrees C.
Size of vector is @a nxyz, where @a nxyz is the number
of grid cells in the user's model (@ref GetGridCellCount).
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetPressure, @ref InitialPhreeqc2Module,
@ref InitialPhreeqcCell2Module, @ref SetPressure, @ref GetTemperature.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
std::vector<double> temperature;
temperature.resize(nxyz, 20.0);
phreeqc_rm.SetTemperature(temperature);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers.
 */
	IRM_RESULT                                SetTemperature(const std::vector<double> &t);
/**
Set current simulation time for the reaction module.
@param time             Current simulation time, in seconds.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetTimeStep, @ref SetTimeConversion.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
time += time_step;
status = phreeqc_rm.SetTime(time);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetTime(double time);
/**
Set a factor to convert from seconds to user time units. Factor times seconds produces user time units.

@param conv_factor      Factor to convert seconds to user time units.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetTime, @ref SetTimeStep.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
double time_conversion = 1.0 / 86400;
status = phreeqc_rm.SetTimeConversion(time_conversion);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetTimeConversion(double conv_factor);
/**
Set current time step for the reaction module. This is the length
of time over which kinetic reactions are integrated.

@param time_step        Time step, in seconds.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetTime, @ref SetTimeConversion.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
time_step = 86400.;
status = phreeqc_rm.SetTimeStep(time_step);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetTimeStep(double time_step);
/**
Sets input units for exchangers.
In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
@a SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.

If a single EXCHANGE definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of exchangers will be the same regardless of porosity. 
For option 1, the number of moles of exchangers will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of exchangers will vary directly with rock volume and inversely with porosity.

@param option           Units option for exchangers: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetUnitsExchange, @ref InitialPhreeqc2Module, @ref InitialPhreeqcCell2Module,
@ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetUnitsExchange(1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetUnitsExchange(int option);
/**
Set input units for gas phases.
In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
@a SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single GAS_PHASE definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of a gas component will be the same regardless of porosity. 
For option 1, the number of moles of a gas component will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a gas component will vary directly with rock volume and inversely with porosity.

@param option           Units option for gas phases: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetUnitsGasPhase, @ref InitialPhreeqc2Module, @ref InitialPhreeqcCell2Module,
@ref SetPorosity, @ref SetRepresentativeVolume.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetUnitsGasPhase(1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetUnitsGasPhase(int option);
/**
Set input units for kinetic reactants.

In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
@a SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single KINETICS definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of kinetic reactants will be the same regardless of porosity. 
For option 1, the number of moles of kinetic reactants will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of kinetic reactants will vary directly with rock volume and inversely with porosity.

Note that the volume of water in a cell in the reaction module is equal to the product of
porosity (@ref SetPorosity), the saturation (@ref SetSaturation), and representative volume (@ref
SetRepresentativeVolume), which is usually less than 1 liter. It is important to write the RATES
definitions for homogeneous (aqueous) kinetic reactions to account for the current volume of
water, often by calculating the rate of reaction per liter of water and multiplying by the volume
of water (Basic function SOLN_VOL). 

Rates that depend on surface area of solids, are not dependent
on the volume of water. However, it is important to get the correct surface area for the kinetic
reaction. To scale the surface area with the number of moles, the specific area (m^2 per mole of reactant) 
can be defined as a parameter (KINETICS; -parm), which is multiplied by the number of moles of 
reactant (Basic function M) in RATES to obtain the surface area.

@param option           Units option for kinetic reactants: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetUnitsKinetics, @ref InitialPhreeqc2Module, @ref InitialPhreeqcCell2Module,
@ref SetPorosity, @ref SetRepresentativeVolume, @ref SetSaturation.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetUnitsKinetics(1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetUnitsKinetics(int option);
/**
Set input units for pure phase assemblages (equilibrium phases).
In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
@a SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.

If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of a mineral will be the same regardless of porosity. 
For option 1, the number of moles of a mineral will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a mineral will vary directly with rock volume and inversely with porosity.

@param option           Units option for equilibrium phases: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetUnitsPPassemblage, @ref InitialPhreeqc2Module, @ref InitialPhreeqcCell2Module,
@ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetUnitsPPassemblage(1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetUnitsPPassemblage(int option);
/**
Solution concentration units used by the transport model.
Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
PHREEQC defines solutions by the number of moles of each
element in the solution.

To convert from mg/L to moles
of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
multiplied by the solution volume,
which is the product of porosity (@ref SetPorosity), saturation (@ref SetSaturation),
and representative volume (@ref SetRepresentativeVolume).
To convert from mol/L to moles
of element in the representative volume of a reaction cell, mol/L is
multiplied by the solution volume.
To convert from mass fraction to moles
of element in the representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density
(@ref SetDensity) and
multiplied by the solution volume.

To convert from moles
of element in the representative volume of a reaction cell to mg/L, the number of moles of an element is divided by the
solution volume resulting in mol/L, and then converted to mg/L.
To convert from moles
of element in a cell to mol/L,  the number of moles of an element is divided by the
solution volume resulting in mol/L.
To convert from moles
of element in a cell to mass fraction, the number of moles of an element is converted to kg and divided
by the total mass of the solution.
Two options are available for the volume and mass of solution
that are used in converting to transport concentrations: (1) the volume and mass of solution are
calculated by PHREEQC, or (2) the volume of solution is the product of porosity (@ref SetPorosity),
saturation (@ref SetSaturation), and representative volume (@ref SetRepresentativeVolume),
and the mass of solution is volume times density as defined by @ref SetDensity.
Which option is used is determined by @ref UseSolutionDensityVolume.

@param option           Units option for solutions: 1, 2, or 3, default is 1, mg/L.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref SetDensity, @ref SetPorosity, @ref SetRepresentativeVolume, @ref SetSaturation,
@ref UseSolutionDensityVolume.

@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetUnitsSolution(2);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetUnitsSolution(int option);
/**
Set input units for solid-solution assemblages.
In PHREEQC, solid solutions are defined by moles of each component (@a Mp).
@a SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@ P)*RV.

If a single SOLID_SOLUTION definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of a solid-solution component will be the same regardless of porosity. 
For option 1, the number of moles of a solid-solution component will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of a solid-solution component will vary directly with rock volume and inversely with porosity.

@param option           Units option for solid solutions: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetUnitsSSassemblage, @ref InitialPhreeqc2Module, @ref InitialPhreeqcCell2Module,
@ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetUnitsSSassemblage(1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetUnitsSSassemblage(int option);
/**
Set input units for surfaces.
In PHREEQC input, surfaces are defined by moles of surface sites (@a Mp).
@a SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (@a Mc)
is calculated from the input value (@a Mp).

Options are
0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref SetRepresentativeVolume);
1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref SetPorosity); or
2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.

If a single SURFACE definition is used for cells with different initial porosity, 
   the three options scale quite differently. 
For option 0, the number of moles of surface sites will be the same regardless of porosity. 
For option 1, the number of moles of surface sites will be vary directly with porosity and inversely with rock volume. 
For option 2, the number of moles of surface sites will vary directly with rock volume and inversely with porosity.

@param option           Units option for surfaces: 0, 1, or 2.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref GetUnitsSurface, @ref InitialPhreeqc2Module, @ref InitialPhreeqcCell2Module,
@ref SetPorosity, @ref SetRepresentativeVolume.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetUnitsSurface(1);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT                                SetUnitsSurface(int option);
/**
Set solution concentrations in the reaction cells
based on the vector of aqueous species concentrations (@a species_conc).
This method is intended for use with multicomponent-diffusion transport calculations,
and @ref SetSpeciesSaveOn must be set to @a true.
The list of aqueous species is determined by @ref FindComponents and includes all
aqueous species that can be made from the set of components.
The method determines the total concentration of a component
by summing the molarities of the individual species times the stoichiometric
coefficient of the element in each species.
Solution compositions in the reaction cells are updated with these component concentrations.

@param species_conc     Vector of aqueous species concentrations. Dimension of the array is @a nspecies times @a nxyz,
where  @a nspecies is the number of aqueous species (@ref GetSpeciesCount),
and @a nxyz is the number of user grid cells (@ref GetGridCellCount).
Concentrations are moles per liter.
@retval IRM_RESULT      0 is success, negative is failure (See @ref DecodeError).
@see                    @ref FindComponents, 
@ref GetSpeciesConcentrations, 
@ref GetSpeciesCount,
@ref GetSpeciesD25,
@ref GetSpeciesLog10Gammas, 
@ref GetSpeciesNames, 
@ref GetSpeciesSaveOn, 
@ref GetSpeciesStoichiometry, 
@ref GetSpeciesZ,
@ref SetSpeciesSaveOn.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
status = phreeqc_rm.SetSpeciesSaveOn(true);
int ncomps = phreeqc_rm.FindComponents();
int nspecies = phreeqc_rm.GetSpeciesCount();
std::vector<double> c;
status = phreeqc_rm.GetSpeciesConcentrations(c);
...
SpeciesAdvectCpp(c, bc_conc, nspecies, nxyz, nbound);
status = phreeqc_rm.SpeciesConcentrations2Module(c);
status = phreeqc_rm.RunCells();
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
	IRM_RESULT								  SpeciesConcentrations2Module(std::vector<double> & species_conc);
/**
Determines the volume and density to use when converting from the reaction-cell concentrations
to transport concentrations (@ref GetConcentrations).
Two options are available to convert concentration units:
(1) the density and solution volume calculated by PHREEQC are used, or
(2) the specified density (@ref SetDensity)
and solution volume are determined by the product of
saturation (@ref SetSaturation), porosity (@ref SetPorosity),
and representative volume (@ref SetRepresentativeVolume).
Transport models that consider density-dependent flow will probably use the
PHREEQC-calculated density and solution volume (default),
whereas transport models that assume constant-density flow will probably use
specified values of density and solution volume.
Only the following databases distributed with PhreeqcRM have molar-volume information
needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
Density is only used when converting to or from transport units of mass fraction.

@param tf          @a True indicates that the solution density and volume as
calculated by PHREEQC will be used to calculate concentrations.
@a False indicates that the solution density set by @ref SetDensity and the volume determined by the
product of  @ref SetSaturation, @ref SetPorosity, and @ref SetRepresentativeVolume,
will be used to calculate concentrations retrieved by @ref GetConcentrations.
@see                    @ref GetConcentrations,  @ref SetDensity,
@ref SetPorosity, @ref SetRepresentativeVolume, @ref SetSaturation.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
phreeqc_rm.UseSolutionDensityVolume(false);
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root, workers must be in the loop of @ref MpiWorker.
 */
void UseSolutionDensityVolume(bool tf);
/**
Print a warning message to the screen and the log file.

@param warnstr          String to be printed.
@see                    @ref OpenFiles, @ref ErrorMessage, @ref LogMessage, @ref OutputMessage, @ref ScreenMessage.
@par C++ Example:
@htmlonly
<CODE>
<PRE>
phreeqc_rm.WarningMessage("Parameter is out of range, using default");
</PRE>
</CODE>
@endhtmlonly
@par MPI:
Called by root and (or) workers; only root writes to the log file.
 */
	void                                      WarningMessage(const std::string &warnstr);

	// Utilities
	static std::string                        Char2TrimString(const char * str, size_t l = 0);
	static bool                               FileExists(const std::string &name);
	static void                               FileRename(const std::string &temp_name, const std::string &name,
		                                           const std::string &backup_name);
	static IRM_RESULT                         Int2IrmResult(int r, bool positive_ok);
protected:
	IRM_RESULT                                CellInitialize(
		                                          int i,
		                                          int n_user_new,
		                                          int *initial_conditions1,
		                                          int *initial_conditions2,
		                                          double *fraction1,
		                                          std::set<std::string> &error_set);
	IRM_RESULT                                CheckCells();
	int                                       CheckSelectedOutput();
    //void                                      Collapse2Nchem(double *d_in, double *d_out);
    //void                                      Collapse2Nchem(int *i_in, int *i_out);
	IPhreeqc *                                Concentrations2UtilityH2O(std::vector<double> &c_in,
		                                           std::vector<double> &t_in, std::vector<double> &p_in);
	IPhreeqc *                                Concentrations2UtilityNoH2O(std::vector<double> &c_in,
		                                           std::vector<double> &t_in, std::vector<double> &p_in);
	void                                      Concentrations2Solutions(int n, std::vector<double> &c);
	void                                      Concentrations2SolutionsH2O(int n, std::vector<double> &c);
	void                                      Concentrations2SolutionsNoH2O(int n, std::vector<double> &c);
	void                                      cxxSolution2concentration(cxxSolution * cxxsoln_ptr, std::vector<double> & d, double v, double dens);
	void                                      cxxSolution2concentrationH2O(cxxSolution * cxxsoln_ptr, std::vector<double> & d, double v, double dens);
	void                                      cxxSolution2concentrationNoH2O(cxxSolution * cxxsoln_ptr, std::vector<double> & d, double v, double dens);
    void                                      GatherNchem(std::vector<double> &source, std::vector<double> &destination);
	cxxStorageBin &                           Get_phreeqc_bin(void) {return *this->phreeqc_bin;}
	IRM_RESULT                                HandleErrorsInternal(std::vector< int > & r);
	void                                      PartitionUZ(int n, int iphrq, int ihst, double new_frac);
	void                                      RebalanceLoad(void);
	void                                      RebalanceLoadPerCell(void);
	IRM_RESULT                                RunCellsThread(int i);
	IRM_RESULT                                RunFileThread(int n);
	IRM_RESULT                                RunStringThread(int n, std::string & input);
	IRM_RESULT                                RunCellsThreadNoPrint(int n);
	void                                      Scale_solids(int n, int iphrq, double frac);
	void                                      ScatterNchem(double *d_array);
	void                                      ScatterNchem(int *i_array);
	void                                      ScatterNchem(std::vector<double> &source, std::vector<double> &destination);
	void                                      ScatterNchem(std::vector<int> &source, std::vector<int> &destination);
	IRM_RESULT                                SetChemistryFileName(const char * prefix = NULL);
	IRM_RESULT                                SetDatabaseFileName(const char * db = NULL);
	void                                      SetEndCells(void);
	void                                      SetEndCellsHeterogeneous(void);
	double                                    TimeStandardTask(void);
	IRM_RESULT                                TransferCells(cxxStorageBin &t_bin, int old, int nnew);
	IRM_RESULT                                TransferCellsUZ(std::ostringstream &raw_stream, int old, int nnew);

private:
	//IRM_RESULT                                SetGeneric(std::vector<double> &destination, int newSize, const std::vector<double> &origin, int mpiMethod, const std::string &name, const double newValue = 0.0);
	IRM_RESULT                                SetGeneric(const std::vector<double> &source, std::vector<double> &destination_root, std::vector<double> &destination_worker, int mpiMethod, const std::string &name);
protected:

#if defined(_MSC_VER)
/* disable warning C4251: 'identifier' : class 'type' needs to have dll-interface to be used by clients of class 'type2' */
#pragma warning(disable:4251)
#endif

	bool component_h2o;                      // true: use H2O, excess H, excess O, and charge;
	                                         // false total H, total O, and charge
	std::string database_file_name;
	std::string chemistry_file_name;
	std::string dump_file_name;
	std::string file_prefix;
	cxxStorageBin * phreeqc_bin;
	int mpi_myself;
	int mpi_tasks;
	std::vector <std::string> components;	// list of components to be transported
	std::vector <double> gfw;				// gram formula weights converting mass to moles (1 for each component)
	double gfw_water;						// gfw of water
	bool partition_uz_solids;
	int nxyz;								// number of nodes
	int count_chemistry;					// number of cells for chemistry
	double time;						    // time from transport, sec
	double time_step;					    // time step from transport, sec
	double time_conversion;					// time conversion factor, multiply to convert to preferred time unit for output
	std::vector <double> old_saturation_root;	// saturation fraction from previous step
	std::vector <double> old_saturation_worker;
	std::vector<double> saturation_root;	    // nxyz saturation fraction
	std::vector<double> saturation_worker;	    // nchem on workers saturation fraction
	std::vector<double> pressure_root;			// nxyz on root current pressure
	std::vector<double> pressure_worker;		// nchem on workers current pressure
	std::vector<double> rv_root;		        // nxyz on root representative volume
	std::vector<double> rv_worker;		        // nchem on workers representative volume
	std::vector<double> porosity_root;		    // nxyz porosity
	std::vector<double> porosity_worker;	    // nchem on workers porosity
	std::vector<double> tempc_root;             // nxyz on root temperature Celsius 
	std::vector<double> tempc_worker;		    // nchem on workers temperature Celsius 
	std::vector<double> density_root;			// nxyz density
	std::vector<double> density_worker;			// nchem on workers density
	std::vector<double> solution_volume_root;   // nxyz on root solution volume
	std::vector<double> solution_volume_worker;	// nchem on workers solution_volume 
	std::vector<int> print_chem_mask_root;		// nxyz print flags for output file
	std::vector<int> print_chem_mask_worker;	// nchem print flags for output file
	bool rebalance_by_cell;                 // rebalance method 0 std, 1 by_cell
	double rebalance_fraction;			    // parameter for rebalancing process load for parallel
	int units_Solution;                     // 1 mg/L, 2 mol/L, 3 kg/kgs
	int units_PPassemblage;                 // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	int units_Exchange;                     // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	int units_Surface;                      // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	int units_GasPhase;                     // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	int units_SSassemblage;                 // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	int units_Kinetics;                     // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
	std::vector <int> forward_mapping_root;				    // mapping from nxyz cells to count_chem chemistry cells
	std::vector <std::vector <int> > backward_mapping;	    // mapping from count_chem chemistry cells to nxyz cells
	bool use_solution_density_volume;

	// print flags
	std::vector<bool> print_chemistry_on;	// print flag for chemistry output file
	bool selected_output_on;				// create selected output

	int error_count;
	int error_handler_mode;                 // 0, return code; 1, throw; 2 exit;
	bool need_error_check;
	std::string phreeqcrm_error_string;

	// threading
	int nthreads;
	std::vector<IPhreeqcPhast *> workers;
	std::vector<int> start_cell;
	std::vector<int> end_cell;
	PHRQ_io *phreeqcrm_io;
	bool delete_phreeqcrm_io;

	// mpi
#ifdef USE_MPI
	MPI_Comm phreeqcrm_comm;                                       // MPI communicator
#endif
	int (*mpi_worker_callback_fortran) (int *method);
	int (*mpi_worker_callback_c) (int *method, void *cookie);
	void *mpi_worker_callback_cookie;

	// mcd
	bool species_save_on;
	std::vector <std::string> species_names;
	std::vector <double> species_z;
	std::vector <double> species_d_25;
	std::vector <cxxNameDouble> species_stoichiometry;
	std::map<int, int> s_num2rm_species_num;
	std::vector<double> standard_task_vector;   // root only

	// reactant lists
	std::vector <std::string> ExchangeSpeciesNamesList;
	std::vector <std::string> ExchangeNamesList;
	std::vector <std::string> SurfaceSpeciesNamesList;
	std::vector <std::string> SurfaceTypesList;
	std::vector <std::string> SurfaceNamesList;

	std::vector <std::string> EquilibriumPhasesList;
	std::vector <std::string> GasComponentsList;
	std::vector <std::string> KineticReactionsList;
	std::vector <std::string> SolidSolutionComponentsList;
	std::vector <std::string> SolidSolutionNamesList;
	std::vector <std::string> SINamesList;

private:
	//friend class RM_interface;
	static std::map<size_t, PhreeqcRM*> Instances;
	static size_t InstancesIndex;

#if defined(_MSC_VER)
/* reset warning C4251 */
#pragma warning(default:4251)
#endif

};
#endif // !defined(PHREEQCRM_H_INCLUDED)
