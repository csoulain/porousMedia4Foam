/*! @file IrmResult.h
 *  @brief Enumeration used to return error codes.
*/
#if !defined(IRMRESULT_H_INCLUDED)
#define IRMRESULT_H_INCLUDED
/*! \brief Enumeration for PhreeqcRM function return codes.
*/
typedef enum {
	IRM_OK            =  0,  /*!< Success */
	IRM_OUTOFMEMORY   = -1,  /*!< Failure, Out of memory */
	IRM_BADVARTYPE    = -2,  /*!< Failure, Invalid VAR type */
	IRM_INVALIDARG    = -3,  /*!< Failure, Invalid argument */
	IRM_INVALIDROW    = -4,  /*!< Failure, Invalid row */
	IRM_INVALIDCOL    = -5,  /*!< Failure, Invalid column */
	IRM_BADINSTANCE   = -6,  /*!< Failure, Invalid rm instance id */
	IRM_FAIL          = -7   /*!< Failure, Unspecified */
} IRM_RESULT;
#endif
