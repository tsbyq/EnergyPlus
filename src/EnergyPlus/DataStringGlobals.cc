// EnergyPlus Headers
#include <DataStringGlobals.hh>

namespace EnergyPlus {

namespace DataStringGlobals {

	// MODULE INFORMATION:
	//       AUTHOR         Linda K. Lawrie
	//       DATE WRITTEN   September 1997
	//       MODIFIED       na
	//       RE-ENGINEERED  na

	// PURPOSE OF THIS MODULE:
	// This data-only module is a repository for string variables used in parsing
	// "pieces" of EnergyPlus.

	// METHODOLOGY EMPLOYED:
	// na

	// REFERENCES:
	// na

	// OTHER NOTES:
	// na

	// USE STATEMENTS:
	// None!--This module is USEd by other modules; it should not USE anything.

	// Data
	// -only module should be available to other modules and routines.
	// Thus, all variables in this module must be PUBLIC.

	// MODULE PARAMETER DEFINITIONS:
	Fstring const UpperCase( "ABCDEFGHIJKLMNOPQRSTUVWXYZ�����������������������������" );
	Fstring const LowerCase( "abcdefghijklmnopqrstuvwxyz�����������������������������" );
	Fstring const AccentedUpperCase( "�����������������������������" );
	Fstring const AccentedLowerCase( "�����������������������������" );
	Fstring const AllCase( "����������������������������������������������������������ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" );
#ifdef WINDOWS
	Fstring const pathChar( "\\" );
	Fstring const altpathChar( "/" );
#elif defined LINUX
	Fstring const pathChar( "/" );
	Fstring const altpathChar( "\\" );
#elif defined MAC
	Fstring const pathChar( "/" );
	Fstring const altpathChar( "\\" );
#else
	Fstring const pathChar( "\\" ); // default to Windows
	Fstring const altpathChar( "/" );
#endif
	int const PathLimit( 255 );
	Fstring const CharComma( 1, CHAR( 44 ) ); // comma
	Fstring const CharSemicolon( 1, CHAR( 59 ) ); // semicolon
	Fstring const CharTab( 1, CHAR( 9 ) ); // tab
	Fstring const CharSpace( 1, CHAR( 32 ) ); // space

	// DERIVED TYPE DEFINITIONS
	// na

	// INTERFACE BLOCK SPECIFICATIONS
	// na

	// MODULE VARIABLE DECLARATIONS:
	Fstring ProgramPath( PathLimit ); // Path for Program from Energy+.ini
	Fstring CurrentWorkingFolder( PathLimit ); // Current working directory for run
	Fstring FullName( PathLimit + 15 ); // Full name of file to open, including path
	Fstring IDDVerString( 120 ); // Version information from the IDD (line 1)
	Fstring VerString( 120, "EnergyPlus, Version 8.1" ); // String that represents version information
	Fstring MatchVersion( "8.1.0" ); // String to be matched by Version object
	Fstring CurrentDateTime( 40 ); // For printing current date and time at start of run

	//     NOTICE
	//     Copyright � 1996-2014 The Board of Trustees of the University of Illinois
	//     and The Regents of the University of California through Ernest Orlando Lawrence
	//     Berkeley National Laboratory.  All rights reserved.
	//     Portions of the EnergyPlus software package have been developed and copyrighted
	//     by other individuals, companies and institutions.  These portions have been
	//     incorporated into the EnergyPlus software package under license.   For a complete
	//     list of contributors, see "Notice" located in EnergyPlus.f90.
	//     NOTICE: The U.S. Government is granted for itself and others acting on its
	//     behalf a paid-up, nonexclusive, irrevocable, worldwide license in this data to
	//     reproduce, prepare derivative works, and perform publicly and display publicly.
	//     Beginning five (5) years after permission to assert copyright is granted,
	//     subject to two possible five year renewals, the U.S. Government is granted for
	//     itself and others acting on its behalf a paid-up, non-exclusive, irrevocable
	//     worldwide license in this data to reproduce, prepare derivative works,
	//     distribute copies to the public, perform publicly and display publicly, and to
	//     permit others to do so.
	//     TRADEMARKS: EnergyPlus is a trademark of the US Department of Energy.

} // DataStringGlobals

} // EnergyPlus