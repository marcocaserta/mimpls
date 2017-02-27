/***************************************************************************
 *   copyright (C) 2005 by Marco Caserta                                   *
 *   marco.caserta@itesm.mx                                                *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*! \file options.cpp 
  \brief Read options from command line.

  Options are:
  - -f : problem instance file [default = NONE]
  - -t : wall-clock time limit for execution
  - -h : help (list of all options)
*/

#include <iostream>
#include <cstdlib>
/**********************************************************/
#define   TIME_LIMIT_def  180   //!< default wall-clock time limit
#define   HOP_def  6		//!< default hop value
/**********************************************************/

using namespace std;

extern char* _FILENAME; 	//!< name of the instance file
extern int time_limit;		//!< wall-clock time limit
extern int maxI;		//!< HOP constraint
// options used by random file generator
extern int rgItems;		//!< number of items for random generated instance
extern int rgPeriods;		//!< number of periods for random generated instance

/// Parse command line options
int parseOptions(int argc, char* argv[])
{
   bool setFile = false;
   time_limit   = TIME_LIMIT_def;
   maxI         = HOP_def;
   cout <<endl << "CLSP - DW v1.0 -- MC 2011(c)" << endl;
   if (argc == 1)
   {
      cout << "No options specified. Try -h " << endl;
      return -1;
   }  
  
 
   int i = 0;
   while (++i < argc)
   {
      const char *option = argv[i];
      if (*option != '-')
	 return i;
      else if (*option == '\0')
	 return i;
      else if (*option == '-')
      {
	 switch (*++option)
	 {
	    case '\0':
	       return i + 1;
	    case 'f':
	       _FILENAME = argv[i+1];
	       setFile = true;
	       i++;
	       break;
	    case 't':
	       time_limit = atol(argv[i+1]);
	       i++;
	       break;
	    case 'i':
	       maxI = atol(argv[i+1]);
	       i++;
	       break;
	    case 'h':
	       cout << "OPTIONS :: " << endl;
	       cout << "-f : problem instance file" << endl;
	       cout << "-t : time limit (real)" << endl;
	       cout << "-i : HOP constraint" << endl;
	       cout << endl;
	       return -1;
	 }
      }
   }
 
   if (setFile)
      return 0;
   else
   {
      cout <<"Option -f is mandatory. Try ./dw -h" << endl;
      return -1;
   }
}
