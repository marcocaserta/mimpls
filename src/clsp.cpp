/***************************************************************************
 *   copyright (C) 2011 by Marco Caserta                                   *
 *   marco dot caserta at uni-hamburg dot de                               *
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

/*! \mainpage Capacitated Lot Sizing Problem with Setups (CLST)

  Algorithm for Multi Item Multi Period Lot Sizing Problem with Setups solved
  using Dantzig-Wolfe decomposition in the context of a corridor method
  algorithm.

  \author Marco Caserta 2011 (c) (<a href="mailto:marco.caserta@uni-hamburg.de">marco.caserta@uni-hamburg.de</a>) 
  \version v. 1.0.0
  \date Begins: 10.03.11
  \date Ends:

  The model is:

  \f[ 
  \mbox{(MCLS):} \quad \left|
  \begin{array}{rll}
  \min  & \displaystyle\sum_{i=1}^P\sum_{t=1}^T \left( c_{jt}y_{jt}+f_{jt}x_{jt} +h_{jt}s_{jt}   \right)\\
  \mbox{s.t.}   & \displaystyle\sum_{i=1}^P \left( a_{jt}x_{jt}+m_{jt}y_{jt}
  \right) \leq P_t  & t=1,2,\ldots, T \\
  &  y_{jt}+s_{j,t-1}-s_{jt}=d_{jt} & j=1,2,\ldots,P;
  t=1,2,\ldots, T\\    
  &  y_{jt}\leq M x_{jt} & j=1,2,\ldots, P; t=1,2,\ldots, T\\
  & y_{jt}, s_{jt}\geq 0, \ x_{jt}\in \{0,1\} & j=1,2,\ldots, P;
  t=1,2,\ldots, T.\\
  \end{array}
  \right.
  \f]

  In this model, \f$P\f$ and \f$T\f$ stand for the number of items and the
  number of periods in the planning horizon, respectively. The model
  presents two groups of continuous variables, which account for the
  production and the inventory level of each item in each period
  (\f$y_{jt}\f$ and \f$s_{jt}\f$, respectively). In addition, a group of binary
  variables \f$x_{jt}\f$ is used to indicate whether any positive quantity of
  item \f$j\f$ is produced during period \f$t\f$. Parameters \f$c_{jt}\f$,
  \f$f_{jt}\f$, and \f$h_{jy}\f$ stand for unitary production cost,
  setup cost, and inventory cost for item \f$i\f$ during period \f$t\f$,
  respectively. Production capacity during period \f$t\f$, e.g., available
  production time, is indicated by \f$b_t\f$. Finally, \f$a_{jt}\f$ and \f$m_{jt}\f$
  represent the 
  unitary production time and the setup time for product \f$j\f$ during
  period \f$t\f$, respectively. The first set of constraints accounts for
  the production capacity limit over the periods and ensures that the
  overall production does not require more than the available capacity.
  The second set of constraints, inventory balance constraints,
  ensures that, for each period and each item, demands are
  satisfied. Finally, the third set of constraints ensures that
  appropriate setup costs are paid if some positive production is taken. 


  \file dw.cpp
  \brief General Implementation of the CLSP

*/
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

typedef IloArray <IloNumVarArray> TwoD;
#include <limits> 


#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <map>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cassert>
#include <vector>

#include "timer.h"
#include "options.h"

double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM  = 2147483647;
const double ZERO      = 0.0e0;
const double EPSI      = 0.00000001;


//=============== VARIABLES DECLARATION ======================================= 
struct INSTANCE {
    double *** dcum;
    double  ** d;
    double  ** c;
    double  ** f;
    double  ** h;
    double  ** a;
    double  ** m;
    double  ** max_prod;
    double   * cap;

    int nI;
    int nP;
};

INSTANCE inp;


timer tTime;			//!< Object clock to measure REAL and VIRTUAL (cpu) time
int time_limit;			//!< Wall-clock time limit
char * _FILENAME;		//!< Instance name file
int maxI;
//=========== END VARIABLES DECLARATION ======================================= 


//=============== FUNCTIONS DECLARATION ======================================= 
void read_problem_data(INSTANCE & inp);
void print_options(INSTANCE & inp);

// forward declaration
double solve_MIP(INSTANCE & inp);
//=============== END FUNCTIONS DECLARATION =================================== 


/************************ MAIN PROGRAM ******************************/
/// Main program
/************************ MAIN PROGRAM ******************************/
int main(int argc, char *argv[])
{


    int err = parseOptions(argc, argv);
    if ( err != 0)
    {
        if (err != -1)
            cout << "Error argument " << err+1 << endl;
        exit(1);
    }

    tTime.resetTime();		// start clock
    read_problem_data(inp);		// read instance data
    print_options(inp);

    double z_opt = solve_MIP(inp);	// solve CLSP to optimality using cplex
    cout << "z* = " << z_opt << " found in " <<  tTime.elapsedTime(timer::REAL) << " seconds." << endl;
    ofstream fout("result.csv", ios::out);
    fout << inp.nI << "\t" << inp.nP << "\t" << maxI << "\t" << z_opt << "\t" 
        << tTime.elapsedTime(timer::REAL) << endl; 
    fout.close();

    exit(123);

}


/// Read Trigeiro instances.
/** Describe the format here.
*/
void read_problem_data(INSTANCE & inp)
{
    double a1, h1, m1, f1;
    int temp;

    ifstream fmimpls(_FILENAME, ios::in);
    if (!fmimpls) {cerr << "Cannot open file " << _FILENAME << endl; exit(1);}
    fmimpls >> inp.nI >> inp.nP;
    fmimpls >> temp;   

    inp.cap = new double[inp.nP];
    for (int t = 0; t < inp.nP; t++)
        inp.cap[t] = temp;

    inp.d = new double*[inp.nI];
    inp.c = new double*[inp.nI];
    inp.f = new double*[inp.nI];
    inp.h = new double*[inp.nI];
    inp.a = new double*[inp.nI];
    inp.m = new double*[inp.nI];


    for (int j = 0; j < inp.nI; j++)
    {
        inp.c[j] = new double[inp.nP];
        inp.f[j] = new double[inp.nP];
        inp.h[j] = new double[inp.nP];
        inp.a[j] = new double[inp.nP];
        inp.m[j] = new double[inp.nP];
        inp.d[j] = new double[inp.nP];
    }

    for (int i = 0; i < inp.nI; i++)
    {
        fmimpls >> a1 >> h1 >> m1 >> f1;
        //cout  << " " << a1 << " " <<  h1 << " " <<  m1 << " " <<  f1 << endl;
        for  (int t = 0; t < inp.nP; t++) 
        {
            //a[i][t] = a1;
            inp.a[i][t] = 1.0;
            inp.h[i][t] = h1;
            inp.m[i][t] = m1;
            inp.f[i][t] = f1;
            inp.c[i][t] = 0.0;
        }
    }

    for  (int t = 0; t < inp.nP; t++) 
        for (int i = 0; i < inp.nI; i++)
            fmimpls >> inp.d[i][t];

    fmimpls.close();

    // compute maximum production per item and period
    inp.max_prod = new double*[inp.nI];
    for (int j = 0; j < inp.nI; j++)
        inp.max_prod[j] = new double[inp.nP];

    for (int t = 0; t < inp.nP; t++)
        for (int j = 0; j < inp.nI; j++)
            inp.max_prod[j][t] = ((inp.cap[t] - inp.m[j][t])/inp.a[j][t]);
}

/// Print instance info and algorithmic parameters.
void print_options(INSTANCE & inp)
{
    cout << "-------------------------------------" << endl;
    cout << "- OPTIONS : " << endl;
    cout << "-------------------------------------" << endl;
    cout << "  DATA FILE   =  " << _FILENAME        << endl;
    cout << "  TIME LIMIT  =  " << time_limit       << endl;
    cout << "  N. ITEMS    =  " << inp.nI           << endl;
    cout << "  N. PERIODS  =  " << inp.nP           << endl;
    cout << "-------------------------------------" <<  endl << endl;   
}


