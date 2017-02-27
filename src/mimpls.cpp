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

/*! \file clsp.cpp
  \brief Solving MIP using cplex.
  */

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

typedef IloArray <IloNumVarArray> TwoD;

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
extern INSTANCE inp;

extern int maxI;		// hop constraint

/// Define math model for the CLSP.
void define_ilo_model(INSTANCE & inp, IloModel model, TwoD y_ilo, TwoD x_ilo, 
        TwoD s_ilo, IloNumVarArray sI)
{
    IloEnv env = model.getEnv();

    // capacity constraints
    for (int t = 0; t < inp.nP; t++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.nI; j++)
            sum += (inp.a[j][t]*x_ilo[j][t] + inp.m[j][t]*y_ilo[j][t]);

        model.add(sum <= inp.cap[t]);
        sum.end();
    }

    // demand constraints
    for (int j = 0; j < inp.nI; j++)
    {
        // first period
        model.add(x_ilo[j][0] + sI[j] - inp.d[j][0] - s_ilo[j][0] == 0);

        // t = 2, ..., nPeriods
        for (int t = 1; t < inp.nP-1; t++)
            model.add(x_ilo[j][t] + s_ilo[j][t-1] - inp.d[j][t] - s_ilo[j][t] == 0); 

        // last period (ending inventory equal to zero)
        model.add(x_ilo[j][inp.nP-1] + s_ilo[j][inp.nP-2] - inp.d[j][inp.nP-1] == 0);
    }

    // x--y constraints
    for (int j = 0; j < inp.nI; j++)
        for (int t = 0; t < inp.nP; t++)
            model.add(x_ilo[j][t] <= inp.max_prod[j][t]*y_ilo[j][t]);


    // hop constraints
    // in each period, the level of inventory must be <= the sum of the
    // demands of the following \c maxI periods 
    for (int j = 0; j < inp.nI; j++)
        for (int t = 0; t <= inp.nP - maxI ; t++)
        {
            IloExpr sum(env);
            for (int k = t+1; k <= t + maxI; k++)
                sum += inp.d[j][k];

            //cout << "Hop(" << j << "," << t << ") = " << sum << endl;
            model.add(s_ilo[j][t] - sum <= 0);
            // int abc;
            // cin >> abc;
        }

    // objective function
    IloExpr totCost(env);
    for (int j = 0; j < inp.nI; j++)
    {
        totCost += sI[j]*10000;
        for (int t = 0; t < inp.nP; t++)
            totCost += (inp.f[j][t]*y_ilo[j][t] + inp.c[j][t]*x_ilo[j][t] + inp.h[j][t]*s_ilo[j][t]);
        //totCost += (f[j][t]*y_ilo[j][t] + h[j][t]*s_ilo[j][t]);
    }

    model.add(IloMinimize(env, totCost));
}

/// Solve MIP model
double solve_CLSP(IloModel model, IloCplex cplex)
{
    try
    {
        IloEnv env = model.getEnv();
        //cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::MIPInterval,5000);
        cplex.setParam(IloCplex::MIPDisplay,2);

        // Optimize the problem and obtain solution.
        if ( !cplex.solve() ) 
        {
            env.error() << "Failed to optimize MIP" << endl;
            throw(-1);
        }
        return cplex.getObjValue();
    }
    catch (...) 
    {
        cout << "Exception caught " << endl;
        return -100;
    }
}


void print_results(INSTANCE & inp, IloCplex cplex, TwoD y_ilo, TwoD x_ilo, 
        TwoD s_ilo, IloNumVarArray sI)
{
    double z_cplex = cplex.getObjValue();
    cout << "z* = " << z_cplex << endl;

    cout << setw(8) << "t";
    for (int t = 0; t < inp.nP; t++)      
        cout << setw(4) << t+1;
    cout << endl;

    cout << "=======================================================================" << endl;

    for (int j = 0; j < inp.nI; j++)
    {
        cout << setw(8) << "d_t";
        for (int t = 0; t < inp.nP; t++)      
            cout << setw(4) << inp.d[j][t];
        cout << endl;
        cout << "---------------------------------------------------------------------" << endl;

        cout << "(" << setw(2) << j+1 <<") :: ";

        for (int t = 0; t < inp.nP; t++)
            cout << setw(4) << (int)cplex.getValue(y_ilo[j][t]);
        cout << endl;
        cout << setw(8) << ".";
        for (int t = 0; t < inp.nP; t++)
            if (cplex.getValue(x_ilo[j][t]) < 0.0001)
                cout << setw(4) << "0";
            else
                cout << setw(4) << cplex.getValue(x_ilo[j][t]);
        cout << endl;
        cout << setw(4) << ".";
        cout << setw(4) << cplex.getValue(sI[j]);
        for (int t = 0; t < inp.nP; t++)
            if (cplex.getValue(s_ilo[j][t]) < 0.0001)
                cout << setw(4) << "0";
            else
                cout << setw(4) << cplex.getValue(s_ilo[j][t]);
        cout << endl;
        cout << "---------------------------------------------------------------------" << endl;

    }

    // verify capacity
    for (int t = 0; t < inp.nP; t++)
    {
        double totCap = 0.0;
        for (int j = 0; j < inp.nI; j++)
            if (cplex.getValue(y_ilo[j][t]) > 0.1)
                totCap += (inp.m[j][t] + inp.a[j][t]*cplex.getValue(x_ilo[j][t]));
        cout << "t = " << setw(3) << t+1 << " :: " << setw(5) << totCap << "/" << setw(5) << inp.cap[t] << endl;
    }

    // verify z value
    double z = 0.0;
    for (int j = 0; j < inp.nI; j++)
        for (int t = 0; t < inp.nP; t++)
        {
            if(cplex.getValue(y_ilo[j][t]) > 0.1)
                z += inp.f[j][t];
            if(cplex.getValue(s_ilo[j][t]) > 0.1)
                z+= cplex.getValue(s_ilo[j][t])*inp.h[j][t];
        }
    cout << "z verified = " << z << endl;
}

/// Solve MIP to optimality using cplex
double solve_MIP(INSTANCE & inp)
{
    IloEnv env;			//!< cplex variables
    IloNumVarArray sI(env, inp.nI, 0, IloInfinity, ILOFLOAT);
    TwoD y_ilo(env, inp.nI);
    TwoD x_ilo(env, inp.nI);
    TwoD s_ilo(env, inp.nI);
    for (int j = 0; j < inp.nI; j++)
    {
        y_ilo[j] = IloNumVarArray(env, inp.nP, 0, 1, ILOINT);
        x_ilo[j] = IloNumVarArray(env, inp.nP, 0, IloInfinity, ILOFLOAT);
        s_ilo[j] = IloNumVarArray(env, inp.nP, 0, IloInfinity, ILOFLOAT);
    }

    IloModel model(env);
    define_ilo_model(inp, model, y_ilo, x_ilo, s_ilo, sI);
    IloCplex cplex(model);

    double z = solve_CLSP(model, cplex);
    print_results(inp, cplex, y_ilo, x_ilo, s_ilo, sI);

    env.end();

    return z;
}






