#include "main.h"

int main() {
	// Load an instance:
	// - either from a file by specifying: filename (including path)
	// Standard instances
	string filename = "InstancesBPPMCF/Dataset 1/70-8/70-8-1.txt";
	// string filename = "InstancesBPPMCF/Dataset 2/120-2/120-2-1.txt";
	// string filename = "InstancesBPPMCF/Dataset 3/10-100-4/10-100-4-1.txt";
	// string filename = "InstancesBPPMCF/Dataset 4/50-400-3/50-400-3-1.txt";
	// string filename = "InstancesBPPMCF/Triplets/t60_00.txt";
	Instance inst = readInstance(filename);
	//- or define yourself by specifying: name, I, B, C, W, items (in format {c,w} per item)
	// Instance inst = createInstance("Test_instance", 10, 4, 3, 6, { {0,4},{0,3},{0,1},{1,3},{1,2},{1,2},{1,1},{2,3},{2,2},{2,1} });

	// Do pre-processing of the instance
	inst.preprocessing();

	// Find the minimum number of required bins
	// - either read the minimum number of required bins from a file
	string filenameBmin = "InstancesBPPMCF/minNumberOfBinsPerInstance.txt";
	inst.findMinNumberOfBins(filenameBmin);
	// - or find the minimum number of bins using the reflect model
	// inst.Bmin = solveReflectCSP(inst, inst.colorlessItems).UB;

	// Set the number of bins equal to the minimum number of required bins
	if (inst.Bmin > 0) { inst.B = inst.Bmin; } // change the number of bins

	// Find the L2-lower bounds
	inst.findL2();

	// Print information about the instance
	inst.print();

	// -------------------------------------------------------------------------

	//// Find the L^* bounds using Reflect
	//double start = getCPUTime();
	//inst.LBtot = 0;
	//for (int c = 0; c < inst.C; c++) {
	//	inst.LBs[c] = solveReflectCSP(inst, inst.groupedItems[c], inst.LBs[c], c).UB;
	//	inst.LBtot += inst.LBs[c];
	//}
	//double timeLB = start - getCPUTime();

	// Apply BPP-LB
	vector<int> LBs;
	Solution solBPPLB = solveBPPLB(inst, LBs);
	solBPPLB.print(true);

	// --> update the L2 bounds to the L* bounds
	inst.LBs = LBs;
	inst.LBtot = 0; for (int c = 0; c < inst.C; c++) { inst.LBtot += LBs[c]; }

	// Apply BPP-UB
	Solution solBPPUB = solveBPPUB(inst);
	solBPPUB.print(true);

	// Apply the tabu search metaheuristic TS starting from the feasible BPP-UB solution
	Solution solTS = solveTS(inst, solBPPUB);
	solTS.print(true);

	// --------------------------------------------------------------
}

/* Guide to setting up a project with Gurobi in Windows:
right - click on the project name in the solution explorer panel, then select properties.
set platform to x64
set configuration to all configurations
under debugging set environment to path = $(path); $(gurobi_home)\lib
under c / c++ > general > additional include directories, add: $(gurobi_home)\include
(under c / c++ > precompiled headers > precompiled header, select not using precompiled headers)
under linker > general > additional library directories, add: $(gurobi_home)\lib
set configuration to debug
under linker > input > additional dependencies, add gurobi91.lib; gurobi_c++mdd2019.lib
set configuration to release
under linker > input > additional dependencies, add gurobi91.lib; gurobi_c++md2019.lib
also set Solution platforms to x64 (top of Visual Studio, not inside properties)

copy-paste this:
path = $(path); $(gurobi_home)\lib; $(gurobi_home)\include; $(gurobi_home)\lib; gurobi91.lib;gurobi_c++mdd2019.lib; gurobi91.lib;gurobi_c++md2019.lib
OR for Gurobi 10: path = $(path); $(gurobi_home)\lib; $(gurobi_home)\include; $(gurobi_home)\lib; gurobi100.lib;gurobi_c++mdd2017.lib; gurobi100.lib;gurobi_c++md2017.lib
*/