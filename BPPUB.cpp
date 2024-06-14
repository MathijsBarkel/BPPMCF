#include "BPPUB.h"

Solution solveBPPUB(const Instance& inst, bool stopEarly, bool exactRecoloring, double timeLimit) {
	// This function applies the algorithm BPP-UB to the given instance

	double start = getCPUTime();   // starting time
	Solution sol;				   // initialize Solution struct
	sol.method = "BPP-UB";		   // method declaration
	if (not stopEarly) { sol.method += "(NoEarlyStop)"; }
	if (exactRecoloring) { sol.method += "(ExactRecoloring)"; }

	// set a random seed
	srand(53);                   // option 1: every time the same
	// srand((unsigned)time(NULL));	// option 2: every time different

	// set the number of tries that the second phase is attempted
	int nTries = 100;

	vector<vector<int>>	wq = inst.colorlessItems; // combine all items (ignoring color)

	// Create graph
	Graph G = graphConstructionCSP(wq, inst.W);
	// G.print(); // print the graph

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();							// remove Gurobi message
	env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	vector<GRBVar> f(G.A.size());

	// declaration of the variables for the model
	// for a in A: fa = # arc a is used
	int idx;
	for (int a = 0; a < G.A.size(); a++) {
		idx = G.A[a][2];
		if (idx != -1) {
			f[a] = model.addVar(0, wq[idx][1], 0, GRB_INTEGER);
		}
		else {
			f[a] = model.addVar(0, inst.B, 0, GRB_INTEGER);
		}
	}
	model.update();

	// declare linear expressions
	vector<GRBLinExpr> fIn(inst.W + 1, 0);    		// the amount of flow entering each vertex
	vector<GRBLinExpr> fOut(inst.W + 1, 0);   		// the amount of flow leaving each vertex
	vector<GRBLinExpr> typeUsed(wq.size(), 0); 	// the amount of arcs used of each item type 

	// calculate the linear expressions	
	for (int a = 0; a < G.A.size(); a++) {   				// loop over all arcs
		fIn[G.A[a][1]] += f[a];        						// inflow
		fOut[G.A[a][0]] += f[a];       						// outflow
		if (G.A[a][2] >= 0) typeUsed[G.A[a][2]] += f[a];	// number of items used of certain type
	}
	model.update();

	// set the objective: minimize the number of flow going out of 0
	model.setObjective(fOut[0], GRB_MINIMIZE);
	// constraints 1: flow conservation
	for (int v = 1; v < inst.W; v++) {       	// loop over all vertices
		if (G.V[v])
			model.addConstr(fIn[v] == fOut[v]); // inflow = outflow
	}
	model.addConstr(fOut[0] == fIn[inst.W]);	// total flow
	// constraints 2: item type quantities
	for (int j = 0; j < wq.size(); j++) {			// loop over all item types
		model.addConstr(typeUsed[j] == wq[j][1]);	// demand met
	}

	model.update();
	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// (optionally:) stop early if the required number of bins is reached
	if (stopEarly) { model.getEnv().set(GRB_DoubleParam_BestObjStop, inst.B); }

	// find the optimal solution
	model.optimize();

	// store the results in a Solution object
	sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
	sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
	sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients	
	sol.feas = -1;
	sol.opt = false;
	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);

	// if the instance is infeasible
	if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		sol.feas = 0;
		sol.opt = false;
	}
	// if a solution has been found within the time limit
	else if (model.get(GRB_IntAttr_SolCount) >= 1) {
		int nBinsUsed = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		// if the LB is more than the amount of bins --> infeasible instance
		if (sol.LB > inst.B) {
			sol.feas = 0;
			sol.opt = false;
		}
		// if time ran out and the LB is at most the amount of bins, but the UB is more, we don't know if the instance is feasible
		else if (nBinsUsed > inst.B) {
			sol.feas = -1;
			sol.opt = false;
		}
		// if the solution is optimal for the CSP instance and the amount of bin used is at most the amount of bins available
		else if (nBinsUsed <= inst.B) {
			// decompose the solution in a smart way to find the amount of color fragmentation

			if (not exactRecoloring) { // heuristic recoloring through greedy flow decomposition
				sol.feas = 1;
				// first store the used arcs by tail
				int fval;
				vector<vector<vector<int>>> AByTailOriginal(inst.W + 1);
				for (int a = 0; a < G.A.size(); a++) {
					fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
					if (fval > 0) {
						AByTailOriginal[G.A[a][0]].push_back({ G.A[a][1], G.A[a][2], fval }); // add arc to set of arcs with same tail
					}
				}

				// for every type of arc we want the available colors and their quantities
				// (as an intermediate step: every arc type belongs to an item size, which belongs to some colors)
				// colorsPerSize_wc gives the number of items of size w and color c
				vector<vector<int>> colorsPerSizeOriginal(inst.W + 1, vector<int>(inst.C, 0));
				// inst.countedItems contains for every item type the color, size and quantity (c,w,q)
				for (int j = 0; j < inst.countedItems.size(); j++) {
					colorsPerSizeOriginal[inst.countedItems[j][1]][inst.countedItems[j][0]] = inst.countedItems[j][2];
				}
				// colorsPerType_jc gives the number of items of index j and color c
				vector<vector<int>> colorsPerTypeOriginal(wq.size());
				for (int j = 0; j < wq.size(); j++) {
					colorsPerTypeOriginal[j] = colorsPerSizeOriginal[wq[j][0]];
				}

				// declare variables
				vector<vector<vector<int>>> AByTail;
				vector<vector<int>> colorsPerType;
				vector<vector<vector<vector<int>>>> currentPacking, bestPacking;
				vector<vector<vector<int>>> currentBin;
				vector<vector<int>> goodCombinations, badCombinations;
				vector<int> a;
				int tail, head, idx, loss, chosenComb, chosenColor, chosenArc;
				vector<bool> usedColors;
				vector<vector<int>> itemsInBinPerColor;
				int curVal, bestVal = int('inf');

				// try decomposing the flow nTries times
				for (int t = 0; t < nTries; t++) {
					// reset the solution
					currentPacking.clear();
					curVal = 0;
					AByTail = AByTailOriginal;
					colorsPerType = colorsPerTypeOriginal;

					// decompose the flow one-by-one
					for (int b = 0; b < nBinsUsed; b++) { // loop over the to-be-filled bins
						usedColors.clear(); usedColors.resize(inst.C, false);
						itemsInBinPerColor.clear(); itemsInBinPerColor.resize(inst.C);
						currentBin.clear(); currentBin.resize(inst.C);
						tail = 0;				 // start a new path at 0
						while (tail != inst.W) { // while the terminanal vertex hasn't been reached
							loss = -1;			 // contains index of loss arc from current tail, or -1 if no such arc exists
							goodCombinations.clear(); badCombinations.clear();
							for (int a_idx = 0; a_idx < AByTail[tail].size(); a_idx++) { // check each possible outgoing arc
								idx = AByTail[tail][a_idx][1];			 // find the index of the arc
								if (idx != -1) {						 // if the arc is not a loss arc
									for (int c = 0; c < inst.C; c++) {   // check the possible colors
										if (colorsPerType[idx][c] > 0) { // if the color is available
											if (usedColors[c]) {		 // if the color was already used
												goodCombinations.push_back({ a_idx, c }); // add to good combinations
											}
											else {						 // if the color wasn't used yet 
												badCombinations.push_back({ a_idx, c }); // add to bad combinations
											}
										}
									}
								}
								else {				// if the arc is a loss arc
									loss = a_idx;	// save the index of the arc
								}
							}

							// if a good move is possible (i.e., using an arc of an already used color)
							if (goodCombinations.size() > 0) {
								chosenComb = rand() % goodCombinations.size();
								chosenArc = goodCombinations[chosenComb][0];
								chosenColor = goodCombinations[chosenComb][1];
							}
							// otherwise, use a loss arc if possible
							else if (loss != -1) {
								chosenArc = loss;
								chosenColor = -1;
							}
							// otherwise, make a bad move (i.e., using an arc of a new color)
							else {
								chosenComb = rand() % badCombinations.size();
								chosenArc = badCombinations[chosenComb][0];
								chosenColor = badCombinations[chosenComb][1];
							}

							// take the chosen arc + remove afterwards
							head = AByTail[tail][chosenArc][0];
							idx = AByTail[tail][chosenArc][1];
							AByTail[tail][chosenArc][2]--;
							if (AByTail[tail][chosenArc][2] == 0) {
								AByTail[tail].erase(AByTail[tail].begin() + chosenArc);
							}
							tail = head;	// move to the new point

							if (chosenColor != -1) {
								usedColors[chosenColor] = true;    // update list containing colors used in path
								colorsPerType[idx][chosenColor]--; // remove one arc of that size and color
								itemsInBinPerColor[chosenColor].push_back(wq[idx][0]); // add the item
							}

						}

						for (int c = 0; c < inst.C; c++) {			// for every color
							if (itemsInBinPerColor[c].size() > 0) { // if the color is contained in the bin
								curVal++;							// add 1 to the objective value
								currentBin[c] = sortAndCount(itemsInBinPerColor[c]); // group and count the items of the current color in the bin
							}
						}
						currentPacking.push_back(currentBin);		// add the last bin to the set of all bins
					}

					// update the best solution
					if (curVal < bestVal) {
						bestVal = curVal;
						bestPacking = currentPacking;
						// terminate if the solution is proven optimal
						if (bestVal == inst.LBtot) { break; }
					}
				}
				// save the best solution found
				sol.binPacking = bestPacking;
				sol.UB = bestVal;
				sol.LB = inst.LBtot;				// the lower bound is simply set to the L2 bound
				if (sol.UB == sol.LB) { sol.opt = true; }
			}

			else { // optimal recoloring 
				
				// Find the solution graph
				Graph GSol;
				GSol.V.resize(inst.W + 1, false); // initialize graph
				int fval, tail, head;
				for (int a = 0; a < G.A.size(); a++) {
					fval = ceil(f[a].get(GRB_DoubleAttr_X) - EPSILON); // find variable value
					if (fval > 0) {
						tail = G.A[a][0]; head = G.A[a][1];
						GSol.V[tail] = true; GSol.V[head] = true;
						GSol.A.push_back({ tail, head, fval });
					}
				}
				Solution sol2 = recolorGraph(inst, GSol, timeLimit + start - getCPUTime());
				sol2.print(false);
				sol.UB = sol2.UB;
				sol.LB = inst.LBtot;
				sol.LPrel = sol2.LB;    // This is not really the LPrel, instead it is the lower bound found by the graph recoloring problem
				sol.timeP += sol2.timeP;
				sol.binPacking = sol2.binPacking;
				sol.feas = sol2.feas;
				if (sol.UB == sol.LB) { sol.opt = true; }
			}
		}
	}

	sol.timeT = getCPUTime() - start;	// save the total time
	return sol;							// return the solution
}

Solution recolorGraph(const Instance& inst, Graph& G, double timeLimit) {

	double start = getCPUTime();		// starting time
	Solution sol;						// initialize Solution struct
	sol.method = "graphRecoloring";		// method declaration

	// G.print(); // print the graph

	// store the arcs by tail, count the number of flows and save the starting arcs
	vector<vector<int>> arcsByTail(inst.W);
	int nFlows = 0;
	vector<int> startingArcs;
	int tail;
	for (int a = 0; a < G.A.size(); a++) {
		tail = G.A[a][0];
		arcsByTail[tail].push_back(a);
		if (tail == 0) { 
			nFlows += G.A[a][2]; 
			for (int k = 0; k < G.A[a][2]; k++) {
				startingArcs.push_back(a);
			}
		}
	}
	// print2DVector(arcsByTail, "arcsByTail", "tail");

	// find the possible colors per size
	vector<vector<int>> colorsPerSize(inst.W + 1);
	for (int j = 0; j < inst.J; j++) {
		colorsPerSize[inst.countedItems[j][1]].push_back(inst.countedItems[j][0]);
	}
	// print2DVector(colorsPerSize, "colorsPerSize", "size");

	// find a mapping from an items color and size to its index
	vector<vector<int>> cw2j(inst.C, vector<int>(inst.W + 1, -1));
	int c, w;
	for (int j = 0; j < inst.J; j++) {
		c = inst.countedItems[j][0]; w = inst.countedItems[j][1];
		cw2j[c][w] = j;
	}

	// create the subgraph for every bin
	vector<Graph> Gsubs(nFlows);
	int a, head;
	for (int b = 0; b < nFlows; b++) {

		// initialize the graph
		Gsubs[b].V.resize(inst.W + 1, false); Gsubs[b].V[0] = true;

		// add the first arc
		a = startingArcs[b];
		tail = G.A[a][0]; head = G.A[a][1]; w = head - tail;
		for (const int& c: colorsPerSize[w]) {
			Gsubs[b].A.push_back(G.A[a]);
			Gsubs[b].A.back()[2] = c;
		}
		Gsubs[b].V[head] = true;

		// add the other arcs
		for (int v = 1; v < inst.W; v++) {									// loop over the vertices
			if (Gsubs[b].V[v]) {											// if the vertex exists
				for (const int& a : arcsByTail[v]) {						// loop over the arc indices leaving the vertex
					tail = G.A[a][0]; head = G.A[a][1]; w = head - tail;	// decompose the arc
					for (const int& c : colorsPerSize[w]) {					// loop over all possible arc colors
						Gsubs[b].A.push_back(G.A[a]);
						Gsubs[b].A.back()[2] = c;
					}
					Gsubs[b].V[head] = true;
				}
			}
		}
	}

	// create a model
	GRBEnv env = GRBEnv();              	// create an environment
	removeLine();                       	// remove Gurobi message
	// env.set(GRB_IntParam_LogToConsole, 0);  // turn off console output
	env.set(GRB_IntParam_Seed, 53);			// set a random seed to be used by Gurobi
	GRBModel model = GRBModel(env);         // create a new model

	// declaration of the variables for the model
	// f[b][a] = 1 iff bin b contains arc a, 0 otherwise
	// z[b][c] = 1 iff bin b contains color c, 0 otherwise
	vector<vector<GRBVar>> f(nFlows);
	vector<vector<GRBVar>> z(nFlows, vector<GRBVar>(inst.C));
	for (int b = 0; b < nFlows; b++) {						// loop over all bins
		f[b].resize(Gsubs[b].A.size());
		for (int a = 0; a < Gsubs[b].A.size(); a++) {		// loop over all arcs
			f[b][a] = model.addVar(0, 1, 0, GRB_BINARY);
		}
		for (int c = 0; c < inst.C; c++) {					// loop over all colors
			z[b][c] = model.addVar(0, 1, 0, GRB_BINARY);
		}
	}

	model.update();

	// declare linear expressions
	vector<vector<GRBLinExpr>> fIn(nFlows, vector<GRBLinExpr>(inst.W + 1, 0));    		// the amount of flow entering each vertex for each bin
	vector<vector<GRBLinExpr>> fOut(nFlows, vector<GRBLinExpr>(inst.W + 1, 0));   		// the amount of flow leaving each vertex for each bin
	vector<GRBLinExpr> typeUsed(inst.J, 0); 											// the number of arcs used of each item type 
	vector<vector<vector<GRBLinExpr>>> crossingArcsPerColor(nFlows, vector<vector<GRBLinExpr>>(inst.C, vector<GRBLinExpr>(inst.W, 0)));	// the number of arcs in bin b of color c crossing vertex v
	vector<GRBLinExpr> binsPerColor(inst.C, 0);											// the number of bins containing each color
	GRBLinExpr colorFragmentations = 0;													// the total number of color fragmentations

	// calculate the linear expressions	
	int j;
	for (int b = 0; b < nFlows; b++) {								// loop over all bins
		for (int a = 0; a < Gsubs[b].A.size(); a++) {   			// loop over all arcs
			tail = Gsubs[b].A[a][0]; head = Gsubs[b].A[a][1];		// decompose the arc
			c = Gsubs[b].A[a][2]; w = head - tail; j = cw2j[c][w];	// find more arc info
			fIn[b][head] += f[b][a];        						// inflow
			fOut[b][tail] += f[b][a];       						// outflow
			typeUsed[j] += f[b][a];									// number of items used of certain type
			for (int v = tail; v < head; v++) {						// for all vertices from the tail to (excl.) the head
				if (Gsubs[b].V[v]) {								// only if the vertex exists
					crossingArcsPerColor[b][c][v] += f[b][a];		// the number of arcs in bin b of color c crossing vertex v
				}
			}
			
		}
		for (int c = 0; c < inst.C; c++) {							// loop over all colors
			binsPerColor[c] += z[b][c];								// number of bins containing each color
			colorFragmentations += z[b][c];							// the total number of color fragmentations
		}
	}
	model.update();

	// set the objective: minimize the total number of color fragmentations
	model.setObjective(colorFragmentations, GRB_MINIMIZE);

	// constraints 1: flow conservation
	for (int b = 0; b < nFlows; b++) {						// loop over all bins
		model.addConstr(fOut[b][0] == 1);					// exactly one unit of outflow from 0
		for (int v = 1; v < inst.W; v++) {       			// loop over all vertices
			if (Gsubs[b].V[v]) {							// if the vertex exists
				model.addConstr(fIn[b][v] >= fOut[b][v]);	// inflow >= outflow
			}
		}
	}

	// constraints 2: item type quantities
	for (int j = 0; j < inst.J; j++) {								// loop over all item types
		model.addConstr(typeUsed[j] == inst.countedItems[j][2]);	// demand met
	}
	
	//// constraints 3: behavior z (weaker version)
	//for (int b = 0; b < nFlows; b++) {					// loop over all bins
	//	for (int a = 0; a < Gsubs[b].A.size(); a++) {	// loop over all arcs
	//		c = Gsubs[b].A[a][2];						// find the color of the arc
	//		model.addConstr(f[b][a] <= z[b][c]);		// arc only allowed if color is allowed
	//	}
	//}

	// constraints 3: behavior z (stronger version)
	for (int b = 0; b < nFlows; b++) {											// loop over all bins
		for (int c = 0; c < inst.C; c++) {										// loop over all colors
			for (int v = 0; v < inst.W; v++) {									// loop over all vertices
				if (crossingArcsPerColor[b][c][v].size() > 0) {					// if there exists at least one arc in graph b with color c crossing vertex v
					model.addConstr(crossingArcsPerColor[b][c][v] <= z[b][c]);	// at most one of those arcs is allowed, but only if the corresponding z-variable is 1
				}
			}
		}
	}

	// constraints 4 (optional): L2-based constraints
	for (int c = 0; c < inst.C; c++) {
		if (inst.LBs[c] > 0) model.addConstr(binsPerColor[c] >= inst.LBs[c]);
	}

	model.update();
	sol.timeP = getCPUTime() - start;     // save the preprocessing time

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_IntParam_Method, 2); // use Interior Point Methods
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(timeLimit - sol.timeP, double(1)));

	// find the optimal solution
	model.optimize();

	// store the results in a Solution object
	sol.Nvar = model.get(GRB_IntAttr_NumVars);       // number of variables
	sol.Nconstr = model.get(GRB_IntAttr_NumConstrs); // number of constraints
	sol.Ncoeff = model.get(GRB_IntAttr_NumNZs);      // number of non-zero coefficients	
	sol.feas = -1;
	sol.opt = 0;
	sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
	sol.UB = 1000000;

	// if the instance is infeasible
	if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		sol.feas = 0;
	}

	// if a solution has been found within the time limit
	else if (model.get(GRB_IntAttr_SolCount) >= 1) {
		sol.feas = 1;
		sol.UB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		// if the solution is optimal
		if (sol.LB == sol.UB) {
			sol.opt = 1;
		}

		// store the items assigned to each bin (sorted by color)
		vector<vector<int>> binItemSizesPerColor;
		vector<vector<vector<int>>> currentBin;
		for (int b = 0; b < nFlows; b++) { // loop over the to-be-filled bins
			binItemSizesPerColor.resize(inst.C);	// the bin contains a seperate vector per color

			for (int a = 0; a < Gsubs[b].A.size(); a++) {	// loop over all arcs
				if (ceil(f[b][a].get(GRB_DoubleAttr_X) - EPSILON) == 1) { // if the arc is used
					w = Gsubs[b].A[a][1] - Gsubs[b].A[a][0];	// find the item size of the arc
					c = Gsubs[b].A[a][2];						// find the color of the arc
					binItemSizesPerColor[c].push_back(w);		// add the item to the bin
				}
			}
			currentBin.resize(inst.C);
			for (int c = 0; c < inst.C; c++) {
				if (binItemSizesPerColor[c].size() > 0) {
					currentBin[c] = sortAndCount(binItemSizesPerColor[c]); // group and count the items in the bin
				}
			}
			sol.binPacking.push_back(currentBin);		// add the last bin to the set of all bins
			binItemSizesPerColor.clear();				// start a new empty bin
			currentBin.clear();							// start a new empty bin
		}
	}
	sol.timeT = getCPUTime() - start; // save the total time

	return sol;                                 // return the solution
}
