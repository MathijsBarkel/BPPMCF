This repository contains the tested instances and code for all algorithms discussed in the paper "Bounds and Heuristic Algorithms for the Bin Packing Problem with Minimum Color Fragmentation" by Mathijs Barkel, Maxence Delorme, Enrico Malaguti and Michele Monaci. 
The code related to the paper "Pseudo-Polynomial Formulations for the Bin Packing Problem with Minimum Color Fragmentation" (by the same authors) can be found in repository "BPPMCF2".

All algorithms are coded in C++ and part of our methods require the commercial solver Gurobi (we used version 10.0.1). 
The files have the following contents:
- main.h/cpp                | example of front-end code for calling our methods to solve a given instance
- BPPLB.h/cpp               | code required for the algorithm BPP-LB
- BPPUB.h/cpp               | code required for the algorithms BPP-UB and BPP-UB*
- TS.h/cpp                  | code required for the algorithm TS
- AFCSP.h/cpp               | code required for the arcflow formulation (which is used by our algorithms)
- helper_functions.h/cpp    | code containing miscellaneous simple/supportive functions.

Moreover, "InstancesBPPMCF.zip" contains a txt-file for each of our test instances. These are spread over 5 folders, each corresponding to a different data set:
"Dataset 1", "Dataset 2", "Dataset 3" and "Dataset 4" correspond to datasets D1, D2, D3 and D4 as introduced by Mehrani et al. (2022) (see https://github.com/saharnazmehrani/BPPMCF-IJOC), and "Triplets" corresponds to D5*.
Our new instances D1*, D2*, D3* and D4* are obtained by setting the number of bins to the minimum number of required bins (see main.cpp),
where the minimum number of required bins per instance is saved in the file "minNumberOfBinsPerInstance.txt".

Each instance file is structured as follows:
For Datasets 1,2,3 and 4:
- The first line is always a 1
- The second line contains the number of bins (B)
- The third line contains the bin capacity (W)
- Then a blank line, followed by a 0-matrix of size BxW, followed by three more blank lines
- The next line contains the number of colors (C)
- The next line contains the number of items (I)
- Then follows a blank line
- Finally there is a matrix of size Ix2, where each row i gives the color c and the weight w of item i
For Dataset 5 (Triplets):
- The first line contains the number of items I
- The remaining lines contain a matrix of size Ix2, where each row i gives the color c and the weight w of item i.
- For these instances we always have B = I/3, W = 1000 and C = 3.

Finally, "Best LBs and UBs.xlsx" contains the best known lower bounds and upper bounds per instance.

For questions on the code, please send an e-mail to mathijsbarkel3[at]gmail[dot]com.

References:
- Barkel, M., Delorme, M., Malaguti, E., and Monaci, M. (2025). Bounds and heuristic algorithms for the bin packing problem with minimum color fragmentation. European Journal of Operational Research, 320(1):57–68.
- Mehrani, S., Cardonha, C., and Bergman, D. (2022). Models and algorithms for the bin-packing problem with minimum color fragmentation. INFORMS Journal on Computing, 34:1070–1085.
