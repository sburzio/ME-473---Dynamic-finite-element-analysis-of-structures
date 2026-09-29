The codes are divided for each material, here is the structure inside each material:

ENGINEERED WOOD:
- homework2_datagen_sweep.m: Matlab code to generate the data from the 2 variable chosen
- homework2_data_2.mat: Data manually saved from the datagen_sweep code
- homework2_ana.m: Analysis of the swept data -> creates the 4plot subplot
- homework2_ana_optimal.m: Analysis of the optimal bridge ->create all the other figures
- formStiffness2dtruss.m, formLumpedMass2Dtruss.m, formConsistentMass2DTruss.m, computeFrequenciesAndModes.m: are the helper function provided from the MATLAB Drive
- draw2Dtruss.m: is the modified version of the draw2Dtruss.m provided from the MATLAB Drive

STRUCTURAL STEEL: 
- Project2_wrapper_SS.m: overall wrapper file used to perform analysis. Open and press run to generate all figures and analysis used in the report. 
- Project2_calc_SS.m: code which performs actual FEM analysis. This code is called by the wrapper code to analyze each individual bridge design. Calculates displacements and natural frequencies.
- Howe.m: function which takes as an input anchor points for a truss, its maximum height and its complexity. It outputs the nodes and connectivity table for a Howe truss with these parameters.
- colormapplot.m: function which creates plots of truss structures with colors of stress in each member
- labelpoints.m: function to add labels to points on a plot
- noNumPlot.m: function which generates a simple plot of nodes and bars (no numbering). Inputs are nodes coordinates and connectivity table
- simpleplot.m: function which generates a plot of nodes and bars with the option to add numbering for nodes and/or bars. Inputs are nodes coordinates and connectivity table.

ALUMINUM ALLOY:
- Aluminum_alloy_1.m: this code allows the user to change the height h and shape factor of the bridge n and look at the displacement plot and selected mode shape (see code)
- Aluminum_alloy_2.m: this code is used to study the impact of h and n on the overall bridge performances. It consists of two for-loops sweeping from two ranges of h and n.

