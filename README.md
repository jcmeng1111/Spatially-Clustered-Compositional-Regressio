# Spatially Clustered Compositional Regression
Welcome to the GitHub repository for Spatially Clustered Compositional Regression: A Nonparametric Bayesian Approach. This work introduces a compositional regression model with spatially clustered coefficients to assess the varying importance of compositional predictors across spatial locations, all within a nonparametric Bayesian framework.
# Simulation Results
All simulation code from the paper is located in the Simulation folder. The simulations are categorized as follows:
Simulation 1: Cluster Setting 1 + Parameter Setting 1
Simulation 2: Cluster Setting 2 + Parameter Setting 1
Simulation 3: Cluster Setting 1 + Parameter Setting 2
Simulation 4: Cluster Setting 2 + Parameter Setting 2
Workflow:
Run the code in the Model_selection file to determine the optimal value of lambda.
Apply the selected lambda in the corresponding simulation code.
# Real data analysis
Real_d=3.R: Performs real data analysis without any spatial constant factor.
ReadData_d=3_eta.R: Performs real data analysis with the unemployment rate included as a spatial constant factor.


