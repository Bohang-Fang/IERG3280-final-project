# IERG6130-final-project
Report figure.zip consolidates all the figures in my report. The SIR_simulation code provides an example of the dynamic process of the SIR model (ODE equations) and incorporates immunity loss, including the SIRS model, along with code that reflects my proposed assumption. In this model, beta represents the infection rate, gamma represents the recovery rate, delta corresponds to the rate of immunity loss, N denotes the total population, and I0 represents the initial state of infected individuals. The time period considered is 300 days, and these parameters can be modified to generate different plots.

The SI, SIS, SEIR, and SEIRS models are all based on ODEs, and the simulation method is similar to the SIR model.

Cstar.m, figure3ks.m, and lamda_s.m reproduce the experimental results presented in paper [5] in my report. figure3ks.m illustrates the results depicted in figure 11 of my report, while lamda_s.m displays the results shown in figures 9 and 10. Before running the codes in Cstar.m, figure3ks.m, and lamda_s.m, you need to install the CVX packages in Matlab.
