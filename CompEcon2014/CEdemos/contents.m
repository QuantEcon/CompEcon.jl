%% Contents
%
% Demonstration programs to accompany
%   Lecture Notes in Computational Economic Dynamics
%   Mario J. Miranda, Instructor
%   The Ohio State University
%   
% Based on programs documented in book
%   Applied Computational Economics and Finance
%   Mario J. Miranda & Paul L. Fackler
%   2002, MIT Press, Cambridge MA
%
% These programs are continually revised
%   This version 4/9/2014
%
% Copyright 2014 Mario J. Miranda and Paul L Fackler

% Introduction 
demintro01 % Inverse Demand Problem
demintro02 % Rational Expectations Agricultural Market Model
 
% Mathematics Review 
demmath01  % Taylor Approximations
demmath02  % Computing Function Inner Products, Norms & Metrics
demmath03  % Continuous Distribution CDFs and PDFs
demmath04  % Standard Copulas
demmath05  % Illustrates Implicit Function Theorem
demmath06  % Simulate Simple Markov Chain
demmath07  % Operations with Markov Chains
 
% Linear Equations 
demlin01   % Solving Linear Equations by Different Methods
demlin02   % Ill-Conditioning of Vandermonde Matrix
demlin03   % Sparse Linear Equations
 
% Nonlinear Equations 
demslv01   % Compute Root of f(x)=exp(-x)-1
demslv02   % Compute Root of Rosencrantz Function
demslv03   % Compute Fixedpoint of f(x) = x^0.5
demslv04   % Compute Fixedpoint of f(x1,x2)= [x1^2+x2^3;x1*x2-0.5]
demslv05   % Cournot Equilibrium Model
demslv06   % Illustrates Nonlinear Equation Methods
demslv07   % Linear Complementarity Problem Methods
demslv08   % Nonlinear Complementarity Problem Methods
demslv09   % Hard Nonlinear Complementarity Problem with Billup's Function
demslv10   % Illustrates Linear Complementarity Problem
demslv11   % Nonlinear Complementarity Problem Methods
demslv12   % Convergence Rates for Nonlinear Equation Methods
demslv13   % Simple Nonlinear Complementarity Problem
demslv14   % Spatial Equilibrium Model
 
% Finite-Dimensional Optimization
demopt01   % Maximization via Golden Search
demopt02   % Changes in Nelder-Mead Simplex
demopt03   % Nelder-Mead Simplex Method
demopt04   % Maximization of Banana Function, Various Methods
demopt05   % Optimization with qnewton
demopt06   % KKT Conditions for Constrained Optimization Problems
demopt07   % Bound Constrained Optimization via Sequential LCP
demopt08   % Optimization with fmincon
 
% Quadrature
demqua01   % Equidistributed Sequences on Unit Square in R^2
demqua02   % Compute Expectation of Function of Random Normal Vector
demqua03   % Area under 1-D and 2-D Curves, Various Methods
demqua04   % Area under Normal PDF Using Simpson's Rule
demqua05   % Willingness to Pay, Expected Utility Model
demqua06   % Area under a Curve
demqua07   % Integration Using Trapezoidal Rule
demqua08   % Integration Using Simpson's Rule
demqua09   % Chebychev and Legendre Quadrature Nodes and Weights
demqua10   % Monte Carlo Simulation of Time Series
 
% Numerical Differentiation
demdif01   % Finite Difference Hessian Evaluation
demdif02   % Error in Finite Difference Differentiation
 
% Function Approximation
demapp00   % Approximation Using CompEcon Toolbox
demapp01   % Approximating Functions on R
demapp02   % Approximating Functions on R^2
demapp03   % Basis Functions and Standard Nodes for Major Approximation Schemes
demapp04   % Uniform-Node and Chebychev-Node Polynomial Approximation of Runge's Function
demapp05   % Chebychev Polynomial and Spline Approximation of Various Functions
demapp06   % Chebychev and Cubic Spline Derivative Approximation Errors
demapp07   % Solve Cournot Oligopoly Model via Collocation
demapp08   % Compute Function Inverse via Collocation
demapp09   % Linear Spline Approximation
demapp10   % Monopolist's Effective Supply Function
 
% Discrete Time Discrete State Dynamic Programming
demddp01   % Mine Management Model
demddp02   % Asset Replacement Model
demddp03   % Asset Replacement Model With Maintenance
demddp04   % Binomial American Put Option Model
demddp05   % Water Management Model
demddp06   % Bioeconomic Model
demddp07   % Renewable Resource Model
demddp08   % Job Search Model
demddp09   % Deterministic Cow Replacement Model
demddp10   % Stochastic Cow Replacement Model
demddp11   % Stochastic Optimal Growth Model

bigtic = tic;
 
% Discrete Time Continuous State Dynamic Programming
demdp00    % Timber Harvesting (Simple)
demdp01    % Timber Harvesting
demdp02    % Asset Replacement
demdp03    % Industry Entry-Exit
demdp04    % Job Search
demdp05    % American Put Option Pricing
demdp06    % Ramsey Deterministic Economic Growth
demdp06st  % Ramsey Stochastic Economic Growth
demdp07    % Stochastic Economic Growth
demdp08    % Public Renewable Management
demdp09    % Private Non-Renewable Resource Management
demdp10    % Water Resource Management
demdp11    % Monetary Policy
demdp12    % Production Management
demdp13    % Inventory Management
demdp14    % Livestock Feeding (Euler Conditions)
demdp15    % Saving with Transactions Costs
demdp16    % Linear-Quadratic Problem
demdp17    % Miscellaneous Lecture Note Figures
demdp19    % Credit with Strategic Default (Broyden)
demdp19alt % Credit with Strategic Default (dpsolve)
demdp20    % Lifecycle Consumption-Savings (Nonlinear Complementarity)
demdp20alt % Lifecycle Consumption-Savings (Backward Recursion)
demdp21new % Fertility-Savings Model (UNDER DEVELOPMENT)

toc(bigtic)

% Rational Expectations Models
demrem01   % Asset Pricing Model
demrem02   % Commodity Storage Model
demrem03   % Government Price Support Model
 
% Dynamic Games
demgame01  % Production Capacity Game
demgame02  % Income Redistribution Game
demgame03  % Marketing Board Game

% Ordinary Differential Equations
demode01   % Stability of Linear ODEs
demode02   % Generic Nonlinear ODE Example
demode03   % IVP Linear ODE Example
demode04   % Non-IVP Linear ODE Example
demode05   % Commodity Storage Model
demode06   % Predator-Prey Model
demode07   % Commercial Fisheries Model
demode10   % Lorentz Strange Attractor
demode11   % Tobin's Q
demode12   % Regional Migratio Model

% Continuous Time Deterministic Optimal Control
demode08   % Renewable Resource Model
demode09   % Optimal Growth Model

% Continuous Time Stochastic Dynamic Programming 
demcont01  % Ito Process
demcont02  % Renewable Resource Model (stochastic))
demcont03  % Optimal Growth Model (stochastic)
demcont04  % Fish Harvest Model (stochastic))
demcont05  % Production-Adjustment Model (stochastic))
demcont06  % Odd Nonrenewable Resource Model (deterministic)
demcont07  % Simple Log-Linear Example (stochastic)
demcont08  % Renewable Resource Model (deterministic)
 
% Financial Asset Pricing
demfin01   % Cox-Ingersoll-Ross Bond Pricing Model
demfin02   % Black-Scholes Option Pricing Model
demfin04   % American Put Option Pricing
demfin05   % Barrier Option Pricing
demfin06   % Compound Bond Option Pricing
demfin07   % Asian Option Pricing
demfin08   % Affine Asset Pricing
demfin09   % Financial Asset Calibration
 
% Regime Switching Problems
demrs01    % Asset Abandonment Model
demrs02    % Fish Harvest Model
demrs03    % Dixit Entry/Exit Model
 
% Impulse Control Problems
demic01    % Asset Replacement Model
demic02    % Timber Harvesting Model
demic03    % Storage Management Model
demic04    % Capacity Choice Model
demic05    % Cash Management Model
demic06    % Fish Harvest Model