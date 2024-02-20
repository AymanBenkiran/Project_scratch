#include <arrayfire.h>
#include <iostream>
#include <cmath>


//Reference (Solves problems from Section 5.3.3.6 in book)
//https://github.com/lbm-principles-practice/code/blob/master/chapter5/poiseuille_BB.m

//Use in PowerShell
//C:\GitHub\Single_phase_LBM_scratch\Project1\X64\debug >

int main() {
    af::setBackend(AF_BACKEND_CPU);
    //af::setBackend(AF_BACKEND_OPENCL);
    //af::setBackend(AF_BACKEND_CUDA);
    af::setDevice(0);
    af::info(); std::cout << std::endl;
    std::cout << "Hello, Visual Studio C++ Project!" << std::endl;

    double scale = 2.0;                       //Set simulation size
    double NX = 5.0 * scale;               //Channel length
    double NY = 5.0 * scale;               //Channel width
    double Nsteps = 1e4 * std::pow(scale, 2.0);  //Number of simulation time steps
    double tau = std::sqrt(3.0 / 16.0) + 0.5; //relaxation time(BGK model) ((tau = sqrt(3 / 16) + 0.5 gives exact solution)
    double omega = 1.0 / tau;
    double u_max = 0.1 / scale;             // maximum velocity
    double nu = (2.0 * tau - 1.0) / 6.0;       //kinematic shear viscosity
    double Re = NY * u_max / nu;         //Reynolds number; scaling parameter in simulation

    //Lattice parameters (Note : zero is last) 
    double NPOP = 9; //number of velocities
    
    // Weights
    af::array w = { 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0 };
    af_print(w);

    // Velocities, x components
    af::array cx = { 1.0, 0.0, -1.0,  0.0, 1.0, -1.0, -1.0,  1.0, 0.0 };

    // Velocities, y components
    af::array cy = { 0.0, 1.0, 0.0, -1.0, 1.0,  1.0, -1.0, -1.0, 0.0 };

    // Node locations
    af::array x = af::iota(NX, af::dtype::f64) - 0.5;   //Warning of loss data (to investigate)
    af::array y = af::iota(NY, af::dtype::f64) - 0.5;   //Warning of loss data (to investigate)

    //// Display the results
    //af_print(x);
    //af_print(y);

    // Pressure conditions
    double gradP      = 8.0 * nu * u_max / std::pow(NY,2.0);
    double rho_outlet = 1.0;
    double rho_inlet  = 3.0 * (NX - 1.0) * gradP + rho_outlet;

    // Analytical solution : Poiseuille velocity
    double ybottom = 0.0;
    double ytop = NY;
    af::array u_analy = -4.0 * u_max / (std::pow(NY, 2.0)) * (y - ybottom) * (y - ytop);

    // Initialize populations
        // Initialize feq array
    af::array feq = af::constant(0.0, NX, NY, NPOP, af::dtype::f64); //Warning of loss data (to investigate)

        // Set values for feq
    for (int k = 0; k < NPOP; ++k) {
        for (int i = 0; i < NX; ++i) {
            for (int j = 0; j < NY; ++j) {
                feq(i, j, k) = w(k); // assuming density equal one and zero velocity initial state
            }
        }
    }

        // Assign feq to f and fprop
    af::array f = feq;
    af::array fprop = feq;

        // Convergence parameters
    double tol = 1e-12; // tolerance to steady state convergence
    double teval = 100.0; // time step to evaluate convergence

    af::array u_old = af::constant(0.0, NX, NY, af::dtype::f64); //Warning of loss data (to investigate)
    af::array u;
    af::array v;

        // Initialize timer
    af::timer myTimer;
    af::sync();
    myTimer = af::timer::start();

    // Main algorithm
    
        //Macroscopic variables
    double rho;
    af::array indices_component1;
    af::array indices_component2;
    af::array selectedCols_component1;
    af::array selectedCols_component2;

    for (int t = 1; t < Nsteps; t++) {
        // Compute macroscopic quantities
        // Density
        rho = af::sum<double>(fprop, 3);

        //  Momentum components
            // u
        indices_component1 = { 1, 5, 8 };
        indices_component2 = { 3, 6, 7 };
        selectedCols_component1 = fprop(af::span, indices_component1);
        selectedCols_component2 = fprop(af::span, indices_component2);
        u = af::sum(selectedCols_component1, 2) - af::sum(selectedCols_component2, 2);

            // v
        indices_component1 = { 2, 5, 6 };
        indices_component2 = { 4, 7, 8 };
        selectedCols_component1 = fprop(af::span, indices_component1);
        selectedCols_component2 = fprop(af::span, indices_component2);
        v = af::sum(selectedCols_component1, 2) - af::sum(selectedCols_component2, 2);
        
        // Check convergence
        if (std::fmod(t, teval) == 1.0) {
            double conv = af::abs(af::mean(u)).scalar<double>() / af::abs(af::mean(u_old - 1.0)).scalar<double>();

            if (conv < tol) {
                break;
            }
            else {
                u_old = u;
            }
        }

        for (int k = 0; k < NPOP; k++) {

            // Compute equilibrium distribution (linear equilibrium with incompressible model)
            feq = af::tile(af::reorder(w, 1, 0), 1, 1, NPOP) * (rho + 3 * (u * af::tile(cx, 1, 1, NPOP) + v * af::tile(cy, 1, 1, NPOP)));
        }
        
        // Collision step
        f = (1 - omega) * fprop + omega * feq;

        // Inlet / Outlet BC : PBBC(w / i = 1 and i = NX outside layers)
        for (int k = 0; k < NPOP; ++k) {
            f(af::seq(1), af::span, af::seq(k, k)) = w(k) * (rho_inlet + 3 * (cx(k) * u(NX - 1, af::span) + cy(k) * v(NX - 1, af::span))) + (f(NX - 1, af::span, af::seq(k, k)) - feq(NX - 1, af::span, af::seq(k, k)));
            f(af::seq(NX), af::span, af::seq(k, k)) = w(k) * (rho_outlet + 3 * (cx(k) * u(2, af::span) + cy(k) * v(2, af::span))) + (f(2, af::span, af::seq(k, k)) - feq(2, af::span, af::seq(k, k)));
        }

        for (int k = 0; k < NPOP; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    // Streaming step (Periodic streaming of whole domain)
                    int newx = 1 + std::fmod((i - 1 + static_cast<int>(cx(k).scalar<float>()) + NX), NX);
                    int newy = 1 + std::fmod((j - 1 + static_cast<int>(cy(k).scalar<float>()) + NY), NY);
                    fprop(newx, newy, k) = f(i, j, k);
                }
            }
        }

        // Boundary condition (bounce-back)
            // Top wall (rest)
        fprop(af::span, NY, 3) = f(af::span, NY, 1);
        fprop(af::span, NY, 6) = f(af::span, NY, 5);
        fprop(af::span, NY, 7) = f(af::span, NY, 8);

            // Bottom wall (rest)
        fprop(af::span, 0, 4) = f(af::span, 0, 2);
        fprop(af::span, 0, 8) = f(af::span, 0, 7);
        fprop(af::span, 0, 5) = f(af::span, 0, 6);

    }
    af::sync(); //See if it should be here or in the for loop
    double nbsSecondElapsed = af::timer::stop(myTimer);

    printf("timing_manualWithError_exemple3() took %g seconds\n", nbsSecondElapsed);
    af_print(feq);

    // Compute error : L2 norm
// (note i = 1 and i = NX are virtual, not fluid, layers; thus not considered in error)
    af::array error = af::constant(0.0, NX);
    for (int i = 1; i < NX - 1; ++i) {
        error(i) = af::sqrt(af::sum((u(i, af::span) - u_analy)) * af::sum((u(i, af::span) - u_analy))) / af::sqrt(af::sum(u_analy * u_analy));
    }
    double L2 = 1.0 / NX * af::sum<double>(error);

    // Accuracy information
    printf(" ----- accuracy information -----\n");
    printf("        L2(u): %g\n", L2);

    return 0;
}

//To use Arrayfire
//https://arrayfire.org/docs/using_on_windows.htm#gsc.tab=0

//To run it in PowerShell
//1- cd C:\GitHub\Single_phase_LBM_scratch\Project1\x64\Debug
//2- .\Project1.exe

/*
% This code accompanies
% The Lattice Boltzmann Method : Principles and Practice
% T.Krüger, H.Kusumaatmaja, A.Kuzmin, O.Shardt, G.Silva, E.M.Viggen
% ISBN 978 - 3 - 319 - 44649 - 3 (Electronic)
% 978 - 3 - 319 - 44647 - 9 (Print)
% http ://www.springer.com/978-3-319-44647-9
    %
    %This code is provided under the MIT license.See LICENSE.txt.
    %
    %Author : Goncalo Silva
    %
    %Example matlab code for computing a Poiseuille flow with BB
    % Solves problems from Section 5.3.3.6 in book*/


/*
    clear all
    close all
    clc

    % simulation parameters
    scale = 1;% set simulation size
    NX = 5 * scale;% channel length
    NY = 5 * scale;% channel width
    NSTEPS = 1e4 * scale ^ 2;% number of simulation time steps
    tau = sqrt(3 / 16) + 0.5;% relaxation time(BGK model) (tau = sqrt(3 / 16) + 0.5 gives exact solution)
    omega = 1 / tau;
u_max = 0.1 / scale;% maximum velocity
nu = (2 * tau - 1) / 6;% kinematic shear viscosity
Re = NY * u_max / nu;% Reynolds number; scaling parameter in simulation


% Lattice parameters; note zero direction is last
NPOP = 9;% number of velocities
w = [1 / 9 1 / 9 1 / 9 1 / 9 1 / 36 1 / 36 1 / 36 1 / 36 4 / 9];% weights
cx = [1 0 - 1  0 1 - 1 - 1  1 0];% velocities, x components
cy = [0 1  0 - 1 1  1 - 1 - 1 0];% velocities, y components


% Node locations
x = (1:NX) - 0.5;
y = (1:NY) - 0.5;

% Pressure conditions
gradP = 8 * nu * u_max / NY ^ 2;
rho_outlet = 1;
rho_inlet = 3 * (NX - 1) * gradP + rho_outlet;

% Analytical solution : Poiseuille velocity
ybottom = 0;
ytop = NY;
u_analy = -4 * u_max / (NY ^ 2).*(y - ybottom).*(y - ytop);

% initialize populations
feq = zeros(NX, NY, NPOP);
for k = 1:NPOP
feq(:, : , k) = w(k);% assuming density equal one and zero velocity initial state
end
f = feq;
fprop = feq;

% convergence parameters
tol = 1e-12;% tolerance to steady state convergence
teval = 100;% time step to evaluate convergence
u_old = zeros(NX, NY);

% initalize clock
tstart = tic;

% Main algorithm
for t = 1:NSTEPS
% Compute macroscopic quantities
% density
rho = sum(fprop, 3);

% momentum components
u = sum(fprop(:, : , [1 5 8]), 3) - sum(fprop(:, : , [3 6 7]), 3);
v = sum(fprop(:, : , [2 5 6]), 3) - sum(fprop(:, : , [4 7 8]), 3);

% check convergence
if mod(t, teval) == 1

conv = abs(mean(u(:)) / mean(u_old(:)) - 1);

if conv < tol
    break
else
u_old = u;
end
end


for k = 1:NPOP
% Compute equilibrium distribution(linear equilibrium with incompressible model)
feq(:, : , k) = w(k) * (rho + 3 * (u * cx(k) + v * cy(k)));
end
% Collision step
f = (1 - omega) * fprop + omega * feq;

% Inlet / Outlet BC : PBBC(w / i = 1 and i = NX outside layers)
for k = 1 : NPOP
f(1, :, k) = w(k) * (rho_inlet + ...
    3 * (cx(k) * u(NX - 1, :) + cy(k) * v(NX - 1, :)))...
    + (f(NX - 1, :, k) - feq(NX - 1, :, k));
f(NX, :, k) = w(k) * (rho_outlet + ...
    3 * (cx(k) * u(2, :) + cy(k) * v(2, :)))...
    + (f(2, :, k) - feq(2, :, k));
end


for k = 1:NPOP
for j = 1 : NY
for i = 1 : NX

% Streaming step(Periodic streaming of whole domain)
newx = 1 + mod(i - 1 + cx(k) + NX, NX);
newy = 1 + mod(j - 1 + cy(k) + NY, NY);
fprop(newx, newy, k) = f(i, j, k);

end
end
end


% Boundary condition(bounce - back)
% Top wall(rest)
fprop(:, NY, 4) = f(:, NY, 2);
fprop(:, NY, 7) = f(:, NY, 5);
fprop(:, NY, 8) = f(:, NY, 6);

% Bottom wall(rest)
fprop(:, 1, 2) = f(:, 1, 4);
fprop(:, 1, 5) = f(:, 1, 7);
fprop(:, 1, 6) = f(:, 1, 8);
end


% Calculate performance information after the simulation is finished
runtime = toc(tstart);

% Compute error : L2 norm
% (note i = 1 and i = NX are virtual, not fluid, layers; thus not considered in error)
error = zeros(NX, 1);
for i = 2:NX - 1
error(i) = (sqrt(sum((u(i, :) - u_analy). ^ 2))). / sqrt(sum(u_analy. ^ 2));
end
L2 = 1 / NX * sum(error);

% Accuracy information
fprintf(' ----- accuracy information -----\n');
fprintf('        L2(u): %g\n', L2);
*/


