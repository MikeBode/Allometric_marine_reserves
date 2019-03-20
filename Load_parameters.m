
load Species_parameters_matlab

% Define the number of age classes, and a corresponding vector
A = 50;
ages = [1:A];

% These are von Bertalaffy parameters and equation, extracted 
% from the vectors for this particular species
L_inf = L_inf_vec(SP);
K = K_vec(SP);
z = z_vec(SP);
a0 = -0.25; % This is a negative constant, indicating that the species are not born with zero length
Length = L_inf.*(1 - exp(-K.*(ages - a0)));

% % Calculate natural mortality and survival based on Charnov et al. 
M_continuous = K*(Length./L_inf).^(-1.5);
a = exp(-M_continuous)';
NM = 1-a;

% Calculate mass from length according to the length-weight parameters
mass     = LW_a_vec(SP).*Length.^LW_b_vec(SP); mass = mass';
mass_inf = LW_a_vec(SP).*L_inf.^LW_b_vec(SP);

% Relative fecundity in each age class. It has a coefficient, but we can subsume it into the free variable of larval mortality
f = mass.^z; 

% Fecundity is zero before a mimumum length
First_repro = min(find(Length > Min_f_length_vec(SP)));
f(1:First_repro-1) = 0;

% Beverton-Holt parameter beta
Beta = 1;

