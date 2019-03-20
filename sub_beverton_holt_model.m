
function [Abundance,Biomass,Yield_abundance,Yield_biomass,n_h,n_r,settlers,recruits] = sub_beverton_holt_model(PR,r,a,H,f,A,mass,V)

% The input variables are:
%    PR   - The proportion of habitat reserved
%    r    - the Beverton-Holt parameters
%    a    - Natural survival rate
%    H    - Fishing mortality outside the reserves
%    f    - Fecundity vector
%    A    - Number of age-classes
%    mass - Mass vector (for biomass yield)
%    V    - Larval mortality rate

% The output variables are:
%    Abundance          - The total abundance across the metapopulation
%    Biomass            - The total biomass across the metapopulation
%    Yield_abundance    - The total yield in raw individual numbers
%    Yield_biomass      - The total yield in biomass
%    n_h                - The density in harvested locations
%    n_r                - The density in reserved locations
%    settlers           - The density of settlers
%    recruits           - The density of recruits (e.g., the settlers, after density-dependent mortality has acted)

% Initialise the reserved (r) and harvested (h) populations
% ** Note that we're modelling the density, not the total populations **
n_r = ones(A,1)./A; n_h = n_r;

% Total initial abundance is the density in reserved area multiplied by the proportion of reserves, plus 
% the density in harvested area multiplied by the proportion of harvested habitat.
N = sum(n_r*PR + n_h*(1-PR));

% Beverton Holt parameters
alpha = r(1);
beta = r(2);

%% Forward simulate the population
delta = 1; t = 1;
while delta > 1e-6 & t < 1e5
   
   % The model dynamics follow the sequence defined in Hastings & Botsford (1999 Science)
   % - Reproduction
   % - Natural mortality
   % - Settlement
   % - Harvest
   
   % Larvae are produced and put into a pool (before mortality)
   settlers = V.*sum(f.*(n_r*PR + n_h*(1-PR))); % Larval Pool with Equal Redistribution
   
   % Recruitment mortality occurs according to a Beverton-Holt model. We haven't explicitly 
   % spread the settlers out among the habitat before we applied this, because that spreading 
   % parameter is assumed to also be inside "V".
   recruits = alpha.*settlers./(1 + beta.*settlers); 

   % Older age classes in all areas experience natural mortality (at rate 1-a)
   n_h = n_h.*a;
   n_r = n_r.*a;

   % Individuals get older
   n_r(2:end) = n_r(1:end-1);
   n_h(2:end) = n_h(1:end-1);
   
   % Recruitment occurs
   n_r(1) = recruits; 
   n_h(1) = recruits; 
   
   % The fished proportion of the habitat (1-PR) experiences fishing mortality at rate H
   Catch = (1-PR)*n_h.*H;
   n_h = n_h.*(1-H); 
   
   % Keep track of the overall population in the vector N
   N = [N sum(n_r*PR + n_h*(1-PR))];
   
   % Exit the loop if we've reached equilibrium (i.e., the total population has ceased to change)
   if t > 2e2; delta = abs(N(end) - N(end-1)); end; t = t+1;
end
Yield_biomass   = sum(Catch.*mass); % Fishing yield (biomass)
Yield_abundance = sum(Catch); % Fishing yield (numbers)
Abundance = N(end); % The most recent abundance

% Calculate total biomass across the system
Biomass = sum((n_r*PR + n_h*(1-PR)).*mass);













