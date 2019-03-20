function SteepnessSearch()

% The goal of this function is to estimate the value of two demographic parameters 
% that create a fishery with a particular steepness value. The parameters we're looking 
% for are:
%     * alpha, the Beverton-Holt parameter. We're going to assume that beta = 1 for simplicity and only worry about alpha
%     * V, an amplification for the fecundity parameter. It can be interpreted as the rate of larval mortality.

% What steepness values are we looking to recreate?
Steepness_goal_VEC = [0.3 0.5 0.7 0.9];

% These are upper and lower bounds for the search function
LB = [-30 -30]; UB = [5 3];

% load BestFitParameters Parameters
load Species_parameters_matlab Name NumSpp

for SP = 1:NumSpp
   disp(['Species ' num2str(SP)])
   disp(['Analysing ' Name{SP}])
   
   % Tailor the original guess for the FMINCON optimisation algorithm
   if     SP == 1;  X0 = [-5 -10]; 
   elseif SP == 9;  X0 = [1 -3]; 
   elseif SP == 11; X0 = [1 1];  
   elseif SP == 12; X0 = [1 1]; 
   elseif SP == 17; X0 = [-5 -3]; 
   elseif SP == 19; X0 = [-5 -16]; 
   elseif SP == 20; X0 = [1 1]; 
   elseif SP == 21; X0 = [1 1]; 
   elseif SP == 24; X0 = [1 1]; 
   elseif SP == 25; X0 = [1 1]; 
   elseif SP == 26; X0 = [1 1]; 
   elseif SP == 27; X0 = [1 -10]; 
   elseif SP == 29; X0 = [1 1]; 
   elseif SP == 30; X0 = [1 1]; 
   else X0 = [-1 -2];
   end
      
   for s = 1:4
      % What steepness are we looking for here?
      Steepness_goal = Steepness_goal_VEC(s);
      
      % Calculate the parameter values that most closely approximate this goal steepness
      P_opt = fmincon(@sub_Calculate_steepness,X0,[],[],[],[],LB,UB,[],[],Steepness_goal,SP,0);

      % What are the corresponding parameters?
      [STR,ST,SB,VB] = sub_Calculate_steepness(P_opt,Steepness_goal,SP,1);
      
      % Save them in a parameter vector
      Parameters(SP,:,s) = [Steepness_goal P_opt ST SB VB]; pause(0.01)
      
      % For computational reasons, initialise the search for the next steepness values 
      % with the optimal solution to the previous steepness values
      X0 = P_opt;
   end
   % Save the best fit values
   save BestFitParameters
end
save BestFitParameters

function [Steepness_residual,Steepness,SteepBiomass,VirginBiomass] = sub_Calculate_steepness(P,Steepness_goal,SP,PT)

% Load parameters for this particular species set using the script "Load_parameters"
Load_parameters

   % Choose the DD parameters to achieve certain levels of steepness
   % "Steepness" is defined as the fraction of virgin recruitment that is attained when
   % current spawning biomass is 20% of virgin spawning biomass

% For accuracy and readability, we're working with the logarithm of the best-fit parameters
r = [exp(P(1)) Beta]; % DD parameters
V = exp(P(2)); % Amplification for fecundity

% Calculate MSY for traditional fisheries
[HarvestRate,Biomass,Recruitment] = sub_TraditionalYield(A,f,a,r,mass,Length,V);

% Where is biomass 20% of virgin?
[BiomassResidual,I] = min(abs(Biomass - 0.2*Biomass(1)));

% Steepness definition
Steepness = Recruitment(I)/Recruitment(1);
Steepness_residual = abs(Steepness_goal - Steepness);% + abs(BiomassResidual);
SteepBiomass = Biomass(I);
VirginBiomass = Biomass(1);

if Biomass(1) < 1e-2
   Steepness_residual = inf;
end

% PLot the results so we can see that everything makes sense
if PT == 1
   
   if Steepness_goal == 0.3
      figure(1), clf; subplot(2,2,1), hold on; box on; FS = 10;
   elseif Steepness_goal == 0.5
      subplot(2,2,2), cla, hold on; box on; FS = 10;
   elseif Steepness_goal == 0.7
      subplot(2,2,3), cla, hold on; box on; FS = 10;
   elseif Steepness_goal == 0.9
      subplot(2,2,4), cla, hold on; box on; FS = 10;
   end
   
   disp(['Steepness = ' num2str(Steepness)])
   disp(['Best r parameter = ' num2str(P(1))])
   disp(['Best V parameter = ' num2str(P(2))])
   disp(['Steep biomass = ' num2str(Biomass(I),3)])
   disp(['Virgin biomass = ' num2str(Biomass(1),3)])
   title(['Sp = ' num2str(SP) '; Steep = ' num2str(Steepness_goal) '; V. biomass = ' num2str(round(Biomass(1)))],'fontweight','normal')
   plot(HarvestRate,Biomass./Biomass(1),'-','linewidth',2)
   plot(HarvestRate,Recruitment./Recruitment(1),'-','linewidth',2)
   plot([HarvestRate(I) HarvestRate(I)],[0 1],'k--','linewidth',1)
   plot([0 1],[Steepness_goal Steepness_goal],'k--','linewidth',1)
   text(0.4,Steepness_goal+0.05,Name{SP},'fontsize',FS)
   plot([0 1],[0.2 0.2],'k--','linewidth',1)
   xlabel('Additional mortality','fontsize',FS)
   ylabel('Percent virgin','fontsize',FS)
   L = legend('Biomass','Recruitment'); set(L,'fontsize',FS,'location','best')

   if Steepness_goal == 0.9
      Make_TIFF(['Figures/Fig_Sp_' num2str(SP)],[0 0 20 20],'-r200')
   end
end

%% ========================================================
%% ======= Calculate optimal traditional management =======
%% ========================================================
function [HVec,B,R] = sub_TraditionalYield(A,f,a,r,mass,Length,V)

% Go through a wide range of harvest rates
HVec = linspace(0,1,200);
for h = 1:length(HVec)
   H = HVec(h);
   
   % Run the model to equilibrium with this set of parameters, but with no MPAs
   [N(h),B(h),Ya(h),Yb(h),n_h,n_r,S,R(h)] = sub_beverton_holt_model(0,r,a,H,f,A,mass,V);

   % If the population is extinct, then any higher harvest rates will also drive it to extinction.
   % We therefore may as well stop the loop here, and assume the remaining values
   if N(h) == 0
      N(h+1:length(HVec)) = 0; B(h+1:length(HVec)) = 0; R(h+1:length(HVec)) = 0;
      break
   end
end




