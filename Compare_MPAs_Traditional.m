function Compare_MPAs_Traditional()

% First initialise the CSV file that's going to keep all the results
load Species_parameters_matlab NumSpp
for Steepness = 1:4
   CELL{1,1 } = 'Species';
   CELL{1,2 } = 'MSYb - Effort';
   CELL{1,3 } = 'MSYb - Mixed';
   CELL{1,4 } = 'MSYb - Scorched';
   CELL{1,5 } = 'MSYb - (Mixed/Effort)%';
   CELL{1,6 } = 'Av. length - Effort';
   CELL{1,7 } = 'Av. length - Mixed (reserved sites)';
   CELL{1,8 } = 'Av. length - Mixed (harvested sites)';
   CELL{1,9 } = 'Av. length - Scorched (reserved sites)';
   CELL{1,10} = 'SSB - Effort';
   CELL{1,11} = 'SSB - Mixed';
   CELL{1,12} = 'SSB - Scorched';
   CELL{1,13} = 'SS abundance - Effort';
   CELL{1,14} = 'SS abundance - Mixed';
   CELL{1,15} = 'SS abundance - Scorched';
   CELL{1,16} = 'Optimal mixed harvest rate';
   CELL{1,17} = 'Optimal mixed reserve proportion';
   CELL{1,18} = 'Abundance at biomass MSY - Effort';
   CELL{1,19} = 'Abundance at biomass MSY - Mixed';
   CELL{1,20} = 'Abundance at biomass MSY - Scorched';
   
   for SP = 1:NumSpp
      CELL = sub_Yield_modelling(SP,Steepness,CELL)
      
      if Steepness == 1
         Make_TIFF(['Figures/OptimalManagement_Sp_' num2str(SP) '_St_03'],[0 0 35 10],'-r100')
      elseif Steepness == 2
         Make_TIFF(['Figures/OptimalManagement_Sp_' num2str(SP) '_St_05'],[0 0 35 10],'-r100')
      elseif Steepness == 3
         Make_TIFF(['Figures/OptimalManagement_Sp_' num2str(SP) '_St_07'],[0 0 35 10],'-r100')
      elseif Steepness == 4
         Make_TIFF(['Figures/OptimalManagement_Sp_' num2str(SP) '_St_09'],[0 0 35 10],'-r100')
      end
      
   end
   
   if Steepness == 1
      cell2csv('AllSpeciesOutputs_St_03.csv',CELL)
   elseif Steepness == 2
      cell2csv('AllSpeciesOutputs_St_05.csv',CELL)
   elseif Steepness == 3
      cell2csv('AllSpeciesOutputs_St_07.csv',CELL)
   elseif Steepness == 4
      cell2csv('AllSpeciesOutputs_St_09.csv',CELL)
   end
   clear CELL
end

%% This is the function that actually runs the model
function CELL = sub_Yield_modelling(SP,Steepness,CELL)

if nargin == 0
   % Which species are we modelling?
   SP = 1;
   % How steep do we want the ecology? Options are [0.3 0.5 0.7 0.9], which we call 1,2,3,4.
   Steepness = 2;
end

% Load the right parameters for this species
% ---- Parameters(SP,:,s) = [Steepness_goal r_opt ST SB VB];
Load_parameters
load BestFitParameters Parameters
R1 = Parameters(SP,2,Steepness);
V = exp(Parameters(SP,3,Steepness));
r = [exp(R1) 1]; % DD parameters

figure(1); clf

% MPA only
H_outside = 1;
[MSY_R,MSYa_R,LFD_R_r,LFD_R_h,SSB_R,SSA_R,ORP_R] = sub_ReserveYield(H_outside,A,f,a,r,mass,Length,V,1);

% Mixed MPA-Effort
H_mixed_opt = fminbnd(@sub_ReserveYield,0,1,[],A,f,a,r,mass,Length,V,0);
[MSY_M,MSYa_M,LFD_M_r,LFD_M_h,SSB_M,SSA_M,ORP_M] = sub_ReserveYield(H_mixed_opt,A,f,a,r,mass,Length,V,1);

% Effort
[MSY_E,MSYa_E,LFD_E,SSB_E,SSA_E,Opt_Effort] = sub_TraditionalYield(A,f,a,r,mass,Length,V);

% Calculate the average lengths
LFD_E = LFD_E./sum(LFD_E);
LFD_M_h = LFD_M_h./sum(LFD_M_h);
LFD_M_r = LFD_M_r./sum(LFD_M_r);
LFD_R = LFD_R_r./sum(LFD_R_r);
AvLength_E = sum(LFD_E.*Length');
AvLength_M_h = sum(LFD_M_h.*Length');
AvLength_M_r = sum(LFD_M_r.*Length');
AvLength_R = sum(LFD_R_r.*Length');

disp(['Av length = ' num2str(AvLength_E,2) ' (effort)'])
disp(['Av length = ' num2str(AvLength_M_r,2) ' (mixed - protected)'])
disp(['Av length = ' num2str(AvLength_M_h,2) ' (mixed - harvested)'])
disp(['Av length = ' num2str(AvLength_R,2) ' (reserves only)'])
disp(['MSY       = ' num2str(MSY_E,4) ' (effort)'])
disp(['MSY       = ' num2str(MSY_M,4) ' (mixed)'])
disp(['MSY       = ' num2str(MSY_R,4) ' (reserves only)'])
disp(['%         = ' num2str(MSY_R./MSY_E,4) ' (reserves vs. effort)'])
disp(['%         = ' num2str(MSY_M./MSY_E,4) ' (mixed vs. effort)'])

CELL{SP+1,1 } = Name{SP};                             %'Species';
CELL{SP+1,2 } = MSY_E;                                %'MSYb - Effort';
CELL{SP+1,3 } = MSY_M;                                %'MSYb - Mixed';
CELL{SP+1,4 } = MSY_R;                                %'MSYb - Scorched';
CELL{SP+1,5 } = 100*(MSY_M/MSY_E);                    %'MSYb - Relative';
CELL{SP+1,6 } = AvLength_E;                           %'Av. length - Effort';
CELL{SP+1,7 } = AvLength_M_r;                         %'Av. length - Mixed (reserved sites)';
CELL{SP+1,8 } = AvLength_M_h;                         %'Av. length - Mixed (harvested sites)';
CELL{SP+1,9 } = AvLength_R;                           %'Av. length - Scorched (reserved sites)';
CELL{SP+1,10 } = SSB_E;                               %'SSB - Effort';
CELL{SP+1,11} = SSB_M;                                %'SSB - Mixed';
CELL{SP+1,12} = SSB_R;                                %'SSB - Scorched';
CELL{SP+1,13} = SSA_E;                                %'SS abundance - Effort';
CELL{SP+1,14} = SSA_M;                                %'SS abundance - Mixed';
CELL{SP+1,15} = SSA_R;                                %'SS abundance - Scorched';
CELL{SP+1,16} = H_mixed_opt;                          %'Optimal mixed harvest rate';
CELL{SP+1,17} = ORP_M;                                %'Optimal mixed reserve proportion';
CELL{SP+1,18} = MSYa_E;                               %'Abundance at biomass MSY - Effort'; 
CELL{SP+1,19} = MSYa_M;                               %'Abundance at biomass MSY - Mixed'; 
CELL{SP+1,20} = MSYa_R;                               %'Abundance at biomass MSY - Scorched'; 


%% ========================================================
%% ======= Calculate optimal MPA-effort management ========
%% ========================================================
function [MSY,MSY_a,LF_r,LF_h,SSB,SSA,Opt_reserve_proportion] = sub_ReserveYield(H_outside,A,f,a,r,mass,Length,V,PT)

% Go through all possible reserve sizes
ReserveSizes = linspace(0,1,100);
for i = 1:length(ReserveSizes)
   
   % What proportion is reserved?
   PR = ReserveSizes(i);
   
   % Run the B-H model to equilibrium with this level of reserves and effort outside
   [N(i),B(i),Ya(i),Yb(i),n_h,n_r,S,R] = sub_beverton_holt_model(PR,r,a,H_outside,f,A,mass,V);
   
   LengthDist_r(:,i) = n_r; % Length distribution in the reserved patches
   LengthDist_h(:,i) = n_h; % Length distribution in the harvested patches (scorched earth)
end
% Extract all the relevant values at the point of MSY
[MSY,I] = max(Yb); % The MSY biomass itself
MSY_a = Ya(I); % The number of individuals which make up the MSY
LF_r = LengthDist_r(:,I); % The length of individuals inside the reserves
LF_h = LengthDist_h(:,I); % The length of individuals in the fished habitat
SSB = B(I); % The equilibrium biomass at MSY
SSA = N(I); % The equilibrium abundance at MSY
Opt_reserve_proportion = ReserveSizes(I); % The optimal reserved proportion

% PLotting
if PT == 1
   CL = [H_outside 0 0]; LW = 2;
   if H_outside == 1; CL = 'g'; LW = 3; end
   subplot(1,3,1); hold on; plot(ReserveSizes,B,'-','linewidth',LW,'color',CL); axis tight
   subplot(1,3,2); hold on; plot(ReserveSizes,Ya,'-','linewidth',LW,'color',CL); axis tight
   subplot(1,3,3); hold on; plot(ReserveSizes,Yb,'-','linewidth',LW,'color',CL); axis tight
else
   % We need to make it negative it because FMINCON is a minimisation function
   MSY = -1*MSY;
end

%% ========================================================
%% ======= Calculate optimal traditional management =======
%% ========================================================
function [MSY,MSY_a,LF,SSB,SSA,Opt_h] = sub_TraditionalYield(A,f,a,r,mass,Length,V)

% What are all the possible harvest rates 
HVec = linspace(0,1,100);

for h = 1:length(HVec)
   
   % Alter the harvest rate
   H = HVec(h);
   
   % Run the B-H model to equilibrium with this level of harvest
   [N(h),B(h),Ya(h),Yb(h),n_h,n_r,S,R(h)] = sub_beverton_holt_model(0,r,a,H,f,A,mass,V);
   
   LengthDist(:,h) = n_h; % Length distribution in the harvested patches
   NormL = n_h./sum(n_h);
   MeanLength(h) = sum(NormL.*Length');
   
   if N(h) == 0
      N(h+1:length(HVec)) = 0; B(h+1:length(HVec)) = 0; R(h+1:length(HVec)) = 0; Ya(h+1:length(HVec)) = 0; Yb(h+1:length(HVec)) = 0; MeanLength(h:length(HVec)) = 0;
      break
   end
end
[MSY,I] = max(Yb);
MSY_a = Ya(I);
LF = LengthDist(:,I);
SSB = B(I);
SSA = N(I);
Opt_h = HVec(I);

subplot(1,3,1); hold on; plot(HVec,B,'b','linewidth',3);
title('Biomass','fontsize',15);
subplot(1,3,2); hold on; plot(HVec,Ya,'b','linewidth',3);
title('Yield (abundance)','fontsize',15);
subplot(1,3,3); hold on; plot(HVec,Yb,'b','linewidth',3);
title('Yield (biomass)','fontsize',15); L = legend('MPA+Scorched','MPAs-Effort','Effort only'); set(L,'fontsize',15,'box','off','location','best')

[~,I] = min(abs(B./B(1)-0.2));
Steepness = R(I)/R(1);
disp(['Steepness = ' num2str(Steepness,2)])
















