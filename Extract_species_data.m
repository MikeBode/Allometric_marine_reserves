clear all

% First extract the values from the appropriate spreadsheet
[D,T] = xlsread('Species_parameters.xlsx');

% How many species are there in this iteration of the spreadsheet?
NumSpp = size(D,1);

% Extract the names of each species and separate the Genus and Species binomial elements
Name = T(2:end,1);
for i = 1:NumSpp
   ThisName = Name{i};
   F = find(ThisName == ' ');
   Genus{i,1} = ThisName(1:F-1);
   Species{i,1} = ThisName(F+1:end);
end

% Extract the von Bertalanffy parameter K
K_vec = D(:,1);

% Extract the von Bertalanffy parameter L_infinity
L_inf_vec = D(:,2);

% Extract the allometric fertility exponent z
z_vec = D(:,4);

% Extract the allometric Length-Weight relationship parameters
LW_a_vec = D(:,6);
LW_b_vec = D(:,7);

% Extract the size at first reproduction
Min_f_length_vec = D(:,9);

% Save these as a Matlab .MAT file
save Species_parameters_matlab