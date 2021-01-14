function [] = postProcess(dir_name,workspace_name,num_beams,cube_dims, pos_var, num_shifts, range_var, en_var)
%%Import and Post-process EventMaps from TOPAS MC Simulation with Binary input files
%%for radiation plan and input data from Matrad
% Input: dir_name       - directory where EventMaps are located [string]
%        workspace_name - path to Matrad Workspace [string]
%        num_beams      - number of beams/angles [scalar]
%        cube_dims      - dimensions of Dose/ct-cube [3x1 vector]
%        pos_var        - variance in set-up/beam positions (due to uncertainty) [3x1 vector]
%        num_shifts     - number of shifts to compute variance [scalar]
%        range_var      - variance in range (due to uncertainty) [scalar]
%        en_var         - beam energy spread (not due to uncertainty)[scalar]
%
% Output: nominal dose, expected dose and variance are saved to current directory

%Import Files
tic
display("Started post processing...")
B=num_beams; %number of beams
testfiledir = dir_name;

for b=1:B
matfiles_ID{b} = dir(fullfile(testfiledir, sprintf('*_ID_*_%i.bin',b)));
matfiles_voxVal{b} = dir(fullfile(testfiledir, sprintf('*_voxVal_*_%i.bin',b)));
matfiles_voxIdx{b} = dir(fullfile(testfiledir, sprintf('*_voxIdx_*_%i.bin',b)));
matfiles_numEntries{b} = dir(fullfile(testfiledir, sprintf('*_numEntries_*_%i.bin',b)));
matfiles_initVert{b} = dir(fullfile(testfiledir, sprintf('*_initVert_*_%i.bin',b)));
matfiles_initEnergy{b} = dir(fullfile(testfiledir, sprintf('*_initEnergy_*_%i.bin',b)));
nfiles{b} = length(matfiles_ID{b});
end


n_beams = zeros(B+1,1);

Data_Events=[];
Data_Histories=[];

for b=1:B
for i=1:nfiles{1,b}
Data_Events_temp=[];
Data_Histories_temp =[];

Data_Events_temp(1,:)=uint16(fread(fopen(fullfile(matfiles_ID{1,b}(i).folder, matfiles_ID{1,b}(i).name)),'int'));
Data_Events_temp(2,:)= uint16(fread(fopen(fullfile(matfiles_numEntries{1,b}(i).folder, matfiles_numEntries{1,b}(i).name)),'int'));
Data_Histories_temp(1,:)= fread(fopen(fullfile(matfiles_voxVal{1,b}(i).folder, matfiles_voxVal{1,b}(i).name)),'double');
Data_Histories_temp(2,:)=uint32(fread(fopen(fullfile(matfiles_voxIdx{1,b}(i).folder, matfiles_voxIdx{1,b}(i).name)),'int'));
Data_Events_temp(3:5,:)= reshape(fread(fopen(fullfile(matfiles_initVert{1,b}(i).folder, matfiles_initVert{1,b}(i).name)),'double'),3,[]);
Data_Events_temp(6,:)=fread(fopen(fullfile(matfiles_initEnergy{1,b}(i).folder, matfiles_initEnergy{1,b}(i).name)),'double');

Data_Events=horzcat(Data_Events,Data_Events_temp);
Data_Histories=horzcat(Data_Histories, Data_Histories_temp);
end
n_beams(b+1) = length(Data_Events) - n_beams(b);
end

n_beams=cumsum(n_beams);
fclose('all');

display("Loaded EventMaps")
display("Number of histories = "+ length(Data_Events));

%% -- Get parameters from input files -- %%
%Load MatRad workspace
load(workspace_name)
display("loaded Workspace")

% Initial positions of particles
StartingPoints=Data_Events(3:6,:)';
Data_Events=Data_Events(1:2,:);

% Set machine parameters
SAD_mm= machine.meta.SAD; 
nozzleAxialDistance_mm=1500; %As set in interface script (matRad_exportMCinputFiles.m, line 37)

%Initialise some variables
mup = [];         % Mean of particle distribution (without uncertainty)
mue = [0 0 0];    % Mean of uncertainty distribution (= 0, so on average there is no error)   
mu_energy = [];    % Vector for mean of energy distribution
K(1)=0;         % K(i+1) is number of pencil beams in beam i
sigmae= ones(1,3).*pos_var; % Variance of uncertain parameter

%Derive means of beam position & energy
for i=1:B
%Source point of particles
stf(i).sourcePoint=stf(i).sourcePoint*100; %for more parallel beams

%get direction of beams (vector from source to target point)
rayPos=reshape([stf(i).ray.rayPos],3,stf(i).numOfRays)';
vec_dirs = rayPos - stf(i).sourcePoint;
vec_dirs = vec_dirs./vecnorm(vec_dirs,2,2);
vec_dirs=[repelem(vec_dirs(:,1),stf(i).numOfBixelsPerRay,1) repelem(vec_dirs(:,2),stf(i).numOfBixelsPerRay,1) repelem(vec_dirs(:,3),stf(i).numOfBixelsPerRay,1)];

%mean of initial positions is given at distance (SAD_mm - nozzleAxialDistance_mm)
%from source point
mu_i = stf(i).sourcePoint + vec_dirs*(SAD_mm - nozzleAxialDistance_mm);

%rotate into beams eye view (one of the axis is parallel to beam direction)
mup_i = rotateAxis(mu_i',stf(i).couchAngle,stf(i).gantryAngle);
if min(var(mup_i,[],2)) > 0
    [~,ix]=min(var(mup_i,[],2));
    rot_idx_bev(i,:)=setdiff([1 2 3],ix);
else
    rot_idx_bev(i,:)=[1 3];
end
mup = horzcat(mup,squeeze(mup_i(rot_idx_bev(i,:),:)));

%get mean of energy distribution
mu_energy=horzcat(mu_energy, stf(i).ray(:).energy);

%Fill in K (number of pencil beams in beam i)
K(i+1)=stf(i).totalNumOfBixels;

%get beam std deviation
sigmae_bev(i) = {rotateAxis(sqrt(sigmae'),stf(i).couchAngle,stf(i).gantryAngle).^2};
end

%Set standard deviations 
%From input:
beamEnergySpread = en_var^2;
sigma_range=range_var;

%Range-energy relation
sigma_energy=(sqrt(beamEnergySpread +(sigma_range/1.77)^2)/100*mu_energy).^2;
sigmae1=zeros(1,1,length(mu_energy)); %energy spread without uncertainty
sigmae2=zeros(1,1,length(mu_energy)); %energy spread with uncertainty
for i=1:length(mu_energy)
sigmae1(:,:,i)=sigma_energy(i);
sigmae2(:,:,i)=beamEnergySpread*(mu_energy(i)^2)/(100^2);
end

%Initialise variables
sigma_beam={}; %beam variance without uncertainty
sigma_total={}; %beam variance with uncertainty
sigmap=zeros(1,2,length(mup));
sigmag=zeros(1,2,length(mup));

%find out which pencil beams are used in computation
[~,idx_ray]=ismember(mu_energy,[machine.data.energy]); 

for b=1:numel(stf)
for i=1:numel(idx_ray)
sigmap(1,:,i)=machine.data(idx_ray(i)).initFocus.sigma(1,:).^2; % "natural" spread of particles
sigmag(1,:,i)=machine.data(idx_ray(i)).initFocus.sigma(1,:).^2+sigmae_bev{b}([3 1])'; % convolution with std. dev. of uncertain parameter 
end
sigma_beam{b}=sigmap;
sigma_total{b}=sigmag;
end


D=num_shifts; % Set number of shifts

% Angles of beams
rot_angles=[[stf.couchAngle]' [stf.gantryAngle]'];

%Define distributions
sampleDist = my_gmm(mup,sigma_beam,rot_angles,mu_energy,sigmae1,rot_idx_bev,resultGUI.w); % Distribution of simulated particles
targetDistNom = my_gmm(mup,sigma_beam,rot_angles,mu_energy,sigmae2,rot_idx_bev,resultGUI.w); % Distribution of particles without uncertainty
targetDistExp = my_gmm(mup,sigma_total,rot_angles,mu_energy,sigmae1,rot_idx_bev,resultGUI.w); %Distribution of particles with uncertainty
                                                                                              
%%Clear matRad  workspace 
clearvars ct cst pln machine dij;
resultGUI=rmfield(resultGUI,{'info','physicalDose_beam1','wUnsequenced','usedOptimizer'});


%% -- Sample random variables -- %%

%Initialise + helper variables
shift=[];
shift_rep=[];
K_sum=cumsum(K);

%Get realisations for random set-up error
shift_rep=getShifts(stf,D,mup,mue,sigmae,K,K_sum,'full');

%%Clear stf
clearvars stf;

%Rotate shifts into right coordinate system (in beams eye views of
%different beams)
for i=1:D
for r=1:B
   tmp = rotateAxis(squeeze(shift_rep(i,K_sum(r)+1:K_sum(r+1),:))',rot_angles(r,1),rot_angles(r,2))';
   shift(i,K_sum(r)+1:K_sum(r+1),:) = tmp(:,[1 3]);
end
end

% Random energy scaling (only needed for range uncertainties)
% shift_energy_prct=normrnd(0,sigma_range/1.77,D,1);
% for i = 1:length(mu_energy)
% shift_energy(:,i)=shift_energy_prct.*mu_energy(i)/100;
% end

%% -- Dose and variance computation -- %%
%Compute weights and accumulate weighted events in hit voxels to Dose and
%variance estimates

%Identify indices of hit voxels and transform from linear to subscript

%This is needed if arrays are too large (instead of line 214)
%{ 
parpool('local',8,'SpmdEnabled',false)
A=tall(Data_Histories(2,:)'+1);
Data_Histories=Data_Histories(1,:);

[i,c,~]=unique(A);
[~,t2]=sort(c,1);
idx_lin=i(t2);
idx_lin=gather(idx_lin);
clearvars i c t2;
[~,subs]=ismember(A,idx_lin);
subs=gather(subs);
clear A;
poolobj = gcp('nocreate');
delete(poolobj);

display("Parallel computations done!") 
%}

[idx_lin,~,subs]=unique(Data_Histories(2,:)'+1,'stable');  
[idx(:,1),idx(:,2),idx(:,3)]=ind2sub(cube_dims,idx_lin);

%Probability of sample in original (simulated) distribution
Sample_prob=sampleDist.prob(StartingPoints,K,n_beams);

%% -- Expected Dose -- %%
%Weight = probability of sample in new/target dist. divided by prob. in
%original distribution
W_exp=targetDistExp.prob(StartingPoints,K,n_beams)./Sample_prob;

%Build up 3D dose cube
W_exp=repelem(W_exp,Data_Events(2,:));
Dose_exp=full(ndSparse.build(idx,accumarray(subs,W_exp.*Data_Histories(1,:)'),cube_dims));

%Save and clear
save('Estimate_ExpDose','Dose_exp')
clearvars W_exp Dose_exp;

%% -- Nominal Dose -- %%
%Weight = probability of sample in new/target dist. divided by prob. in
%original distribution
W_nom=targetDistNom.prob(StartingPoints,K,n_beams)./Sample_prob;

%Build up 3D dose cube
W_nom=repelem(W_nom,Data_Events(2,:));
Dose_nom=full(ndSparse.build(idx,accumarray(subs,W_nom.*Data_Histories(1,:)'),cube_dims));

%Save and clear
save('Estimate_NomDose','Dose_nom')
clearvars W_nom Dose_nom;
display("Nominal and expected dose computed successfully!")

%% -- Variance -- %%
% Initialise some variables
Dose_mean=zeros(length(idx),1);
Var = Dose_mean;
Var_w = zeros(length(Data_Histories),1);
Mean_w = Var_w;

%Create shifted distribution
shiftedDist=my_gmm(mup,sigma_beam,rot_angles,mu_energy,sigmae2,rot_idx_bev,resultGUI.w);

for i=1:D
W_shifts=[];

%Shift distribution according to ith realisation
shiftedDist.mu = mup + squeeze(shift(i,:,:))';

%shiftedDist.mu_energy = mu_energy + shift_energy(i,:); %Only for energy/range uncertainties

%Compute weight
W_shifts=shiftedDist.prob(StartingPoints,K,n_beams)./Sample_prob;
W_shifts=repelem(W_shifts,Data_Events(2,:));

%Accumulate shifted dose
Dose_w=accumarray(subs,W_shifts.*Data_Histories(1,:)');

%Update mean and variance 
Dose_mean=Dose_mean+(Dose_w-Dose_mean)./i;
Var= Var + (Dose_w - Dose_mean).*(Dose_w - Dose_mean);
end

%Build full 3D result matrix for variance
Var_full=full(ndSparse.build(idx,Var,cube_dims))./i;

display("Variance computed successfully!")
save('Estimate_Variance','Var_full')

toc
display("Saved results, exiting...")
end
