function [] = postProcess(dir_name,workspace_name,num_beams,cube_dims, pos_var, num_shifts, range_var, en_var)
%%Import and Post-process EventMaps from TOPAS MC Simulation with Binary input files
%%for radiation plan and input data from Matrad
% Input: dir_name       - directory where EventMaps are located
%        workspace_name - path to Matrad Workspace
%        num_beams      - number of beams/angles
%        cube_dims      - dimensions of Dose/ct-cube
%        pos_var        - variance in set-up/beam positions
%        num_shifts     - number of shifts to compute variance
%        range_var      - variance in range
%        en_var         - beam energy spread (not due to uncertainty)
%
% Output: nominal dose, expected dose and variance are saved to current directory

%Import Files
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

Data_Events_temp(1,:)= fread(fopen(fullfile(matfiles_ID{1,b}(i).folder, matfiles_ID{1,b}(i).name)),'int');
Data_Events_temp(2,:)= fread(fopen(fullfile(matfiles_numEntries{1,b}(i).folder, matfiles_numEntries{1,b}(i).name)),'int');
Data_Histories_temp(1,:)= fread(fopen(fullfile(matfiles_voxVal{1,b}(i).folder, matfiles_voxVal{1,b}(i).name)),'double');
Data_Histories_temp(2,:)=fread(fopen(fullfile(matfiles_voxIdx{1,b}(i).folder, matfiles_voxIdx{1,b}(i).name)),'int');
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

%% Compute weights according to given distributions and initial particle
%% positions
%Define necessary parameters
load(workspace_name)
display("Loaded Workspace")
StartingPoints=Data_Events(3:6,:)';


sigmae= ones(1,3).*pos_var;

%For Matrad Data
SAD_mm= machine.meta.SAD; %stf.SAD;
nozzleAxialDistance_mm=1500; %As set in interface script (matRad_exportMCinputFiles.m, line 37)

mup=[];
mue=[0 0 0];
pos_t=[];
energies=[];
K(1)=0;

%Derive means of beam position & energy
for i=1:numel(stf)
%stf(i).sourcePoint=stf(i).sourcePoint*100; %for more parallel beams

rayPos=reshape([stf(i).ray.rayPos],3,stf(i).numOfRays)';
pos =  [repelem(rayPos(:,1),stf(i).numOfBixelsPerRay,1) repelem(rayPos(:,2),stf(i).numOfBixelsPerRay,1) repelem(rayPos(:,3),stf(i).numOfBixelsPerRay,1)];
vec_dirs = rayPos - stf(i).sourcePoint;
vec_dirs = vec_dirs./vecnorm(vec_dirs,2,2);
vec_dirs=[repelem(vec_dirs(:,1),stf(i).numOfBixelsPerRay,1) repelem(vec_dirs(:,2),stf(i).numOfBixelsPerRay,1) repelem(vec_dirs(:,3),stf(i).numOfBixelsPerRay,1)];

mu_i = stf(i).sourcePoint + vec_dirs*(SAD_mm - nozzleAxialDistance_mm);
mup_i = rotateAxis(mu_i',stf(i).gantryAngle,stf(i).couchAngle);
[~,ix]=min(var(mup_i'));
rot_idx_bev(i,:)=setdiff([1 2 3],ix);
mup = horzcat(mup,squeeze(mup_i(rot_idx_bev(i,:),:)));
energies=horzcat(energies, stf(i).ray(:).energy);
K(i+1)=stf(i).totalNumOfBixels;
end

mu_energy=energies;

%Set standard deviations 
beamEnergySpread = en_var^2;
sigma_range=range_var;
sigma_energy=(sqrt(beamEnergySpread +(sigma_range/1.77)^2)/100*mu_energy).^2;
sigmae1=zeros(1,1,length(mu_energy));
sigmae2=zeros(1,1,length(mu_energy));
for i=1:length(mu_energy)
sigmae1(:,:,i)=sigma_energy(i);
sigmae2(:,:,i)=beamEnergySpread*(mu_energy(i)^2)/(100^2);
end

sigmap=zeros(1,2,length(mup));
sigmag=zeros(1,2,length(mup));

[~,idx_ray]=ismember(energies,[machine.data.energy]); 
for i=1:numel(idx_ray)
sigmap(1,:,i)=machine.data(idx_ray(i)).initFocus.sigma(1,:).^2;
sigmag(1,:,i)=machine.data(idx_ray(i)).initFocus.sigma(1,:).^2+sigmae(1:2);
end
numBixelsPerRay= [0 [stf.numOfBixelsPerRay]]';

D=num_shifts; %% Set number of shifts
M=length(StartingPoints);

rot_angles=[[stf.gantryAngle]' [stf.couchAngle]'];

%Define distributions
sampleDist = my_gmm(mup,sigmag,rot_angles,mu_energy,sigmae1,rot_idx_bev,resultGUI.w);
targetDistNom = my_gmm(mup,sigmap,rot_angles,mu_energy,sigmae2,rot_idx_bev,resultGUI.w);
targetDistExp = my_gmm(mup,sigmag,rot_angles,mu_energy,sigmae1,rot_idx_bev,resultGUI.w);

%% Compute nominal and expected Dose
%Generate random shifts 
shift=[];
shift_rep=[];
K_sum=cumsum(K);

shift_3D=quasimvnrnd(mue,diag(sigmae),D,'halton');

D=1;
for i=1:D
shift_rep(i,:,:) = repelem(shift_3D(i,:),length(mup),1);
end
for i=1:D
for r=1:size(rot_angles,1)
    tmp = rotateAxis(squeeze(shift_rep(i,K_sum(r)+1:K_sum(r+1),:))',rot_angles(r,1),rot_angles(r,2))';
    shift(i,K_sum(r)+1:K_sum(r+1),:) = tmp(:,rot_idx_bev(r,:));
end
end

shift_energy_prct=normrnd(0,sigma_range/1.77,D,1);
for i = 1:length(mu_energy)
shift_energy(:,i)=shift_energy_prct.*mu_energy(i)/100;
end
tic

%Identify indices of hit voxels and transform from linear to subscript
[idx_lin,~,subs]=unique(Data_Histories(2,:)+1,'stable');
[idx(:,1),idx(:,2),idx(:,3)]=ind2sub(cube_dims,idx_lin);

%Compute weights and accumulate weighted events in hit voxels to Dose and
%variance estimates
Sample_prob=sampleDist.prob(StartingPoints,K,n_beams);

W_exp=targetDistExp.prob(StartingPoints,K,n_beams)./Sample_prob;
W_exp=repelem(W_exp,Data_Events(2,:));
Dose_exp=full(ndSparse.build(idx,accumarray(subs,W_exp.*Data_Histories(1,:)'),cube_dims));

W_nom=targetDistNom.prob(StartingPoints,K,n_beams)./Sample_prob;
W_nom=repelem(W_nom,Data_Events(2,:));
Dose_nom=full(ndSparse.build(idx,accumarray(subs,W_nom.*Data_Histories(1,:)'),cube_dims));

display("Nominal and expected dose computed successfully!")

%% Compute Variance
%load(ref_name)
%Var_ref=arr;
Dose_mean=zeros(length(idx),1);
Var =Dose_mean;
Var_w =zeros(length(Data_Histories),1);
Mean_w=Var_w;
%error_conv = zeros(D,1);

%Set shifted distribution and precompute sample distribution
shiftedDist=my_gmm(mup,sigmap,rot_angles,mu_energy,sigmae2,rot_idx_bev,resultGUI.w);

for i=1:D
W_shifts=[];
shiftedDist.mu=mup +squeeze(shift(i,:,:))';
%shiftedDist.mu_energy=mu_energy+shift_energy(i,:); %For energy/range uncertainties

W_shifts=shiftedDist.prob(StartingPoints,K,n_beams)./Sample_prob;
W_shifts=repelem(W_shifts,Data_Events(2,:));

Dose_w=accumarray(subs,W_shifts.*Data_Histories(1,:)');

Dose_mean_old=Dose_mean;
Dose_mean=Dose_mean_old+(Dose_w-Dose_mean_old)./i;
Var= Var + (Dose_w - Dose_mean_old).*(Dose_w - Dose_mean);



%Var_full=full(ndSparse.build(idx,Var,cube_dims))./i;
%error_conv(i)= mean((Var_full-Var_ref).^2,'all');
end

%Build full 3D result matrix for variance
Var_full=full(ndSparse.build(idx,Var,cube_dims))./i;

display("Variance computed successfully!")
save('Estimate_NomDose','Dose_nom')
save('Estimate_ExpDose','Dose_exp')
save('Estimate_Variance','Var_full')
%save('Error_Convergence','error_conv')

display("Saved results, exiting...")
end
