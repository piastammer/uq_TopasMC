function [shift_rep] = getShifts(stf,D,mup,mue,sigmae,K,K_sum,correlation_type)
%Generates shifts in 3D coordinate system for given type of correlation
%between pencil beams

if strcmp(correlation_type,'full')
%% Perfectly correlated
shift_3D=quasimvnrnd(mue,diag(sigmae),D,'halton');

for i=1:D
shift_rep(i,:,:) = repelem(shift_3D(i,:),length(mup),1);
end

elseif strcmp(correlation_type,'ray')
%% Raywise
bixelsPerRay_accum= [0 cumsum([stf.numOfBixelsPerRay])];
shift_3D = zeros(D,3);
shift_rep = zeros(D,length(mup),3);
for i=1:stf.numOfRays
    shift_3D=quasimvnrnd(mue,diag(sigmae),D,'halton');
    for j=1:stf.numOfBixelsPerRay
     shift_rep(:,bixelsPerRay_accum(j)+i,:)=shift_3D;
    end
end

elseif strcmp(correlation_type,'beam')
%% Beamwise
shift_3D = zeros(D,3);
shift_rep = zeros(D,length(mup),3);
for r=1:size(rot_angles,1)
    shift_3D=quasimvnrnd(mue,diag(sigmae),D,'halton');
    for i=1:K(r+1)
     shift_rep(:,K_sum(r)+i,:)=shift_3D;
    end
end

elseif strcmp(correlation_type,'diag')
%% Independent
shift_rep = zeros(D,length(mup),3);
for i=1:length(mup)
    shift_rep(:,i,:)=quasimvnrnd(mue,diag(sigmae),D,'halton');
end

else 
    display('No known correlation type entered, please choose between "full","ray","beam" or "diag"')
end

end

