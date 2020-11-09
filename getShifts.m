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
numOfRays = [0 [stf(:).numOfRays]];
totalNumOfRays=sum(numOfRays);
shift_3D=mvnrnd(mue,diag(sigmae),totalNumOfRays*D);
shift_rep = zeros(D,length(mup),3);

for r=1:numel(stf)
    bixelsPerRay_accum= [0 cumsum([stf(r).numOfBixelsPerRay])];
for i=1:stf(r).numOfRays
    for j=1:stf(r).numOfBixelsPerRay
     shift_rep(:,bixelsPerRay_accum(j)+i,:)=shift_3D((r-1)*D*numOfRays(r)+D*(i-1)+1:(r-1)*D*numOfRays(r)+D*(i-1)+D,:);
    end
end
end

elseif strcmp(correlation_type,'beam')
%% Beamwise
shift_3D=mvnrnd(mue,diag(sigmae),numel(stf)*D);
shift_rep = zeros(D,length(mup),3);
for r=1:numel(stf)
    for i=1:K(r+1)
     shift_rep(:,K_sum(r)+i,:)=shift_3D((r-1)*D+1:(r-1)*D+D,:);
    end
end

elseif strcmp(correlation_type,'diag')
%% Independent
shift_rep = zeros(D,length(mup),3);
shift_3D=mvnrnd(mue,diag(sigmae),length(mup)*D);
for i=1:length(mup)
    shift_rep(:,i,:)=shift_3D((i-1)*D+1:(i-1)*D+D,:);
end

else 
    display('No known correlation type entered, please choose between "full","ray","beam" or "diag"')
end

end

