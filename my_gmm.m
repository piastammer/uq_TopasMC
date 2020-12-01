classdef my_gmm
%Implements a gaussian mixture model where the components are
%defined in different rotated coordinate systems

    properties
        mu     % Kxd matrix with means of components (means of positions)
        sigma  % Kxdxd/Kxdx1 matrix with covariance matrices/diagonal elements (variance of positions)
        angles % Kx2 rotation angles for components (beam angles)
        mu_energy %(Number of bixels)x1 (mean of energies)
        sigma_energy %scalar (variance of energies)
        w %weight vector (weights of components in mixture model)
        idx %indices to keep after rotation
    end
    
    methods
        function obj = my_gmm(m,s,w,m_e,s_e,rot_idx,weights)
            %Constructor that sets parameters
            obj.mu = m;
            obj.sigma = s;
            obj.angles = w;
            obj.mu_energy = m_e;
            obj.sigma_energy = s_e;
            obj.idx = rot_idx;
            if nargin<7
                obj.w = ones(length(m_e),1)./length(m_e);
            else
            obj.w = weights/sum(weights);
            end
        end
        
        function p = prob(obj,x,K,n_beams)
            %Computes probabilities of input points x (Nxd)
            B=size(obj.angles,1);
            N=length(x);
            p=zeros(N,1);
            if obj.sigma_energy > 0
            for i=1:B
                x_rot=rotateAxis(x(n_beams(i)+1:n_beams(i+1),1:3)',obj.angles(i,1),obj.angles(i,2))';
                x_rot2D=x_rot(:,obj.idx(i,:));
                x_rot3D=[x_rot2D x(n_beams(i)+1:n_beams(i+1),4)];
                clear x_rot;
%Computation with matrix multiplication (not possible, matrices get too large)
%                 z_1 = (x_rot2D(:,1) - obj.mu(1,K(i)+1:K(i+1)))*(1./sqrt(obj.sigma(:,1,K(i)+1:K(i+1))));
%                 z_2 = (x_rot2D(:,2) - obj.mu(2,K(i)+1:K(i+1)))*(1./sqrt(obj.sigma(:,2,K(i)+1:K(i+1))));
%                 p_1 = normpdf(z_1);
%                 p_2 = normpdf(z_2);
%                 p(n_beams(i)+1:n_beams(i+1))= sum(mvnpdf(x_rot2D,obj.mu(:,K(i)+1:K(i+1))',obj.sigma(:,:,K(i)+1:K(i+1))) .* obj.w(K(i)+1:K(i+1)));
                for j=1:K(i+1)
                    p(n_beams(i)+1:n_beams(i+1)) = p(n_beams(i)+1:n_beams(i+1)) + mvnpdf(x_rot3D,[obj.mu(:,(i-1)*K(i)+j); obj.mu_energy((i-1)*K(i)+j)]',[obj.sigma{i}(:,:,(i-1)*K(i)+j) obj.sigma_energy(:,:,(i-1)*K(i)+j)])*obj.w((i-1)*K(i)+j);
                end
            end
	    
            else
		for i=1:B
                x_rot=rotateAxis(x(n_beams(i)+1:n_beams(i+1),1:3)',obj.angles(i,1),obj.angles(i,2))';
                x_rot2D=x_rot(:,obj.idx(i,:));
                clear x_rot;

                for j=1:K(i+1)
               	 p(n_beams(i)+1:n_beams(i+1)) = p(n_beams(i)+1:n_beams(i+1)) + mvnpdf(x_rot2D,obj.mu(:,(i-1)*K(i)+j)',obj.sigma{i}(:,:,(i-1)*K(i)+j))*obj.w((i-1)*K(i)+j).*ismember(x(n_beams(i)+1:n_beams(i+1),4),obj.mu_energy((i-1)*K(i)+j));
                end
                end
            end
        end
    end
end

