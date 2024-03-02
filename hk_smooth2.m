function output=hk_smooth2(input,surface,sigma,n_smooth)
%output=hk_smooth2(input,surface,sigma,n_smooth)
%
%The function performs heat kernel smoothing on a triangle mesh [1][2][3].
% 
% input         : Multiple signals to smooth.
% surface       : Structured array consisting of surface.vertices and surface.faces.
%                   The default MATLAB data strcture for isosurface algorithm
%                   is needed
% sigma         : bandwidth. The bandwidth corresponds to diffuion time in a heat equation [3]. 
%                   The bandwidth given in old publications [1] [2] are related to 
%                   new publication [3] as our new_sigma = old_sigma^2/2.
% n_smooth      : number of iterations
%
%
% EXAMPLE: bandwidth sigma=1, number of iteration=50
% output_signal=hk_smooth(input_signal,tri,coord,nbr,1,50);
%
%
% The code is downloaed from http://www.stat.wisc.edu/softwares/hk/hk.html
%
% (C) 2004- Moo K. Chung, mkchung@wisc.edu
%  Department of Biostatisics and Medical Informatics
%  University of Wisconsin, Madison
%
% If you use this code, please reference one of the following papers. The details on
% the mathematical basis of of the algorithm can be found in these papers.
%
% [1] Chung, M.K., Robbins,S., Dalton, K.M., Davidson, Alexander, A.L., R.J., Evans, A.C. 2005. 
%     Cortical thickness analysis in autism via heat kernel smoothing. NeuroImage 25:1256-1265
%     http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf
%
% [2] Chung, M.K., Robbins, S., Evans, A.C. 2005. Unified statistical approach to cortical 
%     thickness analysis. Information Processing in Medical Imaging (IPMI). Lecture Notes in 
%     Computer Science (LNCS) 3565:627-638. Springer-Verlag
%     http://www.stat.wisc.edu/~mchung/papers/IPMI/hk.IPMI.2005.pdf
%
% [3] Chung, M.K. Hartley, R., Dalton, K.M., Davidson, R.J. 2008. Encoding
%     cortical surface by spherical harmonics.  Satistica Sinica 18:1269-1291
%     http://www.stat.wisc.edu/%7Emchung/papers/sinica.2008.pdf

% Update history: 
%2004 Feb 5:
%2007 Sept 18:
%2010 March 20:
%2012 May 29: FINDnbr can compute planer mesh/graph
%2012 Nov 10: function imports multiple signals. 
%2013 Feb 15: If input is [], it smooths surface instead.
%2018 January 1: Fixed inconsistent degree for some nodes. 

% HEAT KERNEL SHAPE

coord=surface.vertices;
tri=surface.faces;
K=inline('exp(-x/(4*sigma))/sum(exp(-x/(4*sigma)))');

% For the convension used in publications [1] and [2], change this line to  
% K=inline('exp(-x/(2*sigma^2))/sum(exp(-x/(2*sigma^2)))');

[nbr, deg] = FINDnbr(surface.faces);

%HEAT KERNEL WEIGHT COMPUTATION
n_vertex=size(nbr,1);
max_degree=max(size(nbr,2));
weight=zeros(n_vertex,max_degree+1); %weight including the current vertex

for i=1:n_vertex
    degree=deg(i);
    current_nbr=nbr(i,nbr(i,:)~=0);
    pipo=[coord(current_nbr,:)-kron(ones(degree,1),coord(i,:))]';
    distance=sum(pipo.*pipo);
    weight(i,1:(1+degree))=K(sigma,[0 distance]);
end;

% ITERATIVE KERNEL SMOOTHING
% increased the index by one to account for 0 indexing in nbr.
% signal(1,:) will be a dummy entry corresponding to 0 indexing.

if isempty(input)
    input = surface.vertices;
end;
    
nos=size(input,2); %number of signals
output=zeros(n_vertex,nos);

for k=1:nos
    signal=zeros(n_vertex+1,1);
    signal(2:n_vertex+1)=input(:,k);
    dummy=ones(1,max_degree);
    for j=1:n_smooth
        Y=[signal signal([dummy; nbr+1])];  %1st row is a dummy entry
        signal(2:n_vertex+1)= sum(weight.*Y(2:n_vertex+1,:),2);
    end;
    output(:,k) =signal(2:n_vertex+1,:);
end;


%------------------------------
% function calls

function [nbr, degree] = FINDnbr(tri)
%[nbr, degree] = FINDnbr(tri)
%
% The function finds the first ring neighbors and degree from a given a triangle mesh.  
% tri    : list of triangle elements
% nbr    : first order neihbor vertex list 
% degree : degree of each vertex
%
%(C) Moo K. Chung   mkchung@wisc.edu
% Department of Statistics and Medical Informatics
% University of Wisconsin-Madison
%
% http://www.stat.wisc.edu/~mchung/softwares/softwares.html
% Update history:
%        Feb 5, 2004; Nov 14, 2007; 
%        July 24, 2010: It works for planer mesh as well
%2018 January 1: Fixed inconsistent degree for some nodes. 

%degree computation

n_tri=size(tri,1);
%n_vertex = (n_tri+4)/2;
% this equation only works for spherical object
%
% V = number of vertices
% F = number of faces
% For spherical mesh we have
% F = 2V - 4: Euler characteristic equation.


% maximum node index becomes the number of vertices
n_vertex=max(max(tri));

if (2*n_vertex-n_tri) ==4
    % this only works for closed spherical mesh
    % it uses the fact the total number of edge is counted twice
    degree=zeros(n_vertex,1);
    for j=1:n_tri
        degree(tri(j,:))=degree(tri(j,:))+1;
    end
    max_degree=max(degree);

    % 1st order neighbor computation

    nbr=zeros(n_vertex,max_degree);

    for i_tri=1:n_tri
        for j=1:3
            cur_point = tri(i_tri,j);
            for k=1:3
                if (j ~= k)
                    nbr_point= tri(i_tri,k);
                    if find(nbr(cur_point,:)==nbr_point)
                        ;
                    else
                        n_nbr = min(find(nbr(cur_point,:) == 0));
                        nbr(cur_point,n_nbr) = nbr_point;
                    end;
                end;
            end;
        end;
    end;
    
else
    %planer mesh
    adj = mesh2adj(tri); % computes the adjaceny matrix
    degree= sum(adj,2); % sum the adjacency matrix to compute the node degree

    max_degree=max(degree);
    nbr=zeros(n_vertex,max_degree);
    for i=1:n_vertex
        nbr(i,1:degree(i))=find(adj(i,:));
    end;
end;

degree = sum(nbr>0,2); %fix inconsistent node degree


%--------------------------------------
function adj = mesh2adj(tri)
%adj = mesh2adj(tri)
%
% MESH2ADJ computes the adjacency matrix of given triangulation
%
% tri    : list of triangle elements that form triangles
% adj    : adjacency matrix
%
%(C) Moo K. Chung   mkchung@wisc.edu
%
% Department of Statistics and Medical Informatics
% University of Wisconsin-Madison
%
%
% Update history:
% July 24, 2010: function created


n_tri=size(tri,1);  % number of triangles
n_vertex = max(max(tri));  % number of nodes

adj=zeros(n_vertex,n_vertex);

for i=1:n_tri
    adj(tri(i,1), tri(i,2))=1;
    adj(tri(i,2), tri(i,3))=1;
    adj(tri(i,3), tri(i,1))=1;
end;


adj=(adj'+adj); % make it symmetric

[I, J] = find(adj);

% renormalize
for i=1:length(I)
    adj(I(i),J(i)) = 1;
end;



% OLD CODE OPTIMIZED FOR MNI FORMAT
%discrete heat kernel
%K=inline('exp(-x/(4*sigma))/sum(exp(-x/(4*sigma)))');
% For publications [1] and [2], change this line to  K=inline('exp(-x/(2*sigma^2))/sum(exp(-x/(2*sigma^2)))');

%kernel weight computation
%n_vertex=size(nbr,1);
%weight=zeros(n_vertex,7);  % excluding the center, there are maximum total of 6 neighbors.
%
%for i=1:n_vertex
    % the first 12 vertices have only 5 neighbors
%     if (i<=12)
%         current_nbr=nbr(i,1:5);
%         pipo=coord(:,current_nbr(1:5))-kron(ones(1,5),coord(:,i));
%         distance=sum(pipo.*pipo);
%         weight(i,:)=[K(sigma,[0 distance]) 0];
%     else
%         % the remaining vertices have 6 neighbors
%         current_nbr=nbr(i,1:6);
%         pipo=coord(:,current_nbr(1:6))-kron(ones(1,6),coord(:,i));
%         distance=sum(pipo.*pipo);
%         weight(i,:)=K(sigma,[0 distance]);
%     end;
% end;
% 
% 
% %iterated smoothing
% signal=input_signal;
% for j=1:n_smooth
%     Y=[signal [[signal(nbr(1:12,1:5)) zeros(12,1)]; signal(nbr(13:n_vertex,:))]]; % signal + neighboring signal
%     signal=sum(weight.*Y,2);  %weighed averaging
% end;



