function output= diffusion_smooth3(input, surf, sigma, A,C)
%
% output= diffusion_smooth3(input, surf, sigma)
%
%
% input         : signal to smooth.
% surf          : structured array consisting of surface.vertices and surface.faces.
%                 The default MATLAB data strcture for isosurface algorithm is used
% sigma         : bandwidth. The bandwidth corresponds to time in a heat equation [3]. 
%                   The bandwidth given in old publications [1] [2] are related to 
%                   new publication [3] as our new_sigma = old_sigma^2/2.
% A,C            : A and C matrices are computed from FEM.m. 
%
%
% (C) 2014 Moo K. Chung
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%  
% email://mchung@stat.wisc.edu
% http://www.stat.wisc.edu/~mchung/softwares/diffusion/diffusion.html
%
%
% If you use this code, please reference [1], [2] or [3]. The details on
% the mathematical basis of of the algorithm can be found in these papers.
% The exact code is based in the simulation study given in [3].
%
% [1] Chung, M.K. 2001. Statistical Morphometry in Neuroanatomy, 
%     PhD Thesis, McGill University.
%     http://www.stat.wisc.edu/~mchung/papers/thesis.pdf
%
% [2] Chung, M.K., Taylor, J. 2004. Diffusion Smoothing on Brain Surface via Finite 
%     Element Method,  IEEE International Symposium on Biomedical Imaging (ISBI). 562.
%     http://www.stat.wisc.edu/~mchung/papers/BMI2004/diffusion_biomed04.pdf
%
% [3] Chung, M.K., Qiu, A., Seo S. Vorperian, H.K. 2015. Unified heat kernel regression 
%     for diffusion, kernel smoothing and wavelets on manifolds and its application to 
%     mandible growth modeling in CT images, Medical Image Analysis. 22:63-76
%     http://www.stat.wisc.edu/%7Emchung/papers/chung.2015.MIA.pdf
%
% Update history. 
% 2019 March 4, additional comments added. 
%
%
%

%Determining the mesh topology
coord=surf.vertices;
tri=surf.faces;
[nbr, deg] = FINDnbr(tri);

%DISCRETE LAPLACE-BELTRAMI OPERATOR ESTIMATION
%based on [1] and [2]. For large matrices, LB-operator should be determined using cotan formualtion

laplace= -inv(A)*C;

% number of iterations (has to be large)
n_time=100;  % 
% step size (has to be small)
delta=sigma/n_time;

%figure_trimesh(surf,input)

curr_signal=input;


for i_time=1:n_time
    maxcurr=max(curr_signal);
    mincurr=min(curr_signal);

    curr_signal= curr_signal + delta*laplace*curr_signal;
    
    ind= find(curr_signal>=maxcurr);
    curr_signal(ind)= maxcurr;
    
    ind= find(curr_signal<=mincurr);
    curr_signal(ind)= mincurr;
    
    %figure; figure_surf(surf, curr_signal); colorbar
end;

output=curr_signal;

%figure_trimesh(surf,curr_signal)
