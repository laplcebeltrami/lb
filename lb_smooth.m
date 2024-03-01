function output=lb_smooth(input,surf,sigma, k, V, D)
%output=lb_smooth(input,surf,sigma, k, V, D)
%
% The function performs heat kernel smoothing using the eigenfunctions 
% of the Laplace-Beltrami operator on a triangle mesh. 
% 
% input       : Signal to smooth.
% surf        : Structured array consisting of surf.vertices and surf.faces.
%               The default MATLAB data strcture for isosurface algorithm
%               is needed
% sigma       : bandwidth. The bandwidth corresponds to diffuion time in a heat equation [3]. 
% k           : number of basis used
% V           : eigenfunctions of the Laplace-Beltrami operator
% D           : diag(D) gives the eigenvalues of the Laplece-Beltrami
%               operator in increasing order.
%
%
% The code was downloaed from http://brainimaging.waisman.wisc.edu/~chung/lb
%
%
% (C) Moo K. Chung, Seongho Seo
%
%  email://mkchung@wisc.edu
%
%  Department of Biostatisics and Medical Informatics
%  University of Wisconsin, Madison
%
%  Department of Brain and Cognitive Sciences
%  Seoul National University
%
% 
% If you use this code, please reference [3]. The details on
% the mathematical basis of of the algorithm can be found in these papers.
%
% [1] Chung, M.K. 2001. Statistical Morphometry in Neuroanatomy, 
%     PhD Thesis, McGill University.
%     http://www.stat.wisc.edu/~mchung/papers/thesis.pdf
%
% [2] Chung, M.K., Taylor, J. 2004. Diffusion Smoothing on Brain Surface via Finite 
%     Element Method,  IEEE International Symposium on Biomedical Imaging (ISBI). 562.
%     http://www.stat.wisc.edu/~mchung/papers/BMI2004/diffusion_biomed04.pdf
%
% [3] Seo, S., Chung, M.K., Vorperian, H. K. Heat kernel smoothing of anatomical
%     manifolds via Laplace-Beltrami eigenfunctions. submitted.
%
% Update history: April 23, 2010.

% At this momement it smooth out surface and input is simply [].

if isempty(input)
    p=surf.vertices;
else
    p=input;
end;

Psi=V(:,1:k); %p = Psi * beta
eigen=diag(D);
W= exp(-1*eigen(1:k)*sigma);

%W= repmat(W', size(Psi,1), 1);
W= repmat(W', size(Psi,1), 1);

% We will also assume sigma=0, in which case we have eigenfunction
% expansion.
beta= inv(Psi'*Psi)*Psi'*p;

phat=(W.*Psi)*beta;

if isempty(input)
    output.vertices=phat;
    output.faces = surf.faces;
else
    output=phat;
end;


%figure;figure_patch(surfhat1,[0.74 0.71 0.61],0.5)

