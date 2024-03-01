function L2=L2norm(vector1, vector2)
%function L2=L2norm(vector1, vector2)
%
% It computes the L2 norm of vector1 (n by 3 matrix) if vector2 is not given. 
% If vector1 and vector2 is given, it computes the L2 distance between two.
%
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



if nargin<=1
    L2=sqrt(sum(vector1.^2,2));
else
        n=size(vector1,1);
        diff= vector1 - kron(ones(n,1), vector2);
        L2=sqrt(sum(diff.^2,2));
end;
    


