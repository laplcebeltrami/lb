 function stat = stat_Fcorrected(stat, sigma, surf);
% pvalue = stat_Fcorrected(stat, sigma, surf);
%
% stat_Fcorrected returns the corrected p-value of F-statistic map with
% given smoothing kernel sigma.
%
% stat    : F-statistics. 
% sigma : smoothing kernel size
% surf    : surface mesh
%
% This function is modified from Keith J. Worsly's original fMRI-STAT
% package using the random field theory (Worsley et al. 
% 1996, Human Brain Mapping, 4:58-73).
% The random field theory threshold is based on the assumption the data is isotropic, 
% i.e. equal fwhm in each direction. If the data is not isotropic, then replacing
% FWHM by the geometric mean of the three x,y,z fwhm's is a very good
% approximation.
%
%
% (C) 2012 Moo K. Chung
% University of Wisconsin-Madison
% mkchung@wisc.edu

t=stat.tstat;
a=stat.df(1);
b=stat.df(2);


lambda=1/(2*sigma^2);
%FWHM=2*sqrt(log(4))*sigma;

mu0=2;  % 2 for 1 surface.
area=GETarea(surf);
mu2=sum(area.faces)/2;

% T-stat
%rho0 = 1 - tcdf(t,df);
%rho2 = 1/FWHM^2*4*log(2)/(2*pi)^(3/2)*gamma(df/2)/((df-1)/2)^(1/2)/gamma((df-1)/2)*...
%    t.*(1+t.^2/(df-1)).^(-(df-2)/2);

rho0 = 1 - fcdf(t,a,b);
%rho2 = 1/FWHM^2*4*log(2)/(2*pi)*exp(gammaln((a+b-2)/2)-gammaln(a/2)-gammaln(b/2))...
%    *(a*t/b).^((a-2)/2).*(1+a*t/b).^(-(a+b-2)/2).*((b-1)*a*t/b-(a-1));

rho2 = 1/(4*sigma^2*pi)*gamma((a+b-2)/2)/ gamma(a/2)/gamma(b/2)...
    *(a*t/b).^((a-2)/2).*(1+a*t/b).^(-(a+b-2)/2).*((b-1)*a*t/b-(a-1));

p=mu0*rho0+rho2*mu2;

%ind=p>0.1;
%p(ind)=0.1;
stat.pcorrected=p; 

%---------------------------------------------------------
function area =GETarea(surf);

coord=surf.vertices;
tri=surf.faces;
n_faces=size(surf.faces,1);
n_vertices=size(surf.vertices,1);

area_f=zeros(n_faces,1);

for i=1:n_faces
    
    p = coord(tri(i,:),:);
    p1 = p(1,:)-p(3,:);
    p2 = p(2,:)-p(3,:);
    q = cross(p1,p2);
    area_f(i)= sqrt(q*q')/2;
end;

area.faces=area_f;
%--------------------
% area at vertex
area_v=zeros(n_vertices,1);


for i=1:size(area_f,1)
    area_v(tri(i,1))=  area_v(tri(i,1)) + area_f(i);
    area_v(tri(i,2))=  area_v(tri(i,2)) + area_f(i);
    area_v(tri(i,3))=  area_v(tri(i,3)) + area_f(i);
end;

% average area of the 1st neighbor.
area.vertices =area_v./3; %degree;









