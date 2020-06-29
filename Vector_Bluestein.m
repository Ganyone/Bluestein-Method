%%%%%%Vector diffraction computation method using Bluestein method
%%%%%%The work is under review by Light: Science and Applications;
%%%%%%Please cite properly if you use the codes for any use
%%%%%%unit: um

function [Ex Ey Ez] = Vector_Bluestein(E,Min,polar,xstart,xend,ystart,yend,z0,Moutx,Mouty)
%%%%%% Definition of input variable:
%%%%%% E----electric field of the incident light beam
%%%%%% Min-----resolution of the input plane in transverse and longitudinal directions
%%%%%% polar----polarization
%%%%%% xstart,xend----computation range in transverse direction
%%%%%% ystart,yend----computation range in longitudinal direction
%%%%%% z0----position shift along the optical axis
%%%%%% Mxout,Myout----resolution of the output plane in transverse and longitudinal directions
%%%%%% Defination of output variable:
%%%%%% Ex Ey Ez----vector components in three orthogonal directions
global lamda k n1 NA fo
R=fo.*NA./n1;% Objective back aperture (radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Debye diffraction calculated using Blustein method
N=(Min-1)/2; 
m=linspace(-Min/2,Min/2,Min);                                               % array in X dimension 
n=linspace(-Min/2,Min/2,Min);                                               % array in Y dimension
[m n] = meshgrid(m,n);                                                      % convert to 2D mesh
th=asin(NA.*sqrt(m.^2 + n.^2)./(N.*n1));                                    % angle of convergance (theta)
thh=th;
th(thh>asin(NA./n1))=0;                                                     % remove parts outside numerical aperture
E(thh>asin(NA./n1))=0;                                                      % remove parts outside numerical aperture
phi = atan2 (n,m);                                                          % define phi azimuthal angle (0-2Pi)
phi(phi<0) = phi(phi<0)+2.*pi;
z=z0;
switch polar
    case('x')
    % % % % % % % X-polarization
     Ex = E./sqrt(cos(th)).*(1+(cos(th)-1).*cos(phi).^2).*exp (i.*k.*n1.*z.*cos (th));
     Ey = E./sqrt(cos(th)).*(cos (phi).*sin (phi).*(cos (th)-1)).*exp (i.*k.*n1.*z.*cos (th));
     Ez = E./sqrt(cos(th)).*(cos (phi).*sin (th)).*exp (i.*k.*n1.*z.*cos (th));
    %  
    case('y')
     % % % % % % % Y-polarization
     Ex = E./sqrt(cos(th)).*((cos(th)-1).*cos(phi).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));
     Ey = E./sqrt(cos(th)).*(1+(cos(th)-1).*sin(phi).^2).*exp (i.*k.*n1.*z.*cos (th));
     Ez = E./sqrt(cos(th)).*(sin (phi).*sin (th)).*exp (i.*k.*n1.*z.*cos (th));

     case('rc')

     % % % % % % % Right circular
     Ex = (1./sqrt(2)).*E./sqrt(cos(th)).*(1+(cos(th)-1).*cos(phi).^2+i.*(cos(th)-1).*cos(phi).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));
     Ey = (1./sqrt(2)).*E./sqrt(cos(th)).*((cos(th)-1).*cos(phi).*sin(phi)+i.*(1+(cos(th)-1).*sin(phi).^2)).*exp (i.*k.*n1.*z.*cos (th));
     Ez = (1./sqrt(2)).*E./sqrt(cos(th)).*(sin(th).*cos(phi)+i.*sin(th).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th)); 

    case('lc') 
     % % % % % % % Left circular
     Ex = (1./sqrt(2)).*E./sqrt(cos(th)).*(i.*(1+(cos(th)-1).*cos(phi).^2)+(cos(th)-1).*cos(phi).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));
     Ey = (1./sqrt(2)).*E./sqrt(cos(th)).*(i.*(cos(th)-1).*cos(phi).*sin(phi)+1+(cos(th)-1).*sin(phi).^2).*exp (i.*k.*n1.*z.*cos (th));
     Ez = (1./sqrt(2)).*E./sqrt(cos(th)).*(i.*sin(th).*cos(phi)+sin(th).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th)); 

     case('ra')
     % % % % % % Radical polarization
     Ex = E./sqrt(cos(th)).*(cos(phi).*(1+(cos(th)-1).*cos(phi).^2)+sin(phi).*(cos(th)-1).*cos(phi).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));
     Ey = E./sqrt(cos(th)).*(cos(phi).*(cos(th)-1).*cos(phi).*sin(phi)+sin(phi).*(1+(cos(th)-1).*sin(phi).^2)).*exp (i.*k.*n1.*z.*cos (th));
     Ez = E./sqrt(cos(th)).*(cos(phi).*sin(th).*cos(phi)+sin(phi).*sin(th).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));

     case('az') 
     % % % % % % % Azimuthal polarization
     Ex = E./sqrt(cos(th)).*(-sin(phi).*(1+(cos(th)-1).*cos(phi).^2)+cos(phi).*(cos(th)-1).*cos(phi).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th));
     Ey = E./sqrt(cos(th)).*(-sin(phi).*(cos(th)-1).*cos(phi).*sin(phi)+cos(phi).*(1+(cos(th)-1).*sin(phi).^2)).*exp (i.*k.*n1.*z.*cos (th));
     Ez = E./sqrt(cos(th)).*(-sin(phi).*sin(th).*cos(phi)+cos(phi).*sin(th).*sin(phi)).*exp (i.*k.*n1.*z.*cos (th)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Ex (real (thh)>asin (NA./n1))=0;                                          
 Ey (real (thh)>asin (NA./n1))=0;
 Ez (real (thh)>asin (NA./n1))=0;

 where=isnan (Ex);
 Ex (where)=1;                                                              % set all Not A Number values to 1
 where=isnan (Ey);
 Ey (where)=1;                                                              % set all Not A Number values to 1
 where=isnan (Ez);
 Ez (where)=1;                                                              % set all Not A Number values to 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% using Bluestein method to replace FFT
fs = fo.*lamda./n1.*(Min-1)./(2.*R);                                        % dimension of the imaging plane
Cystart =fs./2+ystart;                                                      % start position in y direction
Cyend =fs./2+yend;                                                          % end position in y direction
Cxstart = fs./2+xstart;                                                     % start position in x direction
Cxend = fs./2+xend;                                                         % end position in x direction

ay = -exp(j*2*pi*Cystart/fs);                                               % define A (y) according to Equation 13
wy = exp(-j*2*pi*(Cyend-Cystart)/(Mouty*fs));                               % define W (y) according to Equation 14
ax = -exp(j*2*pi*Cxstart/fs);                                               % define A (x) according to Equation 13
wx = exp(-j*2*pi*(Cxend-Cxstart)/(Moutx*fs));                               % define W (x) according to Equation 14


EHold = czt(Ex,Mouty,wy,ay);
lly=linspace(0,Mouty-1,Mouty);lly=lly./Mouty.*(Cyend-Cystart)+Cystart;
Mshifty=floor(-Min/2);                                                      % define Pshift (y) according to Equation 15
Mshifty=repmat(exp(-1i.*2*pi.*lly.*Mshifty/fs),[Min 1]);
EHold=EHold.'.*Mshifty;
EHold = czt(EHold,Moutx,wx,ax);
llx=linspace(0,Moutx-1,Moutx);llx=llx./Moutx.*(Cxend-Cxstart)+Cxstart;
Mshiftx=floor(-Min/2);                                                      % define Pshift (x) according to Equation 15
Mshiftx=repmat(exp(-1i.*2*pi.*llx.*Mshiftx/fs),[Mouty 1]);
Ex=EHold.'.*Mshiftx;

EHold = czt(Ey,Mouty,wy,ay);
EHold=EHold.'.*Mshifty;
EHold = czt(EHold,Moutx,wx,ax);
Ey=EHold.'.*Mshiftx;

EHold = czt(Ez,Mouty,wy,ay);
EHold=EHold.'.*Mshifty;
EHold = czt(EHold,Moutx,wx,ax);
Ez=EHold.'.*Mshiftx;

clear EHold;                    
end