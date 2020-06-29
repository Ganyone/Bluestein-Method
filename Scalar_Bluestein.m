%%%%%%Scalar diffraction computation method using Bluestein method
%%%%%%The work is under review by Light: Science and Applications;
%%%%%%Please cite properly if you use the codes for any use
%%%%%%unit: um

function [gout,pixelout] = Scalar_Bluestein(gin,mxin,myin,pixelin,d,xstart,xend,ystart,yend,mxout,myout)
%%%%%% Defination of input variable:
%%%%%% gin----Complex amplitude of the incident light beam
%%%%%% mxin,myin-----resolution of the input plane in transverse and longitudinal directions
%%%%%% pixelin----pixel size of the incident light beam
%%%%%% d----diffraction distance
%%%%%% xstart,xend----computation range in transverse direction
%%%%%% ystart,yend----computation range in longitudinal direction
%%%%%% mxout,myout----resolution of the output plane in transverse and longitudinal directions
%%%%%% Defination of output variable:
%%%%%% gout----Complex amplitude of the outgoing light beam
%%%%%% pixelout----pixel size of the outgoing light beam

global lamda k

L0=(myin-1)*pixelin;                                                        % dimension of the diffraction plane
x0=linspace(-L0/2,L0/2,mxin);                                               % transverse coordinates (x) in the diffraction plane
y0=linspace(-L0/2,L0/2,myin);                                               % longitudinal coordinates (y) in the diffraction plane
[x0,y0]=meshgrid(x0,y0);

pixelout=(xend-xstart)/(mxout-1);                                           % pixel size of the outgoing light beam
x1=linspace(xstart,xend,mxout);                                             % transverse coordinates (x) in the imaging plane
y1=linspace(ystart,yend,myout);                                             % longitudinal coordinates (y) in the imaging plane
[x1,y1]=meshgrid(x1,y1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          % calculating scalar diffraction below
F0=exp(j*k*d)/(j*lamda*d).*exp(j*k/2/d*(x1.^2+y1.^2));                      % assign exp(ikd)/(i¦Ëd)exp[ik(x2+y2) /2d]; Equation 4 in the paper
F=exp(j*k/2/d*(x0.^2+y0.^2));                                               % assign exp[ik (x02+y02) /2d]; Equation 5 in the paper
gout=gin.*F;
% gout=F0.*fftshift(fft2(fftshift(gout)));                                  % using FFT to calculate the complex amplitude of the outgoing light beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       % using Bluestein method to calculate the complex amplitude of the outgoing light beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            % one-dimensional FFT in one direction
fs = lamda*d/pixelin;                                                       % dimension of the imaging plane
fy1=ystart+fs/2;
fy2=yend+fs/2;
a = -exp(j*2*pi*fy1/fs);                                                    % define A according to Equation 13
w = exp(-j*2*pi*(fy2-fy1)/(myout*fs));                                      % define W according to Equation 14
gout = czt(gout,myout,w,a);%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
lly=linspace(0,myout-1,myout);lly=lly./myout.*(fy2-fy1)+fy1;
Mshifty=floor(-myin/2);
Mshift=repmat(exp(-1i.*2*pi.*lly.*Mshifty/fs),[mxin 1]);                    % define Pshift according to Equation 15
gout=gout.'.*Mshift;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                            % one-dimensional FFT in the other direction
fx1=xstart+fs/2;
fx2=xend+fs/2;
a = -exp(j*2*pi*fx1/fs);                                                    % define A according to Equation 13
w = exp(-j*2*pi*(fx2-fx1)/(mxout*fs));                                      % define W according to Equation 14
gout = czt(gout,mxout,w,a);%%%%%%%%%%%%%%%%%%%%%%%%%%%%
llx=linspace(0,mxout-1,mxout);llx=llx./mxout.*(fx2-fx1)+fx1;
Mshiftx=floor(-mxin/2);
Mshift=repmat(exp(-1i.*2*pi.*llx.*Mshiftx/fs),[myout 1]);                   % define Pshift according to Equation 15
gout=gout.'.*Mshift;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gout=F0.*gout;                                                              % obtain the complex amplitude of the outgoing light beam
end