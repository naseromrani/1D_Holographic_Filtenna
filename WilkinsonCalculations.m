
% // equi-ripple 3dB wideband two section wilkinson power devider design

clear 
close all
clc

% inputs

N=2;%input('Enter the number of sections "N": ');
f0=(10^9)*input('Enter the center frequency of band "f0" in GHz: ');
gama_max=(10^-2)*input('Enter the maximmum reflection coefficient "gama-max" in precent ( % ): ');
Z0=input('Enter the characteristic impedance(Z0) "port 1" in ohm: ');
ZL=input('Enter the impedance of load(ZL) "port 2 or 3" in ohm: ');

% parameters

t=0:0.01*pi:pi;
tm1=acos(1/cosh((1/N)*acosh(log(ZL/Z0)/(2*gama_max))));% // tm1 is in Radian
tm2=pi-tm1;
fm1=2*tm1*f0/pi;% // first frequency of band edge
fm2=2*tm2*f0/pi;% // second frequency of band edge

% calculations

phi=(pi/2)*(1-(1/sqrt(2))*((fm2-fm1)/(fm2+fm1)));
x=cos(t)./cos(tm1);
T2_x=2.*(x.^2)-1;
T2_tm=(1/gama_max)*abs((ZL-Z0)/(ZL+Z0));
A=((ZL-Z0)/(ZL+Z0))*(1/T2_tm);
gama=A.*(exp(-i*N.*t)).*T2_x;%#ok
g0=A/(2*((cos(tm1))^2));
g1=A*((tan(tm1))^2);
g2=g0;
Z1=Z0*exp(2*g0);
Z2=Z1*exp(2*g1);
R1=2*Z1*Z2/sqrt((Z1+Z2)*(Z2-Z1*((cot(phi))^2)));
R2=2*R1*(Z1+Z2)/(R1*(Z1+Z2)-2*Z2);

% plot and display

% plot(2*t*f0*(10^-9)/pi,abs(gama),'color','black','linewidth',1.5)
% title('total reflection coefficient','fontsize',15)
% xlabel('frequency(GHz)')
% ylabel('|gama|')
% grid on
display([R1 R2 Z1 Z2],'[R1 R2 Z1 Z2]')


