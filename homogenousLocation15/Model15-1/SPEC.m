%% Elastic Response Spectra 
% This is a function to generate elastic response specra including Displacement
% Spectrum, Pseudo Acceleration Spectrum and Pseudo Velocity Spectrum which
% are needed in "Response Spectrum Analysis" of Structures. In this function
% to solve "Equation of Motions" for different periods, Newmark Linear Method 
% has been used. 

%% @ Mostafa Tazarv, Carleton University, May 2011

%% SPEC Function Help:

% INPUTS:
% dt:     Time Interval (Sampling Time) of Record
% Ag:     Ground Motion Acceleration in g 
% zet:    Damping Ratio in percent (%); e.g. 5
% g:      Gravitational Constant; e.g. 9.81 m/s/s
% endp:   End Period of Spectra; e.g. 4 sec

% OUTPUTS:
% T:      Period of Structures (sec)
% Spa:    Elastic Pseudo Acceleration Spectrum
% Spv:    Elastic Pseudo Velocity Spectrum
% Sd:     Elastic Displacement Spectrum

function [T,Spa,Spv,Sd]=SPEC(dt,dT,Ag,zet,g,endp)
u=zeros(length(Ag),1);
v=zeros(length(Ag),1);
ac=zeros(length(Ag),1);
Ag=386.4*Ag/9.81;
Ag(end+1)=0;
T(2)=dT;
for j=2:round(endp/dT)                          % equation of motion(Newmark linear method)
    omega(j)=2*pi/T(j);      % Natural Frequency
    m=1;       
    k=(omega(j))^2*m;
    c=2*m*omega(j)*zet/100;
    K=k+3*c/dt+6*m/(dt)^2;
    a=6*m/dt+3*c;
    b=3*m+dt*c/2;    
  for i=1:length(u)-1
     u(1)=0;                      %initial conditions
     v(1)=0;
     ac(1)=0;    
     df=Ag(i+1)-Ag(i)+a*v(i)+b*ac(i);  % delta Force
     du=df/K;
     dv=3*du/dt-3*v(i)-dt*ac(i)/2;
     dac=6*(du-dt*v(i))/(dt)^2-3*ac(i);
     u(i+1)=u(i)+du;
     v(i+1)=v(i)+dv;
     ac(i+1)=ac(i)+dac;     
  end
    Sd(j)=max(abs(u));
    Spv(j)=Sd(j)*omega(j);
    Spa(j)=Sd(j)*(omega(j))^2/g;
    T(j+1)=T(j)+dT;
end
Ag(end)=[];
T(end)=[];
n=8;
T(1:n)=0.0;
Sd(1:n)=0; Spv(1:n)=0;Spa(1:n)=max(abs(Ag))/g;

 