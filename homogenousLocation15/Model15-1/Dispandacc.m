%% Obtain the displacement and accelerations from velocity data 

function [Disp_SW4, Acc_SW4]=Dispandacc(dt,Velocity)
Disp_SW4(1)=0;
Acc_SW4(1)=0;
N=length(Velocity);
for i=1:N-1
    Disp_SW4(i+1)=Disp_SW4(i)+Velocity(i)*dt;
    Acc_SW4(i+1)=(Velocity(i+1)-Velocity(i))/dt;
end
end

 