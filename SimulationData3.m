function [xg,yg,zg,xm,ym,zm,Len,fg,h,Theta,Gama,Fai,Qchange,noiseA,noiseM,A]=SimulationData3(noiseA,noiseM,fg,h,Beta,Qinit)
Beta=Beta/180*pi;
Len=100;
xg=zeros(1,Len);yg=xg;zg=-fg*ones(1,Len);
xm=h*cos(Beta)*ones(1,Len);ym=xg;zm=-h*sin(Beta)*ones(1,Len);
Theta=zeros(1,Len);Gama=zeros(1,Len);Fai=zeros(1,Len);Qchange=zeros(4,Len);
A=unifrnd(-1,1,3,3);A=orth(A);
Xekf=ConversionCtoQ(A);
if ~isempty(Qinit)
    Xekf=Qinit;
    q0=Qinit(4);q1=Qinit(1);q2=Qinit(2);q3=Qinit(3);
    A=[q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q3*q0) 2*(q1*q3-q2*q0);...
       2*(q1*q2-q3*q0) q0^2+q2^2-q1^2-q3^2 2*(q2*q3+q1*q0);...
       2*(q1*q3+q2*q0) 2*(q2*q3-q1*q0) q0^2+q3^2-q2^2-q1^2];
end

for k=1:Len
    tmp=A*[xg(k);yg(k);zg(k)];
    xg(k)=tmp(1);yg(k)=tmp(2);zg(k)=tmp(3);
    tmp=A*[xm(k);ym(k);zm(k)];
    xm(k)=tmp(1);ym(k)=tmp(2);zm(k)=tmp(3);
    Qchange(:,k)=Xekf;
    [Theta,Gama,Fai]=Attitude(A',Theta,Gama,Fai,k);
end
Theta=Theta'./pi.*180;Gama=Gama'./pi.*180;Fai=Fai'./pi.*180;

xg=xg+normrnd(0,noiseA,1,Len);yg=yg+normrnd(0,noiseA,1,Len);zg=zg+normrnd(0,noiseA,1,Len);
xm=xm+normrnd(0,noiseM,1,Len);ym=ym+normrnd(0,noiseM,1,Len);zm=zm+normrnd(0,noiseM,1,Len);

function [Theta,Gama,Fai]=Attitude(C,Theta,Gama,Fai,Index)
Theta(Index)=-asin(C(3,1));   
Gama(Index)=atan2(C(3,2),C(3,3));
Fai(Index)=atan2(C(2,1),C(1,1));

function Q=ConversionCtoQ(C)
C=C./det(C);
q0=0.5*sqrt(1+C(1,1)+C(2,2)+C(3,3));
q3=(C(2,1)-C(1,2))/4/q0;
q2=(C(1,3)-C(3,1))/4/q0;
q1=(C(3,2)-C(2,3))/4/q0;
Q=[q1;q2;q3;q0];
Q=real(Q);
Q=Q./norm(Q);
