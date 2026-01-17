function [Q,C,Xlkf,Plkf,Q_NQDM]=LKF3(Gx,Gy,Gz,xg,yg,zg,xm,ym,zm,fg,h,Beta,SampleTime,Xlkf,Plkf,noiseG,noiseA,noiseM,AdaptiveRegulation,Q_NQDM,SAD,useRegressionModel,R)
fg2=sqrt(xg^2+yg^2+zg^2);h2=sqrt(xm^2+ym^2+zm^2);xg=xg/fg2;yg=yg/fg2;zg=zg/fg2;xm=xm/h2;ym=ym/h2;zm=zm/h2;
dThetaX=Gx*SampleTime;dThetaY=Gy*SampleTime;dThetaZ=Gz*SampleTime;               
Omiga=[0,dThetaZ,-dThetaY,dThetaX;-dThetaZ,0,dThetaX,dThetaY;dThetaY,-dThetaX,0,dThetaZ;-dThetaX,-dThetaY,-dThetaZ,0];
dTheta0=sqrt(dThetaX^2+dThetaY^2+dThetaZ^2);if dTheta0==0, dTheta0=10^-9;end
TrM=cos(dTheta0/2)*eye(4)+sin(dTheta0/2)/dTheta0*Omiga;
Xekf_=TrM*Xlkf;q0=Xekf_(4);q1=Xekf_(1);q2=Xekf_(2);q3=Xekf_(3);
Sigm=[0 -q3 q2;q3 0 -q1;-q2 q1 0];Sigm=Sigm+[q0 0 0;0 q0 0;0 0 q0];Sigm=[Sigm;-[q1 q2 q3]];
linF=eye(4);
a1inv=1/AdaptiveRegulation;if a1inv==0, a1inv=10^-9;end
sigQ=AdaptiveRegulation*(0.5*SampleTime)^2*Sigm*[noiseG 0 0;0 noiseG 0;0 0 noiseG]*Sigm'; 
if SAD==1
    [Qm,~,sigR,~]=QUEST3(xg,yg,zg,xm,ym,zm,fg,h,Beta,noiseA,noiseM,1,R);
elseif SAD==2
    [Qm,~,sigR,~]=QUEST3(xg,yg,zg,xm,ym,zm,fg,h,Beta,noiseA,noiseM,0,R);
elseif SAD==3
    [Qm,~,sigR,~]=AQUA3(xg,yg,zg,xm,ym,zm,fg,h,Beta,noiseA,noiseM,useRegressionModel);
elseif SAD==4
    [Qm,~,sigR,~]=NQDM3(xg,yg,zg,xm,ym,zm,fg,h,Beta,noiseA,noiseM,Q_NQDM,useRegressionModel);Q_NQDM=Qm;
else 
    [Qm,~,~,~]=QUEST(xg,yg,zg,xm,ym,zm,fg,h,Beta,noiseA,noiseM,1);
    sigR=[[noiseA*eye(3) zeros(3,1)];[zeros(1,3) noiseM]];
end
sigR=a1inv*sigR;
[~,d]=max(abs(Qm));
if sign(Xekf_(d))~=sign(Qm(d))
    Qm=-Qm;
end

Pekf_=TrM*Plkf*TrM'+sigQ;
Kk=Pekf_*linF'/(linF*Pekf_*linF'+sigR(:,1:4));
Xlkf=Xekf_+Kk*(Qm-Xekf_);
Plkf=Pekf_-Kk*linF*Pekf_;

Xlkf=Xlkf./norm(Xlkf);
Q=Xlkf(1:4);q0=Q(4);q1=Q(1);q2=Q(2);q3=Q(3);
C=[q0^2+q1^2-q2^2-q3^2 2*(q1*q2-q3*q0) 2*(q1*q3+q2*q0);2*(q1*q2+q3*q0) q0^2+q2^2-q1^2-q3^2 2*(q2*q3-q1*q0);2*(q1*q3-q2*q0) 2*(q2*q3+q1*q0) q0^2+q3^2-q2^2-q1^2];
