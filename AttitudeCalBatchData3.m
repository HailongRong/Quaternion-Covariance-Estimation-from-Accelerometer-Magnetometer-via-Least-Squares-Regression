function [Theta,Gama,Fai,Cchange,Qchange,Covairance,TimeSAD,myrou]=AttitudeCalBatchData3(Algorithm,Gx,Gy,Gz,xg,yg,zg,xm,ym,zm,AdaptiveRegulation,fg,h,Beta,noiseG,noiseA,noiseM,frequency,SAD,useRegressionModel)
SampleTime=1/frequency;
Beta=Beta/180*pi;
R=ComputeRegressionMatrix(noiseA,noiseM,Beta);
LenData=length(xg);
Theta=zeros(1,LenData);Gama=zeros(1,LenData);Fai=zeros(1,LenData);Cchange=cell(1,LenData);Qchange=zeros(4,LenData);
noiseG=noiseG^2;noiseA=noiseA^2/fg^2;noiseM=noiseM^2/h^2;
End1=100;
[Q,C,~,~]=QUEST3(mean(xg(1:End1)),mean(yg(1:End1)),mean(zg(1:End1)),mean(xm(1:End1)),mean(ym(1:End1)),mean(zm(1:End1)),fg,h,Beta,noiseA,noiseM,1,R);
Plkf=10^-9*eye(4);
[Theta,Gama,Fai]=Attitude(C,Theta,Gama,Fai,1);Cchange{1}=C;Qchange(:,1)=Q;Xlkf=Q;Q_NQDM=Q;
TimeSAD=zeros(1,LenData-1);
Covairance=cell(1,LenData-1);
myrou=zeros(2,LenData-1);

for k=2:LenData
    switch Algorithm
        case 'LKF3'
            [Q,C,Xlkf,Plkf,Q_NQDM]=LKF3(Gx(k),Gy(k),Gz(k),xg(k),yg(k),zg(k),xm(k),ym(k),zm(k),fg,h,Beta,SampleTime,Xlkf,Plkf,noiseG,noiseA,noiseM,AdaptiveRegulation(k),Q_NQDM,SAD,useRegressionModel,R);
        case 'EKF3_flops'
            [Q,C,Xlkf,Plkf,myrou]=EKF3_flops(Gx(k),Gy(k),Gz(k),xg(k),yg(k),zg(k),xm(k),ym(k),zm(k),fg,h,Beta,SampleTime,Xlkf,Plkf,noiseG,noiseA,noiseM,AdaptiveRegulation(k),[],k,1);
        case 'QUEST3'
            [Q,C,PQUEST,TimeQUEST]=QUEST3(xg(k),yg(k),zg(k),xm(k),ym(k),zm(k),fg,h,Beta,noiseA,noiseM,useRegressionModel,R);
            Covairance{k-1}=PQUEST;
            TimeSAD(k-1)=TimeQUEST;
            myrou(:,k-1)=[trace(PQUEST);det(PQUEST)];
        case 'AQUA3'
            [Q,C,PAQUA,TimeAQUA]=AQUA3(xg(k),yg(k),zg(k),xm(k),ym(k),zm(k),fg,h,Beta,noiseA,noiseM);
            Covairance{k-1}=PAQUA;
            TimeSAD(k-1)=TimeAQUA;
            myrou(:,k-1)=[trace(PAQUA);det(PAQUA)];
        case 'NQDM3'
            [Q,C,PNQDM,TimeNQDM]=NQDM3(xg(k),yg(k),zg(k),xm(k),ym(k),zm(k),fg,h,Beta,noiseA,noiseM,Q);
            Covairance{k-1}=PNQDM;
            TimeSAD(k-1)=TimeNQDM;
            myrou(:,k-1)=[trace(PNQDM);det(PNQDM)];
        otherwise
            Theta=[];Gama=[];Fai=[];Cchange=[];Qchange=[];return;
    end
    [Theta,Gama,Fai]=Attitude(C,Theta,Gama,Fai,k);Cchange{k}=C;Qchange(:,k)=Q;
end
Theta=Theta'*180/pi;Gama=Gama'*180/pi;Fai=Fai'*180/pi;

function [Theta,Gama,Fai]=Attitude(C,Theta,Gama,Fai,Index)
Theta(Index)=-asin(C(3,1));   
Gama(Index)=atan2(C(3,2),C(3,3));
Fai(Index)=atan2(C(2,1),C(1,1));

function R=ComputeRegressionMatrix(noiseA,noiseM,Beta)
Beta=Beta/pi*180/90;
c12_3=[0.9998   -0.0021    2.6135    0.0330   -3.8506   -0.0883   41.5300   -0.0996  -13.2189    0.2003  -39.6314    0.2834  -17.3110    0.1005   31.0476   -0.2130   82.3518 -0.5307  123.4874];
c12_6=[0.0000    0.0000    0.0067   -0.0000   -0.0089   -0.0000    0.1052    0.0006   -0.0342   -0.0004   -0.0999   -0.0009   -0.0418   -0.0005    0.0815    0.0005    0.2114 0.0017    0.3150];
c78_3=[0.7529    0.0019    2.5899   -0.0062   -3.6488   -0.0176   41.0185    0.0754  -13.1773   -0.0005  -39.0979   -0.0652  -16.7953   -0.0682   31.1546   -0.0265   81.8938 0.0322  122.4874];
c78_6=[-0.0006    0.0000    0.0067   -0.0001   -0.0093    0.0004    0.1061    0.0002   -0.0343   -0.0008   -0.1010   -0.0010   -0.0428   -0.0003    0.0814    0.0009    0.2125 0.0020    0.3171];
c69_3=[-1.0003    0.0006   -2.5905   -0.0068    3.6456    0.0060  -41.0024    0.0657   13.1909   -0.0513   39.0561   -0.1140   16.7309   -0.0559  -31.1789    0.0845  -81.8373 0.2492 -122.3418];
c69_6=[0.0026   -0.0000   -0.0067   -0.0000    0.0092    0.0001   -0.1058   -0.0004    0.0342    0.0001    0.1006    0.0004    0.0425    0.0003   -0.0813   -0.0001   -0.2121 -0.0007   -0.3164];
c19_6=[0.0000    0.0041   -0.0000    0.0024    0.0002    0.0068   -0.0003    0.0026   -0.0000   -0.0027    0.0003   -0.0014    0.0003    0.0048    0.0001    0.0125   -0.0002 0.0191   -0.0006];
c13_3=0.2477;
c13_6=0;
c14_6=0.0026;

A=[1;Beta;Beta.^2;Beta.^3;Beta.^4;Beta.^5;Beta.^6;Beta.^7;Beta.^8;Beta.^9;Beta.^10;Beta.^11;Beta.^12;Beta.^13;Beta.^14;Beta.^15;Beta.^16;Beta.^17;Beta.^18];
% Estimating R12-----------------------------------------------------------
r3=c12_3*A;r6=c12_6*A; 
R12=r3*noiseM.^2+r6*noiseA.^2;
R34=R12;R55=-R12;R1010=-R12;
% Estimating R78-----------------------------------------------------------
r3=c78_3*A;r6=c78_6*A;
R78=r3*noiseM.^2+r6*noiseA.^2;
% Estimating R69-----------------------------------------------------------
r3=c69_3*A;r6=c69_6*A;
R69=r3*noiseM.^2+r6*noiseA.^2;
% Estimating R19-----------------------------------------------------------
r6=c19_6*A;
R19=r6*noiseA.^2;
R26=R19;R710=R19;R810=R19;R39=-R19;R46=-R19;R57=-R19;R58=-R19;
% Estimating R13-----------------------------------------------------------
R13=c13_3*noiseM.^2+c13_6*noiseA.^2;
R24=R13;R66=-R13;R99=-R13;
% Estimating R14-----------------------------------------------------------
R14=c14_6*noiseA.^2;
R23=R14;R77=-R14;R88=-R14;
% Constructing the regression matrix---------------------------------------
R=zeros(10,10);
R(1,2)=R12;R(3,4)=R34;R(5,5)=R55;R(10,10)=R1010;
R(7,8)=R78;R(6,9)=R69;
R(1,3)=R13;R(2,4)=R24;R(6,6)=R66;R(9,9)=R99;
R(1,4)=R14;R(2,3)=R23;R(7,7)=R77;R(8,8)=R88;
R(1,9)=R19;R(2,6)=R26;R(7,10)=R710;R(8,10)=R810;R(3,9)=R39;R(4,6)=R46;R(5,7)=R57;R(5,8)=R58;