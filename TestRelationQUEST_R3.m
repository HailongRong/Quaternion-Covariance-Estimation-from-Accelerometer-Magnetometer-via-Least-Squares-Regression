Len2=100;Len=1000;
VSRMest11=zeros(Len2,Len-1);VSRMest22=VSRMest11;VSRMest33=VSRMest11;VSRMest00=VSRMest11;VSRMest12=VSRMest11;VSRMest13=VSRMest11;VSRMest10=VSRMest11;VSRMest23=VSRMest11;VSRMest20=VSRMest11;VSRMest30=VSRMest11;
VSRMvar11=zeros(Len2,1);VSRMvar22=VSRMvar11;VSRMvar33=VSRMvar11;VSRMvar00=VSRMvar11;VSRMvar12=VSRMvar11;VSRMvar13=VSRMvar11;VSRMvar10=VSRMvar11;VSRMvar23=VSRMvar11;VSRMvar20=VSRMvar11;VSRMvar30=VSRMvar11;
QUESTest11=VSRMest11;QUESTest22=VSRMest11;QUESTest33=VSRMest11;QUESTest00=VSRMest11;QUESTest12=VSRMest11;QUESTest13=VSRMest11;QUESTest10=VSRMest11;QUESTest23=VSRMest11;QUESTest20=VSRMest11;QUESTest30=VSRMest11;
QUESTvar11=VSRMvar11;QUESTvar22=VSRMvar11;QUESTvar33=VSRMvar11;QUESTvar00=VSRMvar11;QUESTvar12=VSRMvar11;QUESTvar13=VSRMvar11;QUESTvar10=VSRMvar11;QUESTvar23=VSRMvar11;QUESTvar20=VSRMvar11;QUESTvar30=VSRMvar11;
AQUAest11=VSRMest11;AQUAest22=VSRMest11;AQUAest33=VSRMest11;AQUAest00=VSRMest11;AQUAest12=VSRMest11;AQUAest13=VSRMest11;AQUAest10=VSRMest11;AQUAest23=VSRMest11;AQUAest20=VSRMest11;AQUAest30=VSRMest11;
AQUAvar11=VSRMvar11;AQUAvar22=VSRMvar11;AQUAvar33=VSRMvar11;AQUAvar00=VSRMvar11;AQUAvar12=VSRMvar11;AQUAvar13=VSRMvar11;AQUAvar10=VSRMvar11;AQUAvar23=VSRMvar11;AQUAvar20=VSRMvar11;AQUAvar30=VSRMvar11;
NQDMest11=VSRMest11;NQDMest22=VSRMest11;NQDMest33=VSRMest11;NQDMest00=VSRMest11;NQDMest12=VSRMest11;NQDMest13=VSRMest11;NQDMest10=VSRMest11;NQDMest23=VSRMest11;NQDMest20=VSRMest11;NQDMest30=VSRMest11;
NQDMvar11=VSRMvar11;NQDMvar22=VSRMvar11;NQDMvar33=VSRMvar11;NQDMvar00=VSRMvar11;NQDMvar12=VSRMvar11;NQDMvar13=VSRMvar11;NQDMvar10=VSRMvar11;NQDMvar23=VSRMvar11;NQDMvar20=VSRMvar11;NQDMvar30=VSRMvar11;
LSMest11=VSRMest11;LSMest22=VSRMest11;LSMest33=VSRMest11;LSMest00=VSRMest11;LSMest12=VSRMest11;LSMest13=VSRMest11;LSMest10=VSRMest11;LSMest23=VSRMest11;LSMest20=VSRMest11;LSMest30=VSRMest11;
LSMvar11=VSRMvar11;LSMvar22=VSRMvar11;LSMvar33=VSRMvar11;LSMvar00=VSRMvar11;LSMvar12=VSRMvar11;LSMvar13=VSRMvar11;LSMvar10=VSRMvar11;LSMvar23=VSRMvar11;LSMvar20=VSRMvar11;LSMvar30=VSRMvar11;
VSRMerror=VSRMvar11;QUESTerror=VSRMvar11;AQUAerror=VSRMvar11;NQDMerror=VSRMvar11;LSMerror=VSRMvar11;
dVSRMerror=VSRMvar11;dQUESTerror=VSRMvar11;dAQUAerror=VSRMvar11;dNQDMerror=VSRMvar11;dLSMerror=VSRMvar11;
C11QUESTin=zeros(1,Len);C12QUESTin=C11QUESTin;C13QUESTin=C11QUESTin;
C21QUESTin=C11QUESTin;C22QUESTin=C11QUESTin;C23QUESTin=C11QUESTin;
C31QUESTin=C11QUESTin;C32QUESTin=C11QUESTin;C33QUESTin=C11QUESTin;
Noise=zeros(Len2,3);

traceLSM2=zeros(Len2,Len-1);detLSM2=traceLSM2;RMSerrorLSM2=traceLSM2;
traceQUEST2=traceLSM2;detQUEST2=traceLSM2;RMSerrorQUEST2=traceLSM2;
traceAQUA2=traceLSM2;detAQUA2=traceLSM2;RMSerrorAQUA2=traceLSM2;
traceNQDM2=traceLSM2;detNQDM2=traceLSM2;RMSerrorNQDM2=traceLSM2;

RecTimeVSRM=zeros(Len2,Len-1);RecTimeQUEST=RecTimeVSRM;RecTimeAQUA=RecTimeVSRM;RecTimeNQDM=RecTimeVSRM;RecTimeLSM=RecTimeVSRM;

TRecQ1QUEST=zeros(Len2,Len);TRecQ2QUEST=TRecQ1QUEST;TRecQ3QUEST=TRecQ1QUEST;TRecQ0QUEST=TRecQ1QUEST;
TRecQ1VSRM=zeros(Len2,Len);TRecQ2VSRM=TRecQ1QUEST;TRecQ3VSRM=TRecQ1QUEST;TRecQ0VSRM=TRecQ1QUEST;
TRecQ1LSM=zeros(Len2,Len);TRecQ2LSM=TRecQ1QUEST;TRecQ3LSM=TRecQ1QUEST;TRecQ0LSM=TRecQ1QUEST;
TRecQ1AQUA=TRecQ1QUEST;TRecQ2AQUA=TRecQ1QUEST;TRecQ3AQUA=TRecQ1QUEST;TRecQ0AQUA=TRecQ1QUEST;
TRecQ1NQDM=TRecQ1QUEST;TRecQ2NQDM=TRecQ1QUEST;TRecQ3NQDM=TRecQ1QUEST;TRecQ0NQDM=TRecQ1QUEST;
noiseG=0.002;noiseA=0.008;noiseM=0.0003;frequency=1/0.0033;h=0.4885;fg=9.7387;Beta=-49.1213; 
AdaptiveRegulation=ones(1,Len);Xinit=[0;1;0;0];UseMag=1;
Gx=zeros(1,Len);Gy=Gx;Gz=Gx;noiseG=0;

for Num=1:Len2
    noiseA=unifrnd(0.001,0.08);noiseM=unifrnd(0.0001,0.004);Beta=unifrnd(-80,80);
    Noise(Num,:)=[noiseA,noiseM,Beta];
    [xg,yg,zg,xm,ym,zm,Len,fg,h,Theta,Gama,Fai,Qchange,noiseA,noiseM,A]=SimulationData3(noiseA,noiseM,fg,h,Beta,[]);
    %%LSRM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [ThetaLSM,GamaLSM,FaiLSM,CchangeLSM,QchangeLSM,myrou,TimeLSM,myrouLSM] = AttitudeCalBatchData3('QUEST3',Gx,Gy,Gz,xg,yg,zg,xm,ym,zm,AdaptiveRegulation,fg,h,Beta,noiseG,noiseA,noiseM,frequency,1,1);
    RecTimeLSM(Num,:)=TimeLSM;
    [~,Ind1]=max(abs(QchangeLSM(:,1)));
    for k=2:Len
        if QchangeLSM(Ind1,1)*QchangeLSM(Ind1,k)<0
            QchangeLSM(:,k)=-QchangeLSM(:,k);
        end
    end
    TRecQ1LSM(Num,:)=QchangeLSM(1,:);TRecQ2LSM(Num,:)=QchangeLSM(2,:);TRecQ3LSM(Num,:)=QchangeLSM(3,:);TRecQ0LSM(Num,:)=QchangeLSM(4,:);
    tmp=cov(QchangeLSM(:,2:Len)');
    traceLSM2(Num,:)=log((myrouLSM(1,:)-trace(tmp)).^2);
    detLSM2(Num,:)=log((myrouLSM(2,:)).^2);
    tmp2=[];
    for k=2:Len
        tmp2=[tmp2 log(sqrt(sum((myrou{k-1}(:)-tmp(:)).^2)))];
    end
    RMSerrorLSM2(Num,:)=tmp2;
    TimeLSM=toc;
    %%QUEST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [ThetQUEST,GamQUEST,FaiQUEST,CchangeQUEST,QchangeQUEST,myrou,TimeQUEST,myrouQUEST] = AttitudeCalBatchData3('QUEST3',Gx,Gy,Gz,xg,yg,zg,xm,ym,zm,AdaptiveRegulation,fg,h,Beta,noiseG,noiseA,noiseM,frequency,1,0);
    RecTimeQUEST(Num,:)=TimeQUEST;
    [~,Ind1]=max(abs(QchangeQUEST(:,1)));
    for k=2:Len
        if QchangeQUEST(Ind1,1)*QchangeQUEST(Ind1,k)<0
            QchangeQUEST(:,k)=-QchangeQUEST(:,k);
        end
    end
    TRecQ1QUEST(Num,:)=QchangeQUEST(1,:);TRecQ2QUEST(Num,:)=QchangeQUEST(2,:);TRecQ3QUEST(Num,:)=QchangeQUEST(3,:);TRecQ0QUEST(Num,:)=QchangeQUEST(4,:);
    tmp=cov(QchangeQUEST(:,2:Len)');
    traceQUEST2(Num,:)=log((myrouQUEST(1,:)-trace(tmp)).^2);
    detQUEST2(Num,:)=log((myrouQUEST(2,:)).^2);
    tmp2=[];
    for k=2:Len
        tmp2=[tmp2 log(sqrt(sum((myrou{k-1}(:)-tmp(:)).^2)))];
    end
    RMSerrorQUEST2(Num,:)=tmp2;
    TimeQUEST=toc;    
%%%%WuQuaternion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [ThetaNQDM,GamaNQDM,FaiNQDM,CchangeNQDM,QchangeNQDM,myrou,TimeNQDM,myrouNQDM] = AttitudeCalBatchData3('NQDM3',Gx,Gy,Gz,xg,yg,zg,xm,ym,zm,AdaptiveRegulation,fg,h,Beta,noiseG,noiseA,noiseM,frequency,1,0);
    RecTimeNQDM(Num,:)=TimeNQDM;
    [~,Ind1]=max(abs(QchangeNQDM(:,1)));
    for k=2:Len
        if QchangeNQDM(Ind1,k)*QchangeNQDM(Ind1,1)<0
            QchangeNQDM(:,k)=-QchangeNQDM(:,k);
        end
    end
    TRecQ1NQDM(Num,:)=QchangeNQDM(1,:);TRecQ2NQDM(Num,:)=QchangeNQDM(2,:);TRecQ3NQDM(Num,:)=QchangeNQDM(3,:);TRecQ0NQDM(Num,:)=QchangeNQDM(4,:);
    tmp=cov(QchangeNQDM(:,2:Len)');
    traceNQDM2(Num,:)=log((myrouNQDM(1,:)-trace(tmp)).^2);%/trace(tmp);
    detNQDM2(Num,:)=log((myrouNQDM(2,:)).^2);%/det(tmp));
    tmp2=[];
    for k=2:Len
        tmp2=[tmp2 log(sqrt(sum((myrou{k-1}(:)-tmp(:)).^2)))];
    end
    RMSerrorNQDM2(Num,:)=tmp2;
    TimeNQDM=toc;
%%%AlgebraicQuaternion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [ThetaAQUA,GamaAQUA,FaiAQUA,CchangeAQUA,QchangeAQUA,myrou,TimeAQUA,myrouAQUA] = AttitudeCalBatchData3('AQUA3',Gx,Gy,Gz,xg,yg,zg,xm,ym,zm,AdaptiveRegulation,fg,h,Beta,noiseG,noiseA,noiseM,frequency,1,0);
    RecTimeAQUA(Num,:)=TimeAQUA;
    [mxQ1,Ind1]=max(abs(QchangeAQUA(:,1)));
    for k=2:Len
        if QchangeAQUA(Ind1,k)*QchangeAQUA(Ind1,1)<0
            QchangeAQUA(:,k)=-QchangeAQUA(:,k);
        end
    end
    TRecQ1AQUA(Num,:)=QchangeAQUA(1,:);TRecQ2AQUA(Num,:)=QchangeAQUA(2,:);TRecQ3AQUA(Num,:)=QchangeAQUA(3,:);TRecQ0AQUA(Num,:)=QchangeAQUA(4,:);
    tmp=cov(QchangeAQUA(:,2:Len)');
    traceAQUA2(Num,:)=log((myrouAQUA(1,:)-trace(tmp)).^2);
    detAQUA2(Num,:)=log((myrouAQUA(2,:)).^2);
    tmp2=[];
    for k=2:Len
        tmp2=[tmp2 log(sqrt(sum((myrou{k-1}(:)-tmp(:)).^2)))];
    end
    RMSerrorAQUA2(Num,:)=tmp2;
    TimeAQUA=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

tmp=[RMSerrorLSM2(:),RMSerrorQUEST2(:),RMSerrorAQUA2(:),RMSerrorNQDM2(:)];subplot(2,2,1);boxplot(tmp,'labels',{'LSRM','QUEST','AQUA','NQDM'});grid on;ylabel('square error of elements')
tmp=[traceLSM2(:),traceQUEST2(:),traceAQUA2(:),traceNQDM2(:)];subplot(2,2,2);boxplot(tmp,'labels',{'LSRM','QUEST','AQUA','NQDM'});grid on;ylabel('square error of trace')
tmp=[detLSM2(:),detQUEST2(:),detAQUA2(:),detNQDM2(:)];subplot(2,2,3);boxplot(tmp,'labels',{'LSRM','QUEST','AQUA','NQDM'});grid on;ylabel('square of determinant')
subplot(2,2,4);bar(["LSRM","QUEST","AQUA","NQDM"],[mean(RecTimeLSM(:)),mean(RecTimeQUEST(:)),mean(RecTimeAQUA(:)),mean(RecTimeNQDM(:))]);grid on;ylabel('time (s)')


