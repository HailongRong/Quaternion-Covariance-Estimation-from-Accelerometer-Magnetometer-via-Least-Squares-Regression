Len1=2;Len2=100;Len3=100;Qinit=[0;0;0;1];
AdaptiveRegulation=ones(1,Len3);
QUESTP11=zeros(Len1,Len2);QUESTP22=QUESTP11;QUESTP33=QUESTP11;QUESTP00=QUESTP11;
QUESTP12=QUESTP11;QUESTP13=QUESTP11;QUESTP10=QUESTP11;
QUESTP23=QUESTP11;QUESTP20=QUESTP11;
QUESTP30=QUESTP11;
Q11QUEST=QUESTP11;Q22QUEST=QUESTP11;Q33QUEST=QUESTP11;Q00QUEST=QUESTP11;
Q12QUEST=QUESTP11;Q13QUEST=QUESTP11;Q10QUEST=QUESTP11;
Q23QUEST=QUESTP11;Q20QUEST=QUESTP11;
Q30QUEST=QUESTP11;
Coef=zeros(Len1,100);
NoiseFre=zeros(Len1,4);
frequency=200;h=0.5;fg=9.8;
Gx=zeros(1,Len3);Gy=Gx;Gz=Gx;noiseG=0;
noiseA=0.008;noiseM=0.001;Beta=0;

for Num1=1:Len1
    Num1
    noiseA=unifrnd(0.001,0.08);noiseM=unifrnd(0.0001,0.004);
%     Beta=unifrnd(-80,80);
    NoiseFre(Num1,:)=[noiseA noiseM Beta frequency];
    for Num2=1:Len2
        Qinit=unifrnd(-1,1,4,1);Qinit=Qinit./norm(Qinit);
        [xg,yg,zg,xm,ym,zm,Len,fg,h,Theta,Gama,Fai,Qchange,noiseA,noiseM,A]=SimulationData3(noiseA,noiseM,fg,h,Beta,[]);
        [ThetaQUEST,GamaQUEST,FaiQUEST,CchangeQUEST,QchangeQUEST,CovarianceQUEST,TimeQUEST] = AttitudeCalBatchData3('QUEST3',Gx,Gy,Gz,xg,yg,zg,xm,ym,zm,AdaptiveRegulation,fg,h,Beta,noiseG,noiseA,noiseM,frequency,1,0);
        [mxQ1,Ind1]=max(abs(QchangeQUEST(:,1)));
        for k=2:Len3
            if QchangeQUEST(Ind1,1)*QchangeQUEST(Ind1,k)<0
                QchangeQUEST(:,k)=-QchangeQUEST(:,k);
            end
        end

        tmpQ1=QchangeQUEST(1,:);tmpQ1=tmpQ1(:);tmpQ2=QchangeQUEST(2,:);tmpQ2=tmpQ2(:);tmpQ3=QchangeQUEST(3,:);tmpQ3=tmpQ3(:);tmpQ0=QchangeQUEST(4,:);tmpQ0=tmpQ0(:);
        tmp=cov([tmpQ1,tmpQ2,tmpQ3,tmpQ0]);
        QUESTP11(Num1,Num2)=tmp(1,1);QUESTP22(Num1,Num2)=tmp(2,2);QUESTP33(Num1,Num2)=tmp(3,3);QUESTP00(Num1,Num2)=tmp(4,4);
        QUESTP12(Num1,Num2)=tmp(1,2);QUESTP13(Num1,Num2)=tmp(1,3);QUESTP10(Num1,Num2)=tmp(1,4);
        QUESTP23(Num1,Num2)=tmp(2,3);QUESTP20(Num1,Num2)=tmp(2,4);
        QUESTP30(Num1,Num2)=tmp(3,4);
        Q11QUEST(Num1,Num2)=mean(tmpQ1.*tmpQ1);Q22QUEST(Num1,Num2)=mean(tmpQ2.*tmpQ2);Q33QUEST(Num1,Num2)=mean(tmpQ3.*tmpQ3);Q00QUEST(Num1,Num2)=mean(tmpQ0.*tmpQ0);
        Q12QUEST(Num1,Num2)=mean(tmpQ1.*tmpQ2);Q13QUEST(Num1,Num2)=mean(tmpQ1.*tmpQ3);Q10QUEST(Num1,Num2)=mean(tmpQ1.*tmpQ0);
        Q23QUEST(Num1,Num2)=mean(tmpQ2.*tmpQ3);Q20QUEST(Num1,Num2)=mean(tmpQ2.*tmpQ0);
        Q30QUEST(Num1,Num2)=mean(tmpQ3.*tmpQ0);
    end
    Input=[Q11QUEST(Num1,:);Q22QUEST(Num1,:);Q33QUEST(Num1,:);Q00QUEST(Num1,:);Q12QUEST(Num1,:);Q13QUEST(Num1,:);Q10QUEST(Num1,:);Q23QUEST(Num1,:);Q20QUEST(Num1,:);Q30QUEST(Num1,:)];
    Output=[QUESTP11(Num1,:);QUESTP22(Num1,:);QUESTP33(Num1,:);QUESTP00(Num1,:);QUESTP12(Num1,:);QUESTP13(Num1,:);QUESTP10(Num1,:);QUESTP23(Num1,:);QUESTP20(Num1,:);QUESTP30(Num1,:)];
    cA=Output*Input'/(Input*Input');
    Coef(Num1,:)=cA(:);
end
