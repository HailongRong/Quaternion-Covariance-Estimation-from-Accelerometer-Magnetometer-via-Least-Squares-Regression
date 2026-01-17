Start=-80;Run=1000;
Beta=(Start:5:80)/90;U=[ones(1,length(Beta));Beta;Beta.^2;Beta.^3;Beta.^4;Beta.^5;Beta.^6;Beta.^7;Beta.^8;Beta.^9;Beta.^10;Beta.^11;Beta.^12;Beta.^13;Beta.^14;Beta.^15;Beta.^16;Beta.^17;Beta.^18];
b12_1=zeros(Run,19);b12_2=b12_1;b78_1=b12_1;b78_2=b12_1;b69_1=b12_1;b69_2=b12_1;b19_2=b12_1;
a13_1=zeros(Run,length(Start:5:80));a13_2=a13_1;a14_2=a13_1;
a12_1=zeros(Run,length(Start:5:80));a12_2=a12_1;a78_1=a12_1;a78_2=a12_1;a69_1=a12_1;a69_2=a12_1;a19_2=a12_1;
lamda=10e-5;
for v=1:Run
    x=1;
    for Beta=Start:5:80
        [v Beta]
        FindRelationQuaternion_P3;
        
        Index=2; %AL(1,2)
        A=[NoiseFre(:,2)'.^2;NoiseFre(:,1)'.^2];B=Coef(:,Index)';c12=B*A'/(A*A');a12_1(v,x)=c12(1);a12_2(v,x)=c12(2);
        Index=68; %AL(7,8)
        B=Coef(:,Index)';c78=B*A'/(A*A');a78_1(v,x)=c78(1);a78_2(v,x)=c78(2);
        Index=59; %AL(6,9)
        B=Coef(:,Index)';c69=B*A'/(A*A');a69_1(v,x)=c69(1);a69_2(v,x)=c69(2);
        Index=3; %AL(1,3)
        B=Coef(:,Index)';c13=B*A'/(A*A');a13_1(v,x)=c13(1);a13_2(v,x)=c13(2);
        Index=9; %AL(1,9)
        A=NoiseFre(:,1)'.^2;B=Coef(:,Index)';c19=B*A'/(A*A');a19_2(v,x)=c19;
        Index=4; %AL(1,4)
        B=Coef(:,Index)';c14=B*A'/(A*A');a14_2(v,x)=c14;
        x=x+1;
    end
    b12_1(v,:)=a12_1(v,:)*U'/(U*U'+lamda*eye(19));b12_2(v,:)=a12_2(v,:)*U'/(U*U'+lamda*eye(19));
    b78_1(v,:)=a78_1(v,:)*U'/(U*U'+lamda*eye(19));b78_2(v,:)=a78_2(v,:)*U'/(U*U'+lamda*eye(19));
    b69_1(v,:)=a69_1(v,:)*U'/(U*U'+lamda*eye(19));b69_2(v,:)=a69_2(v,:)*U'/(U*U'+lamda*eye(19));
    b19_2(v,:)=a19_2(v,:)*U'/(U*U'+lamda*eye(19));
end

