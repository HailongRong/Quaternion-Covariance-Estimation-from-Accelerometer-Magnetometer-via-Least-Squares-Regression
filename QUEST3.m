function [Q,C,P,TimeQUEST]=QUEST3(xg,yg,zg,xm,ym,zm,fg,h,Beta,noiseA,noiseM,useRegressionModel,R)
fg1=sqrt(xg^2+yg^2+zg^2);h1=sqrt(xm^2+ym^2+zm^2);
xg=xg/fg1;yg=yg/fg1;zg=zg/fg1;xm=xm/h1;ym=ym/h1;zm=zm/h1;
b1=[xg;yg;zg];b2=[xm;ym;zm];
b1=b1./sqrt(sum(b1.*b1));b2=b2./sqrt(sum(b2.*b2));
fg=1;h=1;
r1=[0;0;-fg];r2=[h*cos(Beta);0;-h*sin(Beta)];
I=eye(3);
a2=1;a1=1;
m=a1+a2;
Sgm=1/m*(a1*b1'*r1+a2*b2'*r2);
B=1/m*(a1*b1*r1'+a2*b2*r2');
S=B+B';
z=1/m*(a1*[b1(2)*r1(3)-b1(3)*r1(2);b1(3)*r1(1)-b1(1)*r1(3);b1(1)*r1(2)-b1(2)*r1(1)]+a2*[b2(2)*r2(3)-b2(3)*r2(2);b2(3)*r2(1)-b2(1)*r2(3);b2(1)*r2(2)-b2(2)*r2(1)]);
K=[S-Sgm*I,z;z',Sgm];
[V,D]=eig(K);
Dd=diag(D);
Q=V(:,Dd==max(Dd));
if length(Q(1,:))>1, Q=Q(:,1);end
Q=Q./norm(Q);
[~,d]=max(abs(Q));
if sign(Q(d))<0
    Q=-Q;
end
q0=Q(4);q1=Q(1);q2=Q(2);q3=Q(3);
C=[q0^2+q1^2-q2^2-q3^2 2*(q1*q2-q3*q0) 2*(q1*q3+q2*q0);...
   2*(q1*q2+q3*q0) q0^2+q2^2-q1^2-q3^2 2*(q2*q3-q1*q0);...
   2*(q1*q3-q2*q0) 2*(q2*q3+q1*q0) q0^2+q3^2-q2^2-q1^2];

if useRegressionModel
    tic;
    qiqj=[q1*q1;q2*q2;q3*q3;q0*q0;q1*q2;q1*q3;q1*q0;q2*q3;q2*q0;q3*q0];
    tmp=R*qiqj;
    P11=tmp(1);P22=tmp(2);P33=tmp(3);P44=tmp(4);P12=tmp(5);P13=tmp(6);P14=tmp(7);P23=tmp(8);P24=tmp(9);P34=tmp(10);
    P=zeros(4,4);P(1,1)=P11;P(2,2)=P22;P(3,3)=P33;P(4,4)=P44;
    P(1,2)=P12;P(2,1)=P12;
    P(1,3)=P13;P(3,1)=P13;
    P(1,4)=P14;P(4,1)=P14;
    P(2,3)=P23;P(3,2)=P23;
    P(2,4)=P24;P(4,2)=P24;
    P(3,4)=P34;P(4,3)=P34;
    [V,D] = eig(P);D(D<0)=abs(D(1,1))*10e-2*rand;P=V*D*V';
    TimeQUEST=toc;
else
    tic;
    routot=noiseA*noiseM/(noiseA+noiseM);
    normcrossb1b2=(norm(cross(b1,b2)))^-2;
    tmp=(noiseM-routot)*(b1*b1')+(noiseA-routot)*(b2*b2')+routot*(b1'*b2)*((b1*b2')+(b2*b1'));
    Pqq=routot*eye(3)+normcrossb1b2*tmp;
    Pqq=1/4*[Pqq zeros(3,1)];Pqq=[Pqq;zeros(1,4)];
    qopt=[ q0 -q3  q2 q1;
           q3  q0 -q1 q2;
          -q2  q1  q0 q3
          -q1 -q2 -q3 q0];
    P=qopt*Pqq*qopt';
    [U,S,V]=svd(P);S(4,4)=S(3,3)*10e-2*rand;P=U*S*V';
    TimeQUEST=toc;
end