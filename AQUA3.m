function [Q,C,P,TimeAQ]=AQUA3(xg,yg,zg,xm,ym,zm,fg,h,Beta,noiseA,noiseM)
% A Linear Kalman Filter for MARG Orientation Estimation Using the Algebraic Quaternion Algorithm
fg1=sqrt(xg^2+yg^2+zg^2);h1=sqrt(xm^2+ym^2+zm^2);
xg=xg/fg1;yg=yg/fg1;zg=zg/fg1;xm=xm/h1;ym=ym/h1;zm=zm/h1;
xg=-xg;yg=-yg;zg=-zg;
if zg>=0
    y1=sqrt((zg+1)/2);qacc=[y1;-yg/2/y1;xg/2/y1;0];
else
    y2=sqrt((1-zg)/2);qacc=[-yg/2/y2;y2;0;xg/2/y2];
end
q0=qacc(1);q1=qacc(2);q2=qacc(3);q3=qacc(4);
C=[q0^2+q1^2-q2^2-q3^2 2*(q1*q2+q3*q0) 2*(q1*q3-q2*q0);...
   2*(q1*q2-q3*q0) q0^2+q2^2-q1^2-q3^2 2*(q2*q3+q1*q0);...
   2*(q1*q3+q2*q0) 2*(q2*q3-q1*q0) q0^2+q3^2-q2^2-q1^2];
l=C*[xm;ym;zm];
lx=l(1);ly=l(2);
lxy=lx^2+ly^2;
if lx>=0
    qmag=[sqrt(lxy+lx*sqrt(lxy))/sqrt(2*lxy);0;0;ly/sqrt(2)/sqrt(lxy+lx*sqrt(lxy))];
else
    qmag=[ly/sqrt(2)/sqrt(lxy-lx*sqrt(lxy));0;0;sqrt(lxy-lx*sqrt(lxy))/sqrt(2*lxy)];
end
a0=qacc(1);a1=qacc(2);a2=qacc(3);a3=qacc(4);
m0=qmag(1);m1=qmag(2);m2=qmag(3);m3=qmag(4);
q0=a0*m0-a1*m1-a2*m2-a3*m3;
q1=a0*m1+a1*m0+a2*m3-a3*m2;
q2=a0*m2-a1*m3+a2*m0+a3*m1;
q3=a0*m3+a1*m2-a2*m1+a3*m0;

tic;
if zg>=0
    qf1=[m0 0 0 -m3 a0 -a1 -a2 0;0 m0 m3 0 a1 a0 0 a2;0 -m3 m0 0 a2 0 a0 -a1;m3 0 0 m0 0 -a2 a1 a0];
    k=sqrt(1+zg);y=lxy;b1=sqrt(y+lx*sqrt(y));b2=sqrt(y-lx*sqrt(y));
    if lx>=0
        f1f2=1/2/sqrt(2)*[0 0 1/k 0 0 0;0 -2/k yg/k^3 0 0 0;2/k 0 -xg/k^3 0 0 0;0 0 0 0 0 0;0 0 0 ly^2/b1/y lx*ly/b1/y 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 -ly*b1/y^1.5 lx*b1/y^1.5 0];
    else
        f1f2=1/2/sqrt(2)*[0 0 1/k 0 0 0;0 -2/k yg/k^3 0 0 0;2/k 0 -xg/k^3 0 0 0;0 0 0 0 0 0;0 0 0 ly*b2/y^1.5 lx*b2/y^1.5 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 -ly^2/b2/y lx*ly/b2/y 0];
    end
    f2u=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;...
        zm-(2*xg*xm+yg*ym)/k^2 -xg*ym/k^2 xg*(xg*xm+yg*ym)/k^4 1-xg^2/k^2 -xg*yg/k^2 xg;...
        -yg*xm/k^2 zm-(xg*xm+2*yg*ym)/k^2 yg*(xg*xm+yg*ym)/k^4 -xg*yg/k^2 1-yg^2/k^2 yg;...
        -xm -ym zm -xg -yg zg];
else
    qf1=[m0 0 0 -m3 a0 -a1 0 -a3;0 m0 m3 0 a1 a0 -a3 0;0 -m3 m0 0 0 a3 a0 -a1;m3 0 0 m0 a3 0 a1 a0];
    k=sqrt(1-zg);y=lxy;b1=sqrt(y+lx*sqrt(y));b2=sqrt(y-lx*sqrt(y));
    if lx>=0
        f1f2=1/2/sqrt(2)*[0 -2/k -yg/k^3 0 0 0;0 0 -1/k 0 0 0;0 0 0 0 0 0;2/k 0 xg/k^3 0 0 0;0 0 0 ly^2/b1/y -lx*ly/b1/y 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 -ly*b1/y^1.5 lx*b1/y^1.5 0];
    else
        f1f2=1/2/sqrt(2)*[0 -2/k -yg/k^3 0 0 0;0 0 -1/k 0 0 0;0 0 0 0 0 0;2/k 0 xg/k^3 0 0 0;0 0 0 ly*b2/y^1.5 -lx*b2/y^1.5 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 -ly^2/b2/y lx*ly/b2/y 0];
    end
    f2u=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;...
        zm-(2*xg*xm-yg*ym)/k^2 xg*ym/k^2 xg*(-xg*xm+yg*ym)/k^4 1-xg^2/k^2 -xg*yg/k^2 xg;...
        -yg*xm/k^2 zm-(xg*xm-2*yg*ym)/k^2 yg*(-xg*xm+yg*ym)/k^4 -xg*yg/k^2 -1+yg^2/k^2 yg;...
        xm -ym -zm xg -yg -zg];
end
J=qf1*f1f2*f2u;
P=J*[noiseA*eye(3) zeros(3);zeros(3) noiseM*eye(3)]*J';
tmp1=P(2:4,2:4);tmp2=P(2:4,1);tmp3=[P(1,2:4) P(1,1)];
tmp4=[tmp1 tmp2]; P=[tmp4;tmp3];
TimeAQ=toc;
q0=-q0;Q=[q1;q2;q3;q0];
C=[q0^2+q1^2-q2^2-q3^2 2*(q1*q2-q3*q0) 2*(q1*q3+q2*q0);2*(q1*q2+q3*q0) q0^2+q2^2-q1^2-q3^2 2*(q2*q3-q1*q0);2*(q1*q3-q2*q0) 2*(q2*q3+q1*q0) q0^2+q3^2-q2^2-q1^2];
