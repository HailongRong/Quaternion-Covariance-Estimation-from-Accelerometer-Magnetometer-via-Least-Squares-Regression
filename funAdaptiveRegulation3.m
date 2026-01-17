function AdaptiveRegulation=funAdaptiveRegulation3(xg,yg,zg,doorg,ave,fg)
g=sqrt(xg.^2+yg.^2+zg.^2);
g=abs(g-fg);
tmp=g;
for k=1+ave:length(g)
    if mean(tmp(k-ave+1:k))<doorg,g(k)=1;else, g(k)=10e-9;end
end
g(1:ave)=1;
AdaptiveRegulation=g;
return;