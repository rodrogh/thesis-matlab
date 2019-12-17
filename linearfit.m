function [fit, acc, slope,inter,mu]= linearfit(x,y)
% acc = 0;
% inter = 0;
% %New method (simpler)
% slope = (y(end) - y(1) )/ (x(end) - x(1));
% fit = x*slope;

%Fit coefficients
[p,~,mu]= polyfit(x,y,1);
% [p]= polyfit(x,y,1);
%p(x) = p1*x + p2
%p(1) is the slope and p(2) is the intercept of the linear predictor
%(UNSCALED)
%Fit data
fit = polyval(p,x,[],mu);
% fit = polyval(p,x);
%After scalling, calculate the slope of the linear regression. The cvalues
%contained in p are not scaled
slope = mean (diff(fit)./diff(x));
inter = 0;
%Code for obtaining performance of fitting
resid = y - fit;
SSresid = sum(resid.^2);
SStotal = (length(y)-1) * var(y);
%Performance
acc = 1 - SSresid/SStotal;
end