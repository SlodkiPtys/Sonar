%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Program: correct.m
%Description: Matlab program to correct for the reflection coefficients
%Dependencies: Matlab v.7.1.R14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 379.65; %the angular wavenumber in the medium
y0 = 0.30; %the distance to the reflector (snow surface)
reflection = 0.22; %the measured reflection coefficient
thetaG = 0.001; %this is a small angle to the normal
r = 0:600000; %take the evaluation of the integral as far as you can go
integrand = @(r)((exp(j.*k.*sqrt((r.^2)+y0)).*r.*besselj(0,k.*cos(thetaG).*r))./...
sqrt(r.^2+y0.^2));
part1 = real(trapz(integrand(r)));
part2 = (exp(j.*(k.^2-(k.*cos(thetaG).^2).*y0)))./((k.^2)-...
 (k.*cos(thetaG)).^2);
answer = real(-j*part1*(part2.^(-1))*reflection);
