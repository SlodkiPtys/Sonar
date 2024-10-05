%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Program: layers.m
%Description: Matlab program to calculate the SWE of
% the snowpack using reflection coefficients
% and frequency shifts
%Dependencies: Matlab v7.1 R14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The inputs to this program are the frequencies (frequency.txt) at which
%each of the peaks occur in the homodyned response, and the
%reflection coefficients (reflection_coefficients.txt)
f = load('frequency.txt');
reflectioncoeff = load('reflection_coefficients.txt');
N = length(reflectioncoeff);
%allocate memory for the arrays
porosity = zeros(N);
tortuosity = zeros(N);
density = zeros(N);
speed = zeros(N);
%define gamma constant
%gamma = 0.80 for the Saskatchewan sites
%gamma = 0.70 for the Lake O’Hara
gamma = 0.70;
%define the speed in the air layer above the
%surface of the snowpack
c0 = 331;
for i = 1:N-1

 R = reflectioncoeff(i);
 %---Do this only for the first interface air-snow)---
 if (i == 1)
 %define the variables for Newton iteration
 k = 0; %the interator
 x0 = 1; %the starting approximation
 max = 10; %the maximum number of iterations
 %set the starting variables
 p0 = x0;
 p1 = p0;
 while (k < max)
 p0 = p1;
 p1 = p0- (sqrt((gamma + p0 - gamma*p0)/p0) - ...
 (p0*(1 + R))/(1 - R))/(((1 - gamma)/p0 - ...
 (gamma + p0 - gamma*p0)/p0*p0)/(2*sqrt((gamma + ...
 p0 - gamma*p0)/p0)) - (1 + R)/(1 - R));
 k = k+1;
 end
 end

 %---Do this for all of the other interfaces (snow-snow)---
 if (i > 1)
 %define the variables for Newton iteration
 k = 0; %the interator

 x0 = 1; %the starting approximation
 max = 10; %the maximum number of iterations
 p0 = x0;
 p1 = p0;
 phik1 = porosity(i-1);
 alphak1 = tortuosity(i-1);
 while (k < max)
 p0 = p1;
 p1 = p0- ((sqrt((gamma + p0 - gamma*p0)/p0)*phik1)/ ...
 sqrt(alphak1) - (p0*(1 + R))/(1 - R))/...
 ((((1 - gamma)/p0 - (gamma + p0 - gamma*p0)/ ...
 p0*p0)*phik1)/ ...
 (2*sqrt(alphak1)*sqrt((gamma + p0 - gamma*p0)/p0)) ...
 - (1 + R)/(1 - R));
 k = k+1;
 end
 end

 %stuff the values in the vectors
 porosity(i) = p1;
 density(i) = 1000*(1-p1);
 tortuosity(i) = 1-gamma*(1-(1/p1));
 speed(i) = c0/sqrt(tortuosity(i));
end
%---now that you have the speed, calculate the distance---
y = zeros(N);
deltat = 1; %sweep time
B = 19980; %bandwidth (20Hz - 20 kHz)
cumulativetime = 0; %set the variable for the cumulative time
for i = 1:N

 freqshift = f(i);

 %---do this only for the first time---
 %for y0 (the distance from the source to the snowpack)
 if (i==1)
 y(i) = (freqshift*deltat*c0)/(2*B);
 cumulativetime = cumulativetime + (y(i)/c0);
 end
 %--do this all of the other times---
 %for y1,...,y_n
 if (i > 1)
 y(i) = (speed(i-1)/(2*B))*((freqshift*deltat)-2*B*cumulativetime);
 cumulativetime = cumulativetime + y(i)/speed(i-1);
 end
end
%determine the SWE
SWE = 0;
for i = 2:N % zaczynamy od 2, ponieważ ignorujemy y0
    SWE = SWE + y(i)*density(i-1);
end