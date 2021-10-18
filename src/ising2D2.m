% Function to calculate energy and magnetization of the 2D Ising model over
% time. 
% T = Temperature, N = linear lattice size, J = Ising coupling.
function [E,M,specificHeat,magneticSusceptibility]=ising2D2(T,N,J)
%%  Initial configuration
grid = sign(.5-randi([0 1],N,N)); % Random initial configuration
% grid = ones(N,N);
%%  Initialization
t = 1000000;% Number of steps
Elist=zeros(t,1); Mlist=zeros(t,1); 

trials1 = randi(N,t,1); %cheaper to generate all at once
trials2 = randi(N,t,1); %cheaper to generate all at once
%% Calculations for the probability
case1=exp(-8/T);
case2=exp(-4/T);
case3=1; % e^0
%%  Metropolis algorithm
for i=1:t
    x=trials1(i); 
    y=trials2(i);
   % periodic boundrary constraints
    if x~=1; left=grid(y,x-1);else left=grid(y,N);end
    if x~=N; right=grid(y,x+1);else right=grid(y,1);end
    if y~=1; up=grid(y-1,x);else up=grid(N,x);end
    if y~=N; down=grid(y+1,x);else down=grid(1,x);end
    summed = left+right+up+down;
    dE=2* grid(y,x)*(summed);  % change in energy
    % For efficiency the possible probabilities were predetermined.
    switch dE
        case 8
            p = case1;
        case 4
            p = case2;
        case 0 
            p = case3;
    end 
    
    % Acceptance test (including the case dE<0).
    if dE < 0 
        grid(y,x) = -1*grid(y,x); 
    elseif (rand < p)
        grid(y,x) = -1*grid(y,x);
    end
    % Update energy and magnetization.
    Mlist(i) =(sum(sum(grid))/N);
    sumOfNeighbors = ...
      circshift(grid, [ 0  1]) ...
    + circshift(grid, [ 0 -1]) ...
    + circshift(grid, [ 1  0]) ...
    + circshift(grid, [-1  0]);
    Em = -J*grid.*sumOfNeighbors;
    Energy = sum(sum(Em));
    Elist(i) = Energy;
end
%% Display time series of energy and magnetization
Elist(Elist==0)=[];Mlist(Mlist==0)=[];
Mlist=abs(Mlist);
Mlist=Mlist/N^2; Elist=Elist/N^2;    %normalize.
%%  Magnetization and energy density
% Eliminate all configurations before thermalization.
Mlist(1:500000)=[]; Elist(1:500000)=[];
% Average over post thermalization configurations.
M=sum(Mlist)/numel(Mlist);
E=sum(Elist)/numel(Elist);
%% Heat capacity, magnetic susceptibility. 
specificHeat = (mean(Elist.^2)-mean(Elist)^2)/T^2;
magneticSusceptibility = (mean(Mlist.^2)-mean(Mlist)^2)/T;
end
