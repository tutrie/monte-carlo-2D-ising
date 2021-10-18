% For the coupling J = 1 and temperature T = 3 plot the total energy of 
% generated configurations as function of time, and determine the 
% approximate time when thermalization occurs. For this function to work it
% is required that the ising2D function returns lists of the energies and
% the magnetization instead of just the averages at the end.
% T = Temperature, N = linear lattice size, J = Ising coupling.

function thermalization_time = thermalization(T,N,J)
    %% Initialize the energy and magnetization of a 2D Isling model
    
    %T = 3;
    %N = 50;
    %J = 1;
    
    [E,M,Eraw,Mraw,grid] = ising2D(T,N,J,0);
     %[E,M] = ising2D(T,N,J,0);

    % Eliminate all configurations before thermalization.
    Ee = E(500000:end);
    % Find average energy post thermalization.
    thermalList = mean(Ee)*ones(length(E),1);
    
    thermalization_time = find(E<mean(Ee)+std(Ee) & E > mean(Ee)-std(Ee), 1, 'first');  
    
     figure;
     subplot(2,1,2)
     histfit(Eraw(thermalization_time:end))
     title(sprintf('Energy Histogram after Thermalization For T = %d', T))
     subplot(2,1,1)
     histfit(Mraw(thermalization_time:end))
     title(sprintf('Magnetization Histogram after Thermalization For T = %d', T))

     grid
     
     
     figure; 
     hold on
     plot(E,'b')
     plot(thermalList,'r')
     ylabel('Running Average of Energy')
     xlabel('N_{mc}')
     title(sprintf('Energy with equilibrium. For T = %d, N = %d, J = %d', T, N, J))
     hold off
%     histfit(E(thermalization_time:end))


