%% Calculates and plots the thermaliation times as T is increased.

thermal_times = zeros(1,11);
k=1;
% N is increased by .06 from 2 to 2.6.
for N = 2:.06:2.6
    times = zeros(1,10);
    % The average of 10 thermalization times is taken for each T.
    for i = 1:1:10
        times(i) = thermalization(N,50,1);
    end
    
    thermal_times(k) = mean(times);
    k = k+1
end
%%
plot(2:.06:2.6,thermal_times)
xlabel('T')
ylabel('Approximate Thermalization Time(Steps)')
        