
Nelems = [10,20,30,40,50];
figure(2)
plot(Nelems,abs(real(data_tip_stress(1:5)))/abs(real(mean_stress(1))),'-o')
title('Effect of {Nelem} on Root-Normalized Tip Stress')
xlabel('Number of Spar Elements')
ylabel('Root Normalized Tip Stress')