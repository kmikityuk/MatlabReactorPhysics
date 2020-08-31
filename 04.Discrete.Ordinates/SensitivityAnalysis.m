clear all
close all
clc
 
D = [74 86 110 146 170];
T = [7.000213000000000e+02 6.958099000000000e+02 6.929014000000000e+02 6.904866000000000e+02 6.930377000000000e+02];
F = [3.902232000000000e+03 3.892507000000000e+03 3.884975000000000e+03 3.878095000000000e+03 3.885393000000000e+03];

plot(D, T,'b', D, F, 'r')
hold all
for i = 1:5
    plot(D(i),T(i),'b*', D(i),F(i),'r*');
end;
xlabel('Number of Directions');
ylabel('Neutron Flux [1/cm^2s]');
legend('Thermal Flux','Fast Flux');
title('Fast and Thermal Angular Fluxes in Cladding');

