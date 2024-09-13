%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates the plots comparing the acoustic radiation force
% computations from the FEM model and the chain-of-spheres method as
% described in:
%
% H. O. Caldag, S. Yesilyurt, Acoustic radiation forces on 
% magnetically actuated helical swimmers, Physics of Fluids , 32, 
% 092012, 2020. https://doi.org/10.1063/5.0020930
%
% The code generates the left half of Fig. 3 in the article (force compu-
% tations on a nylon helix)
%
% Results from the finite-element simulations are carried out with
% COMSOL Multiphysics. .mat files contain both the FEM results (forces) and
% chain-of-spheres results evaluated with helix_comp_main.m. In each .mat
% file, there is a hel_com_arf array containing the FEM results.
% hel_com_arf contains the chain-of-spheres results. errv stores the
% relative error.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
load('w_nyl_Nsw_data.mat'); % Data for the parameter sweep on number of helical rotations
hel_N_com_arf=hel_com_arf; % Comsol results
hel_N_has_arf=hel_has_arf; % Chain-of-spheres results
x_N_data=N_list; % x-axis of the plots (here it is number of helix rotations)
hel_N_err=errv; % Relative error

load('w_nyl_helmirsw_data.mat'); % Data for the parameter sweep on the minor radius
hel_mir_com_arf=hel_com_arf;
hel_mir_has_arf=hel_has_arf;
x_mir_data=hel_mir_list.*2;
hel_mir_err=errv;

load('w_nyl_helmarsw_data.mat'); % Data for the parameter sweep on the major radius
hel_mar_com_arf=hel_com_arf;
hel_mar_has_arf=hel_has_arf;
x_mar_data=hel_mar_list.*2;
hel_mar_err=(hel_has_arf-hel_com_arf')./(hel_com_arf');

load('w_nyl_hellengthsw_data.mat'); % Data for the parameter sweep on the helix wavelength
hel_length_com_arf=hel_com_arf;
hel_length_has_arf=hel_has_arf;
x_length_data=hel_length_list;
hel_length_err=errv;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
figure;
set(gcf,'Position',[0 0 700 1000]) % Figure positioning and settings

% Plotting starts

subplot(421);
plot(x_N_data,hel_N_com_arf,'ko-','LineWidth',1.5);
hold on
plot(x_N_data,hel_N_has_arf,'ro-','LineWidth',1.5);grid on; axis tight;
ybo=ylim;xbo=xlim;
plot(x_N_data,10000*hel_N_has_arf,'bo--','LineWidth',1.5);
axis([xbo ybo])
xlabel('$N$');ylabel('$F^{rad} [N]$')
legh=legend('FEM Results','Chain-of-Spheres Results','Relative Error');
set(legh,'Position',[0.3440,0.9339,0.3168,0.06562])
set(gca,'FontSize',13)

subplot(422)
errv=(hel_N_has_arf-hel_N_com_arf')./hel_N_com_arf';
plot(x_N_data,errv*100,'bo--','LineWidth',1.5);grid on; axis tight;
xlabel('$N$');ylabel('Relative Error [\%]');set(gca,'FontSize',13)

subplot(423)
plot(x_mir_data,hel_mir_com_arf,'ko-','LineWidth',1.5);
hold on
plot(x_mir_data,hel_mir_has_arf,'ro-','LineWidth',1.5);grid on; axis tight;
xlabel('$d \mathrm{[\mu m]}$');ylabel('$F^{rad} \mathrm{[N]}$');set(gca,'FontSize',13)

subplot(424);
plot(x_mir_data.*10^6,hel_mir_err*100,'bo--','LineWidth',1.5);grid on; axis tight;
xlabel('$d\mathrm{[\mu m]}$');ylabel('Relative Error [\%]');set(gca,'FontSize',13)


subplot(425)
plot(x_mar_data,hel_mar_com_arf,'ko-','LineWidth',1.5);
hold on
plot(x_mar_data,hel_mar_has_arf,'ro-','LineWidth',1.5);grid on; axis tight;
xlabel('$D \mathrm{[\mu m]}$');ylabel('$F^{rad} \mathrm{[N]}$');set(gca,'FontSize',13)

subplot(426);
plot(x_mar_data.*10^6,hel_mar_err*100,'bo--','LineWidth',1.5);grid on; axis tight;
xlabel('$D [\mu m]$');ylabel('Relative Error [\%]');set(gca,'FontSize',13)

subplot(427);
plot(x_length_data,hel_length_com_arf,'ko-','LineWidth',1.5);
hold on
plot(x_length_data,hel_length_has_arf,'ro-','LineWidth',1.5);grid on; axis tight;
xlabel('$\lambda_{h}$');ylabel('$F^{rad} [N]$');set(gca,'FontSize',13)

subplot(428)
plot(x_length_data,hel_length_err*100,'bo--','LineWidth',1.5);grid on; axis tight;
xlabel('$\lambda_{h}$');ylabel('Relative Error [\%]');set(gca,'FontSize',13)



