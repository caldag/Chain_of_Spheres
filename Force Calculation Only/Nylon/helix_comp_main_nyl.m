%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code calculates the acoustic radiation force for a helix made
% of nylon in a standing field in water as described in:
%
% H. O. Caldag, S. Yesilyurt, Acoustic radiation forces on 
% magnetically actuated helical swimmers, Physics of Fluids , 32, 
% 092012, 2020. https://doi.org/10.1063/5.0020930
%
% The code generates one set of data files used to construct the left half
% of Fig. 3 in the article. Currently, it is configured to evaluate the
% acoustic radiation force values for different minor diameters (see Fig.
% 3c). 
%
% To generate the left-hand side of Fig. 3 with all available data,
% run arfcom_comp_all.m.
%
% Results from the finite-element simulations are carried out with
% COMSOL Multiphysics.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Initializations and acoustic energy density
clearvars;
comdata=load('w_nyl_helmirsw_lamb100um.txt');

hel_com_arf=[comdata(:,2);]; % Acoustic radiation force (ARF) values from simulation data
hel_mir_list=[comdata(:,1);]; % Helix minor radii tested in simulations

for jj=1:length(hel_com_arf)
	f=1e6; % Acoustic field frequency, Hz
	c=1480; % Speed of sound (water), m/s
	lamb=c/f; % Acoustic wavelength, m
    
    k=2*pi/lamb; % Acoustic wave number, 1/m
    rholiq=1000; % Water density, kg/m3
    
    pamp=1e5; % Acoustic pressure amplitude, Pa
    A=pamp/(2*rholiq*2*pi*f); % Complex amplitude of the wave
    Ea=1/2*rholiq*k^2*A^2; % Acoustic energy density
    
    %% Geometric parameters of the helix
	
	hel_length=100e-6; % Helix length
	hel_mar=30e-6; % Helix major diameter
	hel_mir=hel_mir_list(jj); % Helix minor diameter
	hel_lamb=hel_length; % Helix wavelength is equal to helix length as it is one complete
	% rotation

	n_sph=3; % Number of spheres to use in chain-of-spheres (CoS) (starting point)
	a=hel_mir; % Initial radius of the spheres in CoS
	mdist=10*a; % Distance between spheres (initially set to a large value to
				% ensure the code enters the while loop below.
	sphloopc=0; % Loop counter as we increase the number of spheres in the representation

	% The loop below aims to maximise the number of spheres used in the approach by minimizing
	% the distance between the spheres (mdist) while placing them along the centerline 
	% of the helix. mdist=2*a is the absolute minimum and corresponds to spheres touching each other.

while mdist>2*a % While the distance is larger than the diameter of the spheres
				
sphloopc=sphloopc+1; % Increase loop counter

% Center coordinates of the helix in Comsol simulation model
hel_z=-hel_length+hel_length/(n_sph*2):hel_length/n_sph:-hel_length/(n_sph*2);
hel_y=hel_mar*cos(2*pi/(hel_lamb).*hel_z+pi/2);
hel_x=hel_mar*sin(2*pi/(hel_lamb).*hel_z+pi/2);
 
for ii=1:length(hel_z)-1 % Computing the distance between each sphere placed alongside the helix centreline
   dist(ii)=sqrt((hel_z(ii+1)-hel_z(ii))^2+(hel_y(ii+1)-hel_y(ii))^2+(hel_x(ii+1)-hel_x(ii))^2); 
end

mdist=mean(dist); % Average distance between sphere centroids
radiuslist(sphloopc)=mdist/2; % Potential sphere radius will be half of the
                              % average distance between the spheres
n_sph=n_sph+1; % Increase number of spheres
clear hel_z hel_y hel_x % Clear helix coordinates for the next loop step
end

selectedind=max(find(radiuslist>a)); % Find the maximum number of spheres that
                                   % can be placed while the distance is
                                   % larger than the radius
n_sph=2+selectedind; % The total number of spheres is 2 added to the number
                     % above (including the first and last spheres)
clear hel_z hel_y hel_x % Clear helix centerline parameters

% These are the center coordinates of the final set of spheres used in the force evaluation
hel_z=-hel_length+hel_length/(n_sph*2):hel_length/n_sph:-hel_length/(n_sph*2);
hel_y=hel_mar*cos(2*pi/(hel_lamb).*hel_z+pi/2);
hel_x=hel_mar*sin(2*pi/(hel_lamb).*hel_z+pi/2);

spharrvol(jj)=length(hel_z)*4/3*pi*a^3; % Volume of the chain of spheres 
% Note that the spheres may geometrically intersect but each is treated                                        
% individually so there is no actual intersection we need to address.

hel_vol(jj)=pi*4*hel_mir^2*hel_length/(cos(atan(2*pi*hel_mar/hel_length)))/4;
% Volume of the helix

volrat=hel_vol(jj)/spharrvol(jj); % Ratio of volumes
a=a.*volrat^(1/3); % Adjusting the sphere radius to match the total chain of 
% spheres volume to the helix volume

%%%%% Uncomment below to plot the coordiantes of the centres of the chain
%%%%% of spheres.
% figure;clf;
% plot3(hel_z,hel_y,hel_x,'-o');
% xlabel('x');ylabel('y');zlabel('z');
% grid on; title('Sphere Centroids');
%%%%%
%%%%%

%% Force function and force evaluation

Yst=hasegawa79_rad_force_fct(k*a); % Radiation force function for a sphere in a
% standing acoustic field

Fac=0; % Acoustic radiation force components

phi=0; % Phase angle for the standing field

for ii=1:length(hel_z) % Evaluating the forces on each sphere
    Fac(ii)=pi*a^2*Ea*Yst*sin(2*k*(hel_z(ii)+phi));
end

% Total radiation force is the summation of forces on each sphere
hel_has_arf(jj)=sum(Fac);

disp(strcat(['Minor radius of helix: ' num2str(hel_mir_list(jj)*1e6) ' um']));
disp(strcat(['Force ratio (Analytic/COMSOL): ' num2str(hel_has_arf(jj)/hel_com_arf(jj))]));
clearvars -except hel_has_arf n_sph_list hel_mir_list hel_com_arf jj spharrvol hel_vol
end


errv=(hel_has_arf-hel_com_arf')./hel_com_arf'; % Relative error percentage

% Plotting the results
% First plot shows the force values, second plot shows the error
figure;
subplot(211)
plot(2*hel_mir_list*10^6,hel_com_arf,'o-','LineWidth',2);
hold on
plot(2*hel_mir_list*10^6,hel_has_arf,'o-','LineWidth',2);grid on; axis tight;
ylabel('Radiation Force [N]')
legend('Comsol Simulation','Chain of Spheres');
set(gca,'FontSize',13)
subplot(212)
plot(2*hel_mir_list*10^6,errv*100,'o-','LineWidth',2);grid on; axis tight;
xlabel('$d [\mu m]$');ylabel('Relative Error [\%]');set(gca,'FontSize',13)
set(gca,'FontSize',13)
set(gcf,'Position',[463 153 572 777])
