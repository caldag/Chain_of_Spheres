% This code calculates the sphere radius in the chain-of-spheres approach
% and then calls hasegawa79_rad_force_fct to evaluate the radiation force
% function for the sphere radius determined.

% Outputs:
% Yst: Radiation force function, standing wave.
% n_sph: Number of spheres used in the chain-of-spheres
% a_a: Sphere radius in the chain

function [Yst, n_sph, a_a] = arf_Yst_in(z,R,f,pamp,hel_length,hel_mar,hel_mir,Nlam,MASSsc,LENGTHsc,TIMEsc)

c=1480/(LENGTHsc*TIMEsc); %Speed of sound, water
lamb=c/f; %Acoustic wavelength, water
k=2*pi/lamb; % Acoustic wave number
rholiq=1000*LENGTHsc^2/(MASSsc)*TIMEsc; %Density, water

A=pamp/(2*rholiq*2*pi*f); % Amplitude of the velocity potential

Ea=1/2*rholiq*k^2*A^2; % Acoustic energy density

hel_lamb=hel_length/Nlam; % Helical wavelength

n_sph=3; % (Initial) number of spheres
a_a=hel_mir; % (Initial) sphere radius
mdist=100; % Mean distance between the spheres (large value to initiate the while loop)
sphloopc=0; % While loop counter

while mdist>2*a_a
    %Check all configurations of sphere placements where the mean distance is larger than the minor diameter
    sphloopc=sphloopc+1;

    % Center coordinates of the helix
    hel_z=-hel_length+hel_length/(n_sph*2):hel_length/n_sph:-hel_length/(n_sph*2);
    hel_y=hel_mar*cos(2*pi/(hel_lamb).*hel_z+pi/2);
    hel_x=hel_mar*sin(2*pi/(hel_lamb).*hel_z+pi/2);

    for ii=1:length(hel_z)-1
        dist(ii)=sqrt((hel_z(ii+1)-hel_z(ii))^2+(hel_y(ii+1)-hel_y(ii))^2+(hel_x(ii+1)-hel_x(ii))^2);
    end

    mdist=mean(dist); %Sphere radius for calculation from small spheres
    distlist(sphloopc)=mdist/2;
    n_sph=n_sph+1;
    clear hel_z hel_y hel_x
end

larger=find(distlist>a_a); %Find the configuration with the distance larger than the minor radius
selectedind=length(larger);
n_sph=3+selectedind-1;
clear hel_z hel_y hel_x

% Sphere placement, finalized
hel_z=-hel_length+hel_length/(n_sph*2):hel_length/n_sph:-hel_length/(n_sph*2);
hel_y=hel_mar*cos(2*pi/(hel_lamb).*hel_z+pi/2);
hel_x=hel_mar*sin(2*pi/(hel_lamb).*hel_z+pi/2);

spharrvol=length(hel_z)*4/3*pi*a_a^3; % Sphere array volume
hel_vol=Nlam*pi*4*hel_mir^2*hel_lamb/(cos(atan(2*pi*hel_mar/hel_lamb)))/4; % Actual helix volume
volrat=hel_vol/spharrvol; % Volume ratio
a_a=a_a.*volrat^(1/3); % Upscale the spheres to match the chain of spheres volume to the helix volume

%% Force function evaluation

Yst=hasegawa79_rad_force_fct(k*a_a,MASSsc,LENGTHsc,TIMEsc);

end
