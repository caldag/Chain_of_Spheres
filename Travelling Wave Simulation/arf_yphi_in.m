% This code calculates the sphere radius in the chain-of-spheres approach
% and then calls hasegawa77_tr_calc to evaluate the radiation force
% function for the sphere radius determined.

% Outputs:
% Y_phi: Radiation force function for the sphere with the radius a_a
% n_sph: Number of spheres in the chain-of-spheres
% a_a: Non-dimensional sphere radius in the chain of spheres approach

function [Y_phi, n_sph, a_a] = arf_yphi_in(z,R,f,pamp,hel_length,hel_mar,hel_mir,Nlam,MASSsc,LENGTHsc,TIMEsc)

rholiq=1000*LENGTHsc^2/(MASSsc)*TIMEsc; % Liquid density (water)
cliq=1480/(LENGTHsc*TIMEsc); % Liquid speed of sound (water)

lamb=cliq/f; %Wavelength
k=2*pi/lamb; %Wave number

hel_lamb=hel_length/Nlam; %Helix wavelength

n_sph=3; %Initial number of spheres
a_a=hel_mir; %Initial sphere radius
mdist=100; %Initial distanc between the spheres (large value to ensure that the while loop runs)
sphloopc=0; %Counter for the while loop below
hel_kw=2*pi/hel_lamb;


while mdist>2*a_a %While the mean distance between the spheres is larger than the minor diameter
    sphloopc=sphloopc+1; %Counter increase
    
    % Sphere placement
    hel_dz=hel_length/n_sph; %Distance between the spheres
    hel_z=-hel_length+hel_length/(n_sph*2):hel_dz:-hel_length/(n_sph*2);
    hel_x=hel_mar/hel_dz/hel_kw*(sin(hel_kw*(hel_z+hel_dz/2)) - sin(hel_kw*(hel_z-hel_dz/2)));
    hel_y=-hel_mar/hel_dz/hel_kw*(cos(hel_kw*(hel_z+hel_dz/2)) - cos(hel_kw*(hel_z-hel_dz/2)));
    
    for ii=1:length(hel_z)-1
        dist(ii)=sqrt((hel_z(ii+1)-hel_z(ii))^2+(hel_y(ii+1)-hel_y(ii))^2+(hel_x(ii+1)-hel_x(ii))^2);
    end %Calculate the distance between each sphere (They are identical but the code still calculates just in case)
    
    mdist=mean(dist); %Mean distance between the spheres
    distlist(sphloopc)=mdist/2; %Half of the distance is the radius of the spheres
    n_sph=n_sph+1; %Increase number of spheres by one
    clear hel_z hel_y hel_x
end

larger=find(distlist>a_a); % Find the sphere radius that is closest to minor radius and larger than minor radius
selectedind=length(larger);
n_sph=3+selectedind-1;
clear hel_z hel_y hel_x

hel_dz=hel_length/n_sph; % Sphere placement finalization
hel_z=-hel_length+hel_length/(n_sph*2):hel_dz:-hel_length/(n_sph*2);
hel_x=hel_mar/hel_dz/hel_kw*(sin(hel_kw*(hel_z+hel_dz/2)) - sin(hel_kw*(hel_z-hel_dz/2)));
hel_y=-hel_mar/hel_dz/hel_kw*(cos(hel_kw*(hel_z+hel_dz/2)) - cos(hel_kw*(hel_z-hel_dz/2)));

spharrvol=length(hel_z)*4/3*pi*a_a^3; % Calculate sphere array volume
hel_vol=Nlam*pi*4*hel_mir^2*hel_lamb/(cos(atan(2*pi*hel_mar/hel_lamb)))/4; % Calculate helix volume
volrat=hel_vol/spharrvol; % Find the ratio
a_a=a_a.*volrat^(1/3); % Scale up spheres such that spharrvol==hel_vol

Y_phi=hasegawa77_tr_calc(k*a_a,MASSsc,LENGTHsc,TIMEsc); % Determine radiation force function for such a sphere

end