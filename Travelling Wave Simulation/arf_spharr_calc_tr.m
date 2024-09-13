% This function calculates the radiation force Frad and radiation torque Tac
% on the helix represented as a chain of spheres.

function [Frad, Tac] = arf_spharr_calc_tr(xp,yp,zp,R,hel_length,hel_mar,n_sph,a,Yphi,pamp,f,t,comv,phi,nlam,MASSsc,LENGTHsc,TIMEsc,timek)

c=1480/(LENGTHsc*TIMEsc); % Speed of sound, water, m/s
lamb=c/f; % Acoustic wavelength
k=2*pi/lamb; % Acoustic wave number
rholiq=1000*LENGTHsc^2/(MASSsc)*TIMEsc; % Liquid density, water, kg/m3

hel_lamb=hel_length/nlam; % Helix geometry, wavelength
hel_kw=2*pi/hel_lamb; % Helix wavenumber

hel_dz=hel_length/n_sph; %Sphere placement
hel_z=-hel_length+hel_length/(n_sph*2):hel_dz:-hel_length/(n_sph*2);
hel_x=hel_mar/hel_dz/hel_kw*(sin(hel_kw*(hel_z+hel_dz/2)) - sin(hel_kw*(hel_z-hel_dz/2)));
hel_y=-hel_mar/hel_dz/hel_kw*(cos(hel_kw*(hel_z+hel_dz/2)) - cos(hel_kw*(hel_z-hel_dz/2)));
% Helix coordinates in the local frame

% Translation and rotation for the current position and orientation
rotpos=[hel_x;hel_y;hel_z;ones(1,length(hel_z))];

Hmat=[R(1,1) R(1,2) R(1,3) xp; R(2,1) R(2,2) R(2,3) yp;R(3,1) R(3,2) R(3,3) zp;0 0 0 1];
%Hmat is the transformation matrix that carries out rotation and
%translation operations together.
rotpos=Hmat*rotpos;
hel_z=rotpos(3,:);
hel_y=rotpos(2,:);
hel_x=rotpos(1,:); % Rotated and translated helix coordinates (centers of spheres)

comvrot=Hmat*[comv;1]; % Obtaining center-of-mass in global frame
comvrot=comvrot(1:3);

Fac=0; % Radiation force
Tac_i=zeros(length(hel_z),3); % Radiation torque

for ii=1:length(hel_z) % Calculate force and torque for each sphere
    pa=-pamp; % Pressure amplitude
    A=pa/(2*rholiq*2*pi*f); % Potential amplitude
    Ea=0.5*rholiq.*k.^2*A^2; % Acoustic energy density
    Fac(ii)=pi*a^2*Ea*Yphi; % Radiation force on the current sphere
    sphpos=rotpos(1:3,ii); % Position of the current sphere
    Tac_i(ii,:)=cross(sphpos-comvrot,[0; 0; Fac(ii)]); % Torque for the current sphere
end

Frad=sum(Fac); % Sum of forces on each sphere
Tac=sum(Tac_i); % Sum of torques due to each sphere
end