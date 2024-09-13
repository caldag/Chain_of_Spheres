% This code calculates the radiation force Frad and radiation torque Tac
% acting on the helix.

function [Frad, Tac] = arf_spharr_calc_st(xp,yp,zp,R,hel_length,hel_mar,n_sph,a,Yst,pamp,f,comv,phi,nlam,MASSsc,LENGTHsc,TIMEsc,timek)

c=1480/(LENGTHsc*TIMEsc); % Speed of sound, water
lamb=c/f; % Acoustic wavelength
k=2*pi/lamb; % Wave number
rholiq=1000*LENGTHsc^2/(MASSsc)*TIMEsc; %Density, water

A=pamp/(rholiq*2*pi*f); % Potential amplitude
Ea=1/2*rholiq*k^2*A^2; % Acoustic energy density

hel_lamb=hel_length/nlam; %Helix geometry, coordinates

%Sphere placement (local frame)
hel_z=-nlam*hel_lamb/2+hel_length/(n_sph*2):nlam*hel_lamb/n_sph:nlam*hel_lamb/2-nlam*hel_lamb/(n_sph*2);
hel_y=hel_mar*cos(2*pi/hel_lamb.*hel_z);
hel_x=hel_mar*sin(2*pi/hel_lamb.*hel_z);

% Rotation and translation to obtain current position and orientation
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

Fac=0; %Radiation force
Tac_i=zeros(length(hel_z),3); %Radiation torque

for ii=1:length(hel_z) % Calculate force and torque for each sphere
    Fac(ii)=pi*a^2*Ea*Yst*sin(2*k*(hel_z(ii)+phi)); 
    %Acoustic radiation force on an individual sphere, standing wave
    sphpos=rotpos(1:3,ii); %Position of that sphere
    Tac_i(ii,:)=cross(sphpos-comvrot,[0; 0; Fac(ii)]); %Torque coming from an individual sphere, standing wave
end

Frad=sum(Fac); % Sum of forces on each sphere
Tac=sum(Tac_i); % Sum of torques due to each sphere
end