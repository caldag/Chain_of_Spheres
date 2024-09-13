% Part of the code package related to the publication:
%
% H. O. Caldag, S. Yesilyurt, Acoustic radiation forces on 
% magnetically actuated helical swimmers, Physics of Fluids , 32, 
% 092012, 2020. https://doi.org/10.1063/5.0020930
%
% This script computes the 3-D trajectory of a slender helix under magnetic
% and standing acoustic fields. Functions called within the script are
% described below:
%
% get_mobility_icra: Computes the resistance matrix. Man & Lauga(2013)
% formulation is simplified as in Caldag & Yesilyurt (2020) ICRA paper.
%
% arf_Yst_in: Determines number of spheres and the diameter of the sphere
% that will represent the helix and the corresponding radiation force func-
% tion for such a sphere. The sphere radius is first set to minor radius of
% the helix, then the value is upscaled such that total volume of the helix
% is equal to the total volume of the chain of spheres.
% This function calls the following function:
%%%%%%%%%% hasegawa79_st: This function calculates the radiation force
%%%%%%%%%% function based on Hasegawa's (1979) formulation. This function
%%%%%%%%%% requires:
%%%%%%%%%%%%%%%%%%%% sphbes1, sphbes2: Spherical Bessel function and Hankel
%%%%%%%%%%%%%%%%%%%% function definitions, respectively.
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Fcalc and alphbetacalc: These scripts calculate the
%%%%%%%%%%%%%%%%%%%% parameters F, alpha and beta variables for the
%%%%%%%%%%%%%%%%%%%% calculation of radiation force function.
%
% arf_spharr_calc_st: Calculates the radiation force and torque.
%
% The helix is assumed to be lying along the z- axis.
% Standing wave is also in z- direction.

% For each new test, the user should change the following parameters:
% paramlist (line 43), testname (line 51), loop parameter (line 56),
% testparam (line 120). One should also comment the default value placed
% inside the loop for the tested parameter (e.g. line 81 here for Ma sweep).

clearvars;
close all;

paramlist=[0.9]; % List of parameters to be tested.
% In the current case, this parameter corresponds to Ma number (Look at the for loop below)

%% Settings
SAVE=0; % Saving the parameters on/off
FIG=1; % Plotting the results on/off
FIGSAVE=0; % Saving the figures on/off

testname='st_fm20water_n3_Masw_200kpa_test'; % File/folder identifier, used while saving
if SAVE || FIGSAVE; mkdir(testname); end % Make directory for saving the results
km = 0; % Time step counter

tic; % Time counter
for Ma=paramlist % Ma is the parameter that is swept
    %% Physical parameters
    
    visc_org=1e-3; % Viscosity (water), dimensional
    MASSsc=visc_org;  % Mass scale is the viscosity
    visc=1; % Non-dimensional viscosity
    
    c0liq=1480; % Speed of sound, water, dimensional
    
    fm_org=20; % Magnetic field rotation frequency
    TIMEsc=2*pi*fm_org; % Time scale is the angular magnetic field frequency
    fm=2*pi*fm_org/TIMEsc; % Non-dimensional magnetic field frequency
    
    th0 = pi/2; % Orientation angle of the magnetization vector
    em0 = [cos(th0);sin(th0);0]; % Magnetization vector
    ez0 = [0;0;1]; % Alignment vector
    rotdir = -1; % Rotation direction of the swimmer
    
    f_a_org=1e6; % Acoustic frequency, dimensional
    f_a=f_a_org/fm_org; % Acoustic frequency, non-dimensional
    pamp_org=200e3; % Acoustic pressure, dimensional
    pamp=pamp_org/(MASSsc*TIMEsc); % Acoustic pressure, non-dimensional
    lam_ac=c0liq/f_a_org;
    phi=0*lam_ac/8; % Acoustic wave phase angle
    
    %     Ma=0.9; %The Mason number
    % (It's commented because it's the sweeping parameter)
    Bm0=1/Ma; % Magnetic field strength
    
    %% Geometric parameters
    
    nlam = 3; % Number of rotations of the tail
    scalecoef=1; % Scaling coefficient (enlarges/shrinks the helix)
    lam_org=370e-6*scalecoef; % Helix wavelength, dimensional
    B_org=250e-6*scalecoef; % Major radius, dimensional
    b_org=40e-6*scalecoef; % Minor radius, dimensional
    htheta=atan((2*pi*B_org/lam_org)); % Helix angle theta (Check Man&Lauga 2013)
    LENGTHsc=lam_org/cos(htheta); % Length scale is the wavelength alongside the centerline
    
    lam=lam_org/LENGTHsc; % Wavelength, non-dimensional
    B=B_org/LENGTHsc; % Major radius, non-dimensional
    b = b_org/LENGTHsc; % Minor radius, non-dimensional
    hel_length=lam*nlam; % Total length of the helix, non-dimensional
    kw = 2*pi/lam; % Helix wave number
    
    ct = 2*pi*visc/(log(2*lam/b)-1/2); % Tangential drag coefficient
    cn = 2*ct; % Normal drag coefficient
    a = 1/sqrt(1+kw*kw*B*B);
    
    zcom = 0;xcom = 0;ycom = 0;
    comv=[xcom;ycom;zcom];  % Center-of-mass in local coordinates
    
    D1=get_mobility_icra(lam,nlam,B,ct); % Mobility matrix evaluation
    RLmat=inv(D1(1:3,1:3)); % Resistance matrix is the inverse
    
    R0=eye(3); % Helix major axis placed along z- axis as in Man & Lauga (2013)
    
    hel_openlength=sqrt((2*pi*B)^2+(nlam*lam)^2); % Total length of the helix when extended
    hel_area=2*pi*b^2+2*pi*b*(hel_openlength); % Helix surface area    
    km = 0; %Loop counter
    
    %% Assignment of the testparam value (Change for each different run!)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    testparam=Ma;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Numerical configuration
    
    t0 = 0*TIMEsc; % Start time, no acoustics
    tF=3*TIMEsc; % End time, no acoustics
    tac=0.2*TIMEsc; % Duration of acoustic actuation
    
    dtna = 1/4000*TIMEsc; % Time step when acoustisc is OFF
    tacst=1/(10*f_a)*TIMEsc; % Time step when acoustics is turned ON (Use lower values if necessary)
    allt = [t0:dtna:tF tF+tacst:tacst:tF+tac]; % Time array
    
    %% Vector initialization for storage
    
    clear UB UL avb avl rball beta theta allR
    K = max(size(allt));
    UB = zeros(3,K); % Linear velocities in the body frame
    UL = zeros(3,K); % Linear velocities in the lab frame
    avb = zeros(3,K); % Angular velocities in the body frame
    avl = zeros(3,K); % Angular velocities in the lab frame
    rball = zeros(3,K); % Cartesian coordinates of the swimmer (lab frame)
    beta = zeros(size(allt)); % Wobbling angle
    
    u = zeros(6,1); % Linear and angular velocity vector (instantaneous)
    U0 = 0; V0 = 0; W0 = 0; vel0 = [U0;V0;W0]; vel1= vel0; vel2 = vel0;
    % Parameters for time-stepping
    R=R0; % Rotation matrix
    wxR1 = zeros(3); wxR2 = zeros(3); % wxR matrices
    rb0 = [0;0;0]; % Initial swimmer position
    
    AB = 1; % Adams-Bashforth stepping
    CN = 0; % Crank-Nicholson stepping
    FE = 0; % Forward-Euler stepping
    k =0;
    zp=zcom; % Current cartesian positions
    yp=ycom;xp=xcom;
    
    %% Radiation force function calculation (remains the same as long as the particle size is fixed)
    [Yst,n_sph,a_a] = arf_Yst_in(rb0(3),R0,f_a,pamp,hel_length,B,b,nlam,MASSsc,LENGTHsc,TIMEsc);
    %% Simulation
    for t = allt % For all time steps
        k = k+1; % Increase step counter
        
        if allt(k)<=tF % Determine time step size
            dt=dtna;
        else
            dt=tacst;
        end
        
        % Rigid body equations are solved in the body frame
        % therefore the angular velocity is calculated in the body frame
        % in this body fixed frame, therefore the angular velocity of this
        % frame is R*w
        
        omega = R*u(4:6);   % rotation of the body wrt fixed frame
        vel0  = R*u(1:3);   % velocity of the body in the lab frame
        if AB|FE
            clear wxR0
            wxR0(:,1) = cross(omega,R0(:,1));
            wxR0(:,2) = cross(omega,R0(:,2));
            wxR0(:,3) = cross(omega,R0(:,3));
        end
        if k==1 & (AB | FE)
            RR = R0  + dt*wxR0;          % Rotation matrix
            rb = rb0 + vel0*dt;          % position
        elseif AB & k==2
            RR = R0 +  dt*(3/2*wxR0-1/2*wxR1);           % Rotation matrix
            rb = rb0 + dt*(3/2*vel0-1/2*vel1);           % position
        elseif AB & k>=3
            RR = R0  + dt*((23/12)*wxR0-(4/3)*wxR1+(5/12)*wxR2); % Rotation matrix
            rb = rb0 + dt*(23/12*vel0-4/3*vel1+5/12*vel2);      % position
        elseif CN
            clear w OHM
            w=omega;
            OHM = [0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
            RN1 = inv(eye(3)-dt/2*OHM);
            RR(:,1) = RN1*(eye(3)+dt/2*OHM)*R0(:,1);
            RR(:,2) = RN1*(eye(3)+dt/2*OHM)*R0(:,2);
            RR(:,3) = RN1*(eye(3)+dt/2*OHM)*R0(:,3);
            rb = rb0 + vel0*dt;
        else
            disp('solver is not specified')
        end
        R = RR;
        
        Bm = Bm0*[cos(t);rotdir*sin(t);0]; % Magnetic field in the lab frame
        Tm = cross(em0,R'*Bm);    % R'*Bm is the magnetic field in the body frame
        Tmarr(:,k)=Tm;
        [Frad(k), Tac]=arf_spharr_calc_st(xp,yp,zp,R,hel_length,B,n_sph,a_a,Yst,pamp,f_a,comv,phi,nlam,MASSsc,LENGTHsc,TIMEsc,k);
        % Computing the acoustic radiation force

        if allt(k)<=tF % No acoustics until t reaches tF
            Tac(:)=0;
            Frad(k)=0;
        end
        
        Fradarr(k,:)=R'*[0;0;Frad(k)];
        Tacarr(k,:)=Tac; % Storing the force and torque components
        
        u = D1*[R'*[0;0;Frad(k)];-Tm+(R'*Tac')]; % Velocity evaluation, local frame
        
        Fpr(:,k)=R*(RLmat*u(1:3)); % Propulsion force
        
        % Numeric assignments for the next step
        if AB|FE; wxR2 = wxR1; wxR1 = wxR0;end
        vel2 = vel1; vel1 = vel0; vel0 = u(1:3);
        R0 = R;
        rb0 = rb;
        
        % Storing the terms
        allR(:,k) = [R(:,1); R(:,2); R(:,3)];
        UB(:,k)    = u(1:3);
        UL(:,k)    = R*u(1:3);
        avl(:,k)   = R*u(4:6);
        avb(:,k)   = u(4:6);
        rball(:,k) = rb;
        beta(k)    = acos(dot(R(:,3),ez0));
        theta(k)   = acos(dot(em0,R'*Bm)/Bm0);
        
        % Updating the current position
        xp=rball(1,k);
        yp=rball(2,k);
        zp=rball(3,k);
    end
    
    km = km+1;
    nt = length(allt);
    scales=[MASSsc LENGTHsc TIMEsc]; % Scales used for non-dimensionalization
    %% Saving and post-processing
    
    if SAVE % Saving velocities, positions, forces, scales
        save(strcat('./', testname, '/tr_UL_', testname,'_',num2str(testparam),'.dat'),'UL','-ascii');
        save(strcat('./', testname, '/tr_avb_',testname,'_',num2str(testparam),'.dat'),'avb','-ascii');
        save(strcat('./', testname, '/tr_Fr_', testname,'_',num2str(testparam),'.dat'),'Fradarr','-ascii');
        save(strcat('./' , testname, '/tr_Fpr_', testname,'_',num2str(testparam),'.dat'),'Fpr','-ascii');
        save(strcat('./' , testname, '/tr_Tac_', testname,'_',num2str(testparam),'.dat'),'Tacarr','-ascii');
        save(strcat('./', testname,  '/tr_rb_', testname,'_',num2str(testparam),'.dat'),'rball','-ascii');
        save(strcat('./', testname,  '/tr_beta_', testname,'_',num2str(testparam),'.dat'),'beta','-ascii');
        save(strcat('./', testname,  '/tr_allR_', testname,'_',num2str(testparam),'.dat'),'allR','-ascii');
        save(strcat('./', testname,  '/tr_scales_', testname,'_',num2str(testparam),'.dat'),'scales','-ascii');
        save(strcat('./', testname,  '/tr_nsphlist_', testname,'_',num2str(testparam),'.dat'),'scales','-ascii');
    end
    
    if FIG
            figure(1) % Velocity gain plot
            lengthbeforeac=length(t0:dtna:tF)-50;
            subplot(311)
            plot(allt(lengthbeforeac:end),UL(1,lengthbeforeac:end)./UL(1,lengthbeforeac),'LineWidth',1.5);ylabel('U/U_n_o_a_c');grid on;hold on;
            subplot(312)
            plot(allt(lengthbeforeac:end),UL(2,lengthbeforeac:end)./UL(2,lengthbeforeac),'LineWidth',1.5);ylabel('V/V_n_o_a_c');grid on;hold on;
            subplot(313)
            plot(allt(lengthbeforeac:end),UL(3,lengthbeforeac:end)./UL(3,lengthbeforeac),'LineWidth',1.5);ylabel('W/W_n_o_a_c');grid on;hold on;
            xlabel('t')
        
            figure(2) % Angular velocity gain plot
            subplot(311)
            plot(allt(lengthbeforeac:end),avb(1,lengthbeforeac:end)./avb(1,lengthbeforeac),'LineWidth',1.5);ylabel('\omega_x in body');grid on;hold on;
            subplot(312)
            plot(allt(lengthbeforeac:end),avb(2,lengthbeforeac:end)./avb(2,lengthbeforeac),'LineWidth',1.5);ylabel('\omega_y in body');grid on;hold on;
            subplot(313)
            plot(allt(lengthbeforeac:end),avb(3,lengthbeforeac:end)./avb(3,lengthbeforeac),'LineWidth',1.5);ylabel('\omega_z in body');grid on;hold on;
            xlabel('t')
        
            figure(3) % Radiation force plot
            subplot(311)
            plot(allt(lengthbeforeac:end),Fradarr(lengthbeforeac:end,1),'LineWidth',1.5);ylabel('F_a_c_,_x');grid on;hold on;
            subplot(312)
            plot(allt(lengthbeforeac:end),Fradarr(lengthbeforeac:end,2),'LineWidth',1.5);ylabel('F_a_c_,_y');grid on;hold on;
            subplot(313)
            plot(allt(lengthbeforeac:end),Fradarr(lengthbeforeac:end,3),'LineWidth',1.5);ylabel('F_a_c_,_z');grid on;hold on;
            xlabel('t')
        
            figure(4); % Radiation torque plot
            subplot(311)
            plot(allt(lengthbeforeac:end),Tacarr(lengthbeforeac:end,1),'LineWidth',1.5);ylabel('T_a_c_,_x');grid on;hold on;
            subplot(312)
            plot(allt(lengthbeforeac:end),Tacarr(lengthbeforeac:end,2),'LineWidth',1.5);ylabel('T_a_c_,_y');grid on;hold on;
            subplot(313)
            plot(allt(lengthbeforeac:end),Tacarr(lengthbeforeac:end,3),'LineWidth',1.5);ylabel('T_a_c_,_z');grid on;hold on;
            xlabel('t')
        
            figure(5); % Wobbling angle plot
            plot(allt(lengthbeforeac:end),beta(lengthbeforeac:end).*180/pi,'LineWidth',1.5);ylabel('\beta');xlabel('t');grid on;hold on;
        
            figure(6); % Trajectory plots
            subplot(311);
            plot(allt(lengthbeforeac:end),rball(1,lengthbeforeac:end),'LineWidth',1.5);ylabel('x');grid on; hold on;
            subplot(312);
            plot(allt(lengthbeforeac:end),rball(2,lengthbeforeac:end),'LineWidth',1.5);ylabel('y');grid on; hold on;
            subplot(313);
            plot(allt(lengthbeforeac:end),rball(3,lengthbeforeac:end),'LineWidth',1.5);ylabel('z');grid on; hold on;xlabel('t');
        
            figure(7); % 3D Trajectory plot
            plot3(rball(3,:),rball(1,:),rball(2,:),'LineWidth',1.5);xlabel('z');ylabel('x');zlabel('y');grid on; hold on;
        
    end
    
end

% Figure fonts/settings and printing
if FIG
    figure(1);set(gcf,'Position',[513 425 706 558])
    subplot(311);set(gca,'FontSize',14);subplot(312);set(gca,'FontSize',14);subplot(313);set(gca,'FontSize',14)
    legh=legend(num2str(paramlist(:)));
    set(legh,'Location','west')
    if FIGSAVE;print(strcat('./', testname, '/tr_UB_', testname,'_'), '-dtiff','-r300');end
    
    figure(2);set(gcf,'Position',[513 425 706 558])
    subplot(311);set(gca,'FontSize',14);subplot(312);set(gca,'FontSize',14);subplot(313);set(gca,'FontSize',14)
    legh=legend(num2str(paramlist(:)));
    set(legh,'Location','west')
    if FIGSAVE;print(strcat('./', testname, '/tr_avB_', testname,'_'), '-dtiff','-r300');end
    
    figure(3);set(gcf,'Position',[513 425 706 558])
    subplot(311);set(gca,'FontSize',14);subplot(312);set(gca,'FontSize',14);subplot(313);set(gca,'FontSize',14)
    legh=legend(num2str(paramlist(:)));
    set(legh,'Location','west')
    if FIGSAVE;print(strcat('./', testname, '/tr_Frad_', testname,'_'), '-dtiff','-r300');end
    
    figure(4);set(gcf,'Position',[513 425 706 558])
    subplot(311);set(gca,'FontSize',14);subplot(312);set(gca,'FontSize',14);subplot(313);set(gca,'FontSize',14)
    legh=legend(num2str(paramlist(:)));
    set(legh,'Location','west')
    if FIGSAVE;print(strcat('./', testname, '/tr_Tac_', testname), '-dtiff','-r300');end
    
    figure(5);set(gcf,'Position',[513 425 706 558]);set(gca,'FontSize',14);
    legh=legend(num2str(paramlist(:)));
    set(legh,'Location','best')
    if FIGSAVE;print(strcat('./', testname, '/tr_beta_', testname), '-dtiff','-r300');end
    
    figure(6);set(gcf,'Position',[513 425 706 558])
    subplot(311);set(gca,'FontSize',14);subplot(312);set(gca,'FontSize',14);subplot(313);set(gca,'FontSize',14)
    legh=legend(num2str(paramlist(:)));
    set(legh,'Location','west')
    if FIGSAVE;print(strcat('./', testname, '/tr_rbsingle_', testname), '-dtiff','-r300');end
    
    
    figure(7);set(gcf,'Position',[513 425 706 558]);set(gca,'FontSize',14);
    legh=legend(num2str(paramlist(:)));
    set(legh,'Location','best')
    if FIGSAVE;print(strcat('./', testname, '/tr_rb3_', testname), '-dtiff','-r300');end
end

toc;