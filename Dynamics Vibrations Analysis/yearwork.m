close all
clear all
clc
%%
while 1
    fprintf('Select 1 for the standard structure\nSelect 2 for the modified structure with central extra beams\nSelect 3 for the modified structure with lateral extra beams');
    select = input('\n');
    if select == 1
        importfile1('standard_mkr.mat');
        break
    elseif select == 2
        importfile1('centralbeams_mkr.mat');
        break
    elseif select == 3
        importfile1('externalbeams_mkr.mat');
        break
    else 
        disp("You do not select neither 1, 2 or 3");
    end
end 

%% Length of the beams
fmax = 5; % range of the frequencies [0,5]Hz
c = 1.7; % value we choose for the parameter c
%omega maxima = c*2*pi*fmax ; where 1.5 ≤ c ≤ 3
disp ("The maximum lengths for the three types of beams are: ")

% red
mR = 60; EAR = 4e8; EJR = 1e8;

LmaxR = sqrt(pi^2*sqrt(EJR/mR)/(c*2*pi*fmax))

%green
mG = 15; EAG = 2e8; EJG = 6e7;

LmaxG = sqrt(pi^2*sqrt(EJG/mG)/(c*2*pi*fmax))

%blue
mB = 20; EAB = 2e8; EJB = 6e7;

LmaxB = sqrt(pi^2*sqrt(EJB/mB)/(c*2*pi*fmax))

%% computation of the damping parametres
% Frequencies of mode vibration under fmax
if select == 1
    in = load("mode_frequencies.mat");
elseif select == 2
    in = load("mode_frequencies_centralbeams.mat");
elseif select == 3
    in = load("mode_frequencies_externalbeams.mat");
end

% Non-dimensional damping for mode vibration
h1 = 7e-3; h2 = 5e-3; h3 = 6e-3;

% From frequencies to pulsation
omega1 = 2*pi*in.freq(1); omega2 = 2*pi*in.freq(2); omega3 = 2*pi*in.freq(3);

% Matrix operation
A = [1/(2*omega1)   omega1/2;
     1/(2*omega2)   omega2/2;
     1/(2*omega3)   omega3/2];
b = [h1; h2; h3];

disp ("The dampings parametres are: ")

x = (inv(A'*A))*A'*b; %least-square solution
alpha = x(1)
beta = x(2)

%% Impulse response
T = input("Select the maximum sampling overall time [6-24]:  ");
while 1
    if isempty(T)
        disp("You do not select a value for T");
        T = input("\nSelect the maximum sampling overall time [6-24]:  ");
    elseif T < 6 || T > 24
        disp("You do not select a valid value for T")
        T = input("\nSelect the maximum sampling overall time [6-24]:  ");
    else
        break
    end
end

disp ("The maximum value of Omegak we can consider is: ")
omega_max = 2*pi*fmax % fmax = 5Hz

% Data
dt = 0.01; % sampling time
sampFreq = 1/dt; % sampling frequency
N = T*sampFreq; % number of samples
t = 0:dt:T; % vector of time

% Creation of the "impulse" signal
T_sin = 2*0.3; % Period of the sin
omega_sin = 2*pi/T_sin; % pulsation of sin
A = 4000; % max value sin
impulse = []; % creation of a vector for the "impulse" signal
for k = 0:dt:T
    if k <= 0.3
        impulse =[impulse,A * sin(omega_sin*k)];
    else
       impulse = [impulse,0];     
    end
end
% Fourier Transform
F_impulse = fft(impulse,N);
freq = (0:length(F_impulse)-1)*sampFreq/length(F_impulse);

% Filtering of the input FT
index = find(freq > 5, 1);
for j = 1:index-1
    F_filtered(j) = F_impulse(j);
    F_filtered(N-j+1) = F_impulse(N-j+1);
end

% Plot
f_filtered = ifft(F_filtered,length(t));

figure('Name','Impulse')
subplot(2,1,1)
plot(t,impulse,t,f_filtered);
title("Comparison in time domain");legend("Original","Filtered"); xlim([0 T])
subplot(2,1,2)
plot(freq,F_impulse,freq,F_filtered);
title("Comparison in Fourier domain");legend("Original","Filtered"); xlim([0 freq(end)])

pause

%% Partitioning of the Matrices M C K
if select == 1
    Mff = M(1:105, 1:105); Mcf = M(106:111, 1:105);
    Cff = R(1:105, 1:105); Ccf = R(106:111, 1:105);
    Kff = K(1:105, 1:105); Kcf = K(106:111, 1:105);
elseif select == 2
    Mff = M(1:105, 1:105); Mcf = M(106:114, 1:105);
    Cff = R(1:105, 1:105); Ccf = R(106:114, 1:105);
    Kff = K(1:105, 1:105); Kcf = K(106:114, 1:105);
elseif select == 3
    Mff = M(1:105, 1:105); Mcf = M(106:111, 1:105);
    Cff = R(1:105, 1:105); Ccf = R(106:111, 1:105);
    Kff = K(1:105, 1:105); Kcf = K(106:111, 1:105);
end

%% Point 4 & Point 5
if select == 1 || select == 3
    n_point = 105; % n° of nodal point
    n_constr = 2; % n° of clamped point at the ground
elseif select == 2
    n_point = 105;
    n_constr = 3;
end

i = sqrt(-1);
Forces = ones(n_point,1); % Vector of input forces
R = zeros(3*n_constr,N); % Vector of constrain forces
Xf = zeros(n_point,N); % Vector of displacement of the point
check = zeros(3*n_constr,N);

for k = 0 : N/2
    omek = k*2*pi/T;
    Forces(77,1) = F_filtered(k+1);
    tf_ff = -omek^2*Mff + i*omek*Cff + Kff;
    tf_cf =-omek^2*Mcf + i*omek*Ccf + Kcf;
    
    if omek <= omega_max % deleting the components with higher values of omega
        Xf(:,k+1) = (tf_ff)^(-1)*Forces;
        R(:,k+1) = tf_cf * Xf(:,k+1);
    end
    
    if k>0 
        Xf(:,N-k+1) = conj(Xf(:,k+1));
        R(:,N-k+1) = conj(R(:,k+1));
    end
end

finale1 = ifft(Xf(76,:),N); % horizontal displacement of A
finale2 = ifft(Xf(77,:),N); % vertical displacement of A
finale3 = ifft(R(1,:),N); % horizontal constaint force on P
finale4 = ifft(R(2,:),N); % vertical constraint force on P
finale5 = ifft(R(3,:),N); % rotational constraint force on P

freq = (0:length(R)-1)*sampFreq/length(R);

figure('Name','A: Time domain')
subplot(2,1,1); plot(t(1:end-1),finale1); title("Case "+ select +": Displacement of xA"); xlim([0 T]);
subplot(2,1,2); plot(t(1:end-1),finale2); title("Case "+ select +": Displacement of yA"); xlim([0 T]);

figure('Name','P: Time domain')
subplot(3,1,1); plot(t(1:end-1),finale3); title("Case "+ select +": Horizontal constraint force"); xlim([0 T]);
subplot(3,1,2); plot(t(1:end-1),finale4); title("Case "+ select +": Vertical constraint force"); xlim([0 T]);
subplot(3,1,3); plot(t(1:end-1),finale5); title("Case "+ select +": Rotational constraint force"); xlim([0 T]);