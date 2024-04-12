clear
clc
close all

%% Parametri modello

mm = 1000; % massa carico [kg]
mc = 300;  % massa carrello [kg]
cm = 100;  % coefficiente attrito carico [N/(m/s)]
cc = 5e3;  % coefficiente attrito carrello [N/(m/s)]
L = 2.5;   % Lunghezza cavo rigido [m]
g = 9.81;  % Accelerazione di gravità [m/s^2]

% Equilibrio
alpha_eq = pi/2;  % Posizione alpha di equilibrio (tutte le altre variabili di stato sono nulle all'equilibrio)

% Tempo di simulazione
dt = 0.01;        % Tempo di campionamento [s]
tempo = 0:dt:30;  % vettore tempo di simulazione [s]


% Inizializzazione Stato
alpha_0 = pi/2;   % Condizione iniziale alpha
alpha_dot_0 = 0;  % Condizione iniziale velocità angolare
xc_0 = 0;         % Condizione iniziale posizione carrello
xc_dot_0 = 0;     % Condizione iniziale velocità carrello


%% SISTEMA LINEARIZZATO 

% Inizializzazione Stato sistema linearizzato
alpha_0_lin = 0;        % Condizione iniziale alpha
alpha_dot_0_lin = 0;    % Condizione iniziale velocità angolare
xc_0_lin = 0;           % Condizione iniziale posizione carrello
xc_dot_0_lin = 0;       % Condizione iniziale velocità carrello


% Matrici rappresentazione state-space modello linearizzato

AL = [0,         0,         0,    1.0000;
      0,         0,    1.0000,         0;
      -32.7000,  0,    -16.6667,       0;
      -17.0040,  0,   -6.6267,   -0.1000];
  
BL = [0,  0, 0.0033, 0.0013]';

CL = [1     0     0     0;
     0     1     0     0];
 
DL = [0;0];

carroponte_lin = ss(AL,BL,CL,DL);

