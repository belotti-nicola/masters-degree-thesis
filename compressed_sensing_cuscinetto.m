clear
close all
clc

%% Acquisizione dati 
folder = [];
file = 'data.mat';
load([folder, file])
t = data.cDAQ_time_acc_enc;
original_signal = data.Acc_Y;
Ts = mean(diff(data.cDAQ_time_acc_enc)); 
fs = 1/Ts;


%% Restrizione intervallo di campionamento

figure('Name','Restrizione intervallo di campionamento')
ta=2.95;tb=3.65;           % istanti di tempo presi in considerazione
ta_ind = find(t == ta, 1); % indice per ta nel vettore dei tempi completo 
tb_ind = find(t == tb, 1); % indice per tb nel vettore dei tempi completo 

% filtraggio definendo una fc e un BW per il passabanda + calcolo inviluppo
fc = 900;
BW = 800; 
filter_order = 100;
vibration = original_signal(ta_ind:tb_ind);
vibration_zeromean = vibration - mean(vibration);
bpf = designfilt('bandpassfir', 'FilterOrder', filter_order, 'CutoffFrequency1', fc-BW/2, ...
   'CutoffFrequency2', fc+BW/2, 'SampleRate', fs);
yBFI = filter(bpf, vibration_zeromean); 
[pEnvBFI_sk, fEnvBFI_sk, yEnvBFI_sk,~] = envspectrum(yBFI, fs, ...
    'FilterOrder', 100, 'Band', [fc-BW/2 fc+BW/2]);

subplot 311
plot(fEnvBFI_sk,pEnvBFI_sk), xlim([0 1000]), ylim([0 .1])
xlabel('Frequency [Hz]'); ylabel('Abs');
title('Risultato Approccio classico')

N = length(yBFI);
n = round(N/2);

ii = randperm(N);  % Genero un vettore 1xN di indici casuali non ripetuti compresi tra 1 e N
ii = ii(1:n);      % e prendo solo i primi n coefficienti
A = dftmtx(N)'/N;  % costruisco la matrice che rende sparso il segnale (inversa di Fourier)
A = A(ii,:);       % Prendo solo n righe casuali delle N possibili
b = yEnvBFI_sk(ii);% e di conseguenza con gli stessi indici costruisco il vettore dei termini noti

%calcolo asse frequenze dft 
L = 1:round(N/2);  
dt = 1/fs; df = 1/(N*dt);
freq = [0:N-1]*df;

subplot 312
x_bpdn = spg_bpdn(A,b,11)/n;
plot(freq(L),abs(x_bpdn(L)),'Linewidth',1.5), hold on
ncomb=5;helperPlotCombs(ncomb, [ 192 ]); grid on;
xlim([0 1000]), ylim([0 .1])
xlabel('Frequency [Hz]'); ylabel('Abs');
title('BPDN')

subplot 313
x_lasso = spg_lasso(A,b,4)/4;
plot(freq(L),abs(x_lasso(L)),'Linewidth',1.5), hold on
ncomb=5;helperPlotCombs(ncomb, [ 192 ]); grid on;
xlim([0 1000]), ylim([0 .1])
xlabel('Frequency [Hz]'); ylabel('Abs');
title('Lasso')


%% Downsample di ordine 2
figure('Name','Downsample di ordine 2')
ta=2.6;tb=4;               % istanti di tempo presi in considerazione
ta_ind = find(t == ta, 1); % indice per ta nel vettore dei tempi completo 
tb_ind = find(t == tb, 1); % indice per tb nel vettore dei tempi completo 
fs_new = fs/2;             % nuova frequenza di campionamento

% filtraggio definendo una fc e un BW per il passabanda + calcolo inviluppo
fc = 650;
BW = 1100; 
filter_order = 100;
vibration = original_signal(ta_ind:tb_ind);
length(vibration)
vibration_zeromean = vibration - mean(vibration);
vibration_zeromean = resample(vibration_zeromean,fs_new,fs);
length(vibration_zeromean)
bpf = designfilt('bandpassfir', 'FilterOrder', filter_order, 'CutoffFrequency1', fc-BW/2, ...
   'CutoffFrequency2', fc+BW/2, 'SampleRate', fs_new);
yBFI = filter(bpf, vibration_zeromean); 
[pEnvBFI_skDS, fEnvBFI_skDS, yEnvBFI_sk, ~] = envspectrum(yBFI, fs_new, ...
    'FilterOrder', 100, 'Band', [fc-BW/2 fc+BW/2]);
subplot 311
plot(fEnvBFI_skDS,pEnvBFI_skDS), xlim([0 1000]), ylim([0 .08])
title('Risultato Approccio classico')


N = length(yEnvBFI_sk);
n = round(N/2);
ii = randperm(N); ii = ii(1:n);
A = dftmtx(N)'/N; A = A(ii,:);
b = yEnvBFI_sk(ii);

%calcolo frequenze dft
L = 1:round(N/2);
dt = 1/fs_new; df = 1/(N*dt);
freq = [0:N-1]*df;
subplot 312
x_bpdnDS = spg_bpdn(A,b,10)/n;
plot(freq(L),abs(x_bpdnDS(L)),'Linewidth',1.5), hold on
ncomb=5;helperPlotCombs(ncomb, [ 192 ]); grid on;
xlim([0 1000]), ylim([0 .08])
xlabel('Frequency [Hz]'); ylabel('Abs');
title('BPDN')

subplot 313
x_lassoDS = spg_lasso(A,b,5)/5;
plot(freq(L),abs(x_lassoDS(L)),'Linewidth',1.5), hold on
ncomb=5;helperPlotCombs(ncomb, [ 192 ]); grid on;
xlim([0 1000])
xlabel('Frequency [Hz]'); ylabel('Abs');
title('Lasso')

