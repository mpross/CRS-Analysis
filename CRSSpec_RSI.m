%% Initial Parameters

lamb = 1064e-9; % Wavelength
fg=0; % Gravitational frequency
% f0 = 17.05e-3; % Resonant frequency
f0 = 17.8e-3;
g=9.81; % Gravitational acceleration
Q=294; % Quality Factor
R=15.24e-2; % Lever-arm 
kb=1.380e-23; 
T=293;    
I=0.0941;

%% Data Handling

% Data loading
data = load('CRS_Data.dat');

% Data to channels parsing
tim = data(:,1);
PD1 = (data(:, 2));
PD2 = (data(:,3));
PD3 = (data(:,4));
PD12 = (data(:, 5));
PD22 = (data(:,6));
PD32 = (data(:,7));
Seis1 = (data(:,8));
Seis2 = (data(:,9));
Seis3 = (data(:,10));

% Sampling frequency
sampF=1/(tim(2)-tim(1));


%% Calibration

[L,originalDistance,ellipseParam,signals] = ellipse_fit_single(PD1,PD2,PD3);
[L2,originalDistance,ellipseParam,signals] = ellipse_fit_single(PD12,PD22,PD32);

L2 = -L2;

% Angle and sum calculation
ang = (L-L2)/(2*R);
sm = (L+L2)/(2*R);

% CRS response inversion filter
CRSInvertFilt = zpk(-2*pi*[pairQ(f0,Q)],-2*pi*[0.01 0.01],1);
CRSInvertFilt = 1*CRSInvertFilt/abs(freqresp(CRSInvertFilt,2*pi*100));

% Applying filter
angfilt = lsim(CRSInvertFilt, ang, tim);
angfilt = angfilt(1e1*sampF:end);
timfilt = tim(1e1*sampF:end);

[b,a] = butter(2,2*1e-1/sampF, 'high');
angPlot = filter(b,a,angfilt);
timPlot = timfilt(1.5e1*sampF:end);
angPlot = angPlot(1.5e1*sampF:end);

%% ASD Calculations

Navg=11;
polyOrder=2;
[A1, ~] = asd2(PD1,1/sampF, Navg, polyOrder, @hann);
[A2, ~] = asd2(PD2,1/sampF, Navg, polyOrder, @hann);
[A3, ~] = asd2(PD3,1/sampF, Navg, polyOrder, @hann);
[A12, ~] = asd2(PD12,1/sampF, Navg, polyOrder, @hann);
[A22, ~] = asd2(PD22,1/sampF, Navg, polyOrder, @hann);
[A32, ~] = asd2(PD32,1/sampF, Navg, polyOrder, @hann);
[AT, F] = asd2(tim,1/sampF, Navg, polyOrder, @hann);
[AP, F] = asd2(L,1/sampF, Navg, polyOrder, @hann);
[AP2, F] = asd2(L2,1/sampF, Navg, polyOrder, @hann);
[AA, F] = asd2(ang,1/sampF, Navg, polyOrder, @hann);
[AS, F] = asd2(sm,1/sampF, Navg, polyOrder, @hann);
[ASe1, F] = asd2(Seis1,1/sampF, Navg, polyOrder, @hann);
[ASe2, F] = asd2(Seis2,1/sampF, Navg, polyOrder, @hann);
[ASe3, F] = asd2(Seis3,1/sampF, Navg, polyOrder, @hann);
[AAF, FA] = asd2(angfilt,1/sampF, Navg, polyOrder, @hann);

%% Noise Estimates

[S, AM2, F2, fVect2] = mccs2(L,[L2,Seis1, Seis2, Seis3], 1/sampF, Navg, polyOrder, @hann);
[S, AM22, F2, fVect2] = mccs2(L,L2, 1/sampF, Navg, polyOrder, @hann);

% Frequency domain CRS Inversion
CRSTrans=-(F2.^2-fg^2)./(F2.^2-f0^2*(1+i/Q)-fg^2);
noise=abs(AM2./CRSTrans)/(2*R);
noise2=abs(AM22./CRSTrans)/(2*R);

noiseModel=CRSTheoryNoiseModel(F2);

f = F2;
Tf=1./abs(f.^2./(f.^2-i*(f0^2/Q)-f0^2));

thermal = Tf.*sqrt(4*kb*T*f0^2./(I*2*pi*f.*Q.*((f.^2-f0^2).^2+f0^4/Q^2)));
damping = Tf.*sqrt(4*kb*T*f0./(I*2*pi*Q.*((f.^2-f0^2).^2+f0^2*f.^2/Q^2)));
dampingAng = sqrt(4*kb*T*f0./(I*2*pi*Q.*((f.^2-f0^2).^2+f0^2*f.^2/Q^2)));
thermalAng = sqrt(4*kb*T*f0^2./(I*2*pi*f.*Q.*((f.^2-f0^2).^2+f0^4/Q^2)));
%% Plots 

% ASD of HoQI Distances
figure(1)
l=loglog(F,AP, F, AP2, F2,AM2, F2, dampingAng*R, F2, thermalAng*R);
xlabel('Frequency (Hz)','Interpreter', 'latex')
ylabel('ASD (m/$\sqrt{Hz}$)','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
ylim([1e-14 1e-6])
xlim([1e-3 1e2])
legend('HoQI 1','HoQI 2','mccs2','Interpreter', 'latex')
grid on

% Time series of filtered angle
figure(2)
l=plot(timPlot, angPlot);
xlabel('Time (s)','Interpreter', 'latex')
ylabel('Angle (rad)','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
% legend('HoQI 1','HoQI 2', 'Interpreter', 'latex')
grid on

% ASD of angle readout
figure(3)
l=loglog( FA,AAF, F2,noise, F2, damping, F2, sqrt(noise.^2+damping.^2), F2, noiseModel);
xlabel('Frequency (Hz)','Interpreter', 'latex')
ylabel('ASD (rad/$\sqrt{Hz}$)','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
ylim([6e-13 1e-7])
xlim([1e-2 1e1])
legend('Inertial Angle','Readout Noise', 'Damping Noise', 'Instrument Noise','Goal','GS13 Noise', 'T240 Noise','T360 Noise',  'Interpreter', 'latex')
grid on

% Time series of raw angle
figure(4)
l=plot(tim, detrend(ang));
xlabel('Time (s)','Interpreter', 'latex')
ylabel('Angle (rad)','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
% legend('HoQI 1','HoQI 2', 'Interpreter', 'latex')
grid on

% Seismometer Subtraction
figure(5)
l=loglog( FA,AAF, F2,noise, F2, noise2, F, ASe1*1e-6, F, ASe2*1e-6, F, ASe3*1e-6);
xlabel('Frequency (Hz)','Interpreter', 'latex')
ylabel('ASD (rad/$\sqrt{Hz}$)','Interpreter', 'latex')
set(l,'LineWidth',1.5);
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
ylim([6e-13 1e-6])
xlim([5e-3 1e1])
legend('Inertial Angle','Readout Noise with sub','Readout Noise','Z','X','Y',  'Interpreter', 'latex')
grid on

% fig=figure(3)
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'CRS_Noise_RSI.pdf','-dpdf','-r1200')