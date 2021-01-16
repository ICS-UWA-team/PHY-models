% Пути
clear; close all;
addpath('./Theory');
addpath('./Watermark/BCH1/mat');
addpath('./Watermark/matlab')

% Параметры для модели канала BCH1                               
% Загрузка импульсной характеристики BCH1 первый элемент
load BCH1_001.mat; 

% Параметры симуляции исходя из полосы пропускания 5 кГц для BCH1    
QAM_ModulationOrder = 4; % QPSK для нашего случая.
Number_of_carriers = 250;
NumberOfSymbolsInTime = 14; % 50
Subcarrier_Spacing = 20; %10, 5
SamplingRate = Number_of_carriers*Subcarrier_Spacing*4;
dt = 1/SamplingRate;
% Simulation_SNR_OFDM_dB = 10;
Simulation_SNR_OFDM_dB = -10:2.5:30; % SNR for OFDM in dB. The average transmit power of all methods is the same! However, the SNR might be different due to filtering (in FOFDM and UFMC) or because a different bandwidth is used (different subcarrier spacing or different number of subcarriers).          
f_central = 35000;
pilot_sep_freq = 1;
pilot_sep_time = 1;

% Number of Monte Carlo repetitions over which we take the average
Simulation_MonteCarloRepetitions = 10;                                  

% Параметры исследуемых многочастотных схем
% OFDM parameters
OFDM_NumberOfSubcarriers     =Number_of_carriers;                                      % Number of subcarriers 
OFDM_NumberOfSymbolsInTime   = NumberOfSymbolsInTime;                                      % Number OFDM symbols in time
OFDM_SubcarrierSpacing       = Subcarrier_Spacing;                                    % Subcarrier spacing (Hz)
OFDM_CyclicPrefixLength      = 1/(NumberOfSymbolsInTime*OFDM_SubcarrierSpacing);           % Length of the cyclic prefix (s)

%% Generate " +Modulation\" Objects
% OFDM Object
OFDM = Modulation.OFDM(...
    OFDM_NumberOfSubcarriers,...                                        % Number of subcarriers
    OFDM_NumberOfSymbolsInTime,...                                      % Number OFDM symbols in time                                                 
    OFDM_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz) 
    SamplingRate,...                                                    % Sampling rate (Samples/s)                                       
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz). Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    OFDM_CyclicPrefixLength, ...                                        % Length of the cyclic prefix (s)                 
    OFDM_CyclicPrefixLength ...                                           % Length of the guard time (s), that is, zeros at the beginning and at the end of the transmission
    );

% Добавление пилотов 
ChannelEstimation_OFDM = ChannelEstimation.PilotSymbolAidedChannelEstimation(...
    'Diamond',...                                                       % Pilot pattern, 'Diamond','Rectangular', 'Custom'
    [...                                                                % Matrix that represents the pilot pattern parameters
    OFDM.Nr.Subcarriers,...                                                 % Number of subcarriers
    pilot_sep_freq; ...                                                                  % Pilot spacing in the frequency domain
    OFDM.Nr.MCSymbols,...                                                   % Number of OFDM Symbols
    pilot_sep_time ...                                                                 % Pilot spacing in the time domain
    ],...                                   
    'MovingBlockAverage' ,...                                           % Interpolation method: 'MovingBlockAverage' (takes the average of a few close pilots to estimate the channel at the data position), 'FullAverage','linear','nearest','natural'
    [pilot_sep_freq OFDM.Nr.MCSymbols] ...                                           % For 'MovingBlockAverage' defines the average region in frequency and time    
    );

NrPilotSymbols_OFDM     = ChannelEstimation_OFDM.NrPilotSymbols;
NrDataSymbols_OFDM      = OFDM.Nr.Subcarriers*OFDM.Nr.MCSymbols-NrPilotSymbols_OFDM;

% Number of samples
N_OFDM  = OFDM.Nr.SamplesTotal;
N = N_OFDM;

% PAM and QAM Object
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');

% Pre-allocate Transmit Power
Ps_OFDM   = zeros(N_OFDM,1);

% Pre-allocate Power Spectral Density
PSD_OFDM  = zeros(N_OFDM,1);

BER_OFDM_Equalized = zeros(length(Simulation_SNR_OFDM_dB), Simulation_MonteCarloRepetitions);

%% Creating cell array of equalizers for each subcarrier
eq_i = comm.LinearEqualizer('Algorithm', 'RLS', 'NumTaps', 2, 'ForgettingFactor', 0.9, 'ReferenceTap', 1);
% %eq_i = comm.DecisionFeedbackEqualizer('Algorithm', 'RLS', 'NumForwardTaps', 4, 'NumFeedbackTaps', 2, 'ForgettingFactor', 0.9, 'ReferenceTap', 1);
eq_array = {eq_i};
for j=2:Number_of_carriers
        eq_array(1, j) = {eq_i};
end

%% Start Simulation
tic
for i_rep = 1:Simulation_MonteCarloRepetitions

    % Binary data
    BinaryDataStream_OFDM  = randi([0 1], OFDM.Nr.Subcarriers  * OFDM.Nr.MCSymbols * log2(QAM.ModulationOrder), 1); 

    % Map  bits to symbols
    xD_OFDM = QAM.Bit2Symbol( BinaryDataStream_OFDM );
    x_OFDM = reshape(xD_OFDM, [OFDM.Nr.Subcarriers, OFDM.Nr.MCSymbols]);

    % Transmitted signal in the time domain
    s_OFDM  =  OFDM.Modulation( x_OFDM );

    % Свертка сигнала с нестационарным каналом в основной полосе 
    % Используется модифицированный интсрумет из Вотермарк replayfilter_bb
    r_OFDM_noNoise  = replayfilter_bb(s_OFDM, SamplingRate, h, fs_tau);
    r_OFDM_noNoise = r_OFDM_noNoise(1:length(s_OFDM));

    % Calculate the transmitted power over time
    Ps_OFDM  = Ps_OFDM  + abs(s_OFDM).^2;

    % Calculat the power spectral density
    PSD_OFDM  = PSD_OFDM  + abs(fft(s_OFDM)/sqrt(N_OFDM)).^2;

    for i_SNR = 1:length(Simulation_SNR_OFDM_dB)
        % Add noise
        r_OFDM = awgn(r_OFDM_noNoise, Simulation_SNR_OFDM_dB(i_SNR), 'measured');

        % Демодуляция сигналов
        y_OFDM  =  OFDM.Demodulation(r_OFDM);

        % Эквализация
%         myEqualizers;     
        y_Equalized_OFDM  = Simple_Equalazer(x_OFDM, y_OFDM, eq_array);

        % Detect symbols (quantization and demapping to bits)
        DetectedBitStream_Equalized_OFDM  = QAM.Symbol2Bit(y_Equalized_OFDM(:));

        % Calculate the BER
        % with pilots
%         BER_OFDM_Equalized(i_SNR, i_rep)   = mean( BinaryDataStream_OFDM(ChannelEstimation_OFDM.PilotMatrix==0)...
%             ~=DetectedBitStream_Equalized_OFDM(ChannelEstimation_OFDM.PilotMatrix==0));
        % without pilots
        BER_OFDM_Equalized(i_SNR, i_rep)   = mean( BinaryDataStream_OFDM~=DetectedBitStream_Equalized_OFDM);
          
    end
        TimeNeededSoFar = toc;
        if mod(i_rep, 5)==0
            disp([int2str(i_rep/Simulation_MonteCarloRepetitions*100) '% Completed. Time Left, approx. ' int2str(TimeNeededSoFar/i_rep*(Simulation_MonteCarloRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(Simulation_MonteCarloRepetitions-i_rep)/3600) 'hour'])
        end
end

% constell = comm.ConstellationDiagram('NumInputPorts', 2);
% constell(x_OFDM((ChannelEstimation_OFDM.PilotMatrix==0)), y_Equalized_OFDM((ChannelEstimation_OFDM.PilotMatrix==0)))

% Take "average"
Ps_OFDM  = Ps_OFDM/Simulation_MonteCarloRepetitions;
PSD_OFDM  = PSD_OFDM/Simulation_MonteCarloRepetitions;
   
%% Plot Stuff
% Define colors for different modulation schemes
ColorOFDM   = [1 0 0];
% Plot BER
figure
semilogy( Simulation_SNR_OFDM_dB , mean(BER_OFDM_Equalized,2)  ,'o-','color',ColorOFDM);  hold on;
