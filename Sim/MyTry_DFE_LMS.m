% Пути
clear; close all;
addpath('./Watermark/BCH1/mat');
addpath('./Watermark/matlab')

% Параметры для модели канала BCH1                               
% Загрузка импульсной характеристики BCH1 первый элемент
load BCH1_001.mat; 

% Параметры симуляции исходя из полосы пропускания 5 кГц для BCH1    
QAM_ModulationOrder = 4; % QPSK для нашего случая.
Number_of_carriers = 250;
NumberOfSymbolsInTime = 50; % 50
Subcarrier_Spacing = 20; %10, 5
SamplingRate = Number_of_carriers*Subcarrier_Spacing*4;
dt = 1/SamplingRate;
% Simulation_SNR_dB = 10;
Simulation_SNR_dB = 0:5:50; % SNR for OFDM in dB. The average transmit power of all methods is the same! However, the SNR might be different due to filtering (in FOFDM and UFMC) or because a different bandwidth is used (different subcarrier spacing or different number of subcarriers).          
f_central = 35000;
pilot_sep_freq = 2;
pilot_sep_time = 1;

% Number of Monte Carlo repetitions over which we take the average
memory_X = {};
memory_Y = {};
memory_Y_eq = {};
Simulation_MonteCarloRepetitions = 5;                                  

% Параметры исследуемых многочастотных схем
% FBMC parameters
FBMC_NumberOfSubcarriers = Number_of_carriers;        
FBMC_NumberOfSymbolsInTime = NumberOfSymbolsInTime;                                    
FBMC_SubcarrierSpacing = Subcarrier_Spacing;                                   
FBMC_PrototypeFilter = 'Hermite-QAM'; % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM.
FBMC_OverlappingFactor = 4; % Overlapping factor, 2,3,4,...

%% Generate " +Modulation\" Objects
% FBMC Object
FBMC = Modulation.FBMC(...
    FBMC_NumberOfSubcarriers,...                                        % Number of subcarriers
    FBMC_NumberOfSymbolsInTime,...                                      % Number FBMC symbols in time
    FBMC_SubcarrierSpacing,...                                          % Subcarrier spacing (Hz)
    SamplingRate,...                                                    % Sampling rate (Samples/s)
    0,...                                                               % Intermediate frequency of the first subcarrier (Hz).  Must be a multiple of the subcarrier spacing
    false,...                                                           % Transmit real valued signal (sampling theorem must be fulfilled!)
    FBMC_PrototypeFilter,...                                            % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM. The data rate of QAM is reduced by a factor of two compared to OQAM, but robustness in doubly-selective channels is inceased
    FBMC_OverlappingFactor, ...                                         % Overlapping factor (also determines oversampling in the frequency domain)                                   
    0, ...                                                              % Initial phase shift
    true ...                                                            % Polyphase implementation
    );
FBMC_BlockOverlapTime = (FBMC.PrototypeFilter.OverlappingFactor-1/2)*FBMC.PHY.TimeSpacing;

% Добавление пилотов 
ChannelEstimation = ChannelEstimation.PilotSymbolAidedChannelEstimation(...
    'Diamond',...                                                       % Pilot pattern, 'Diamond','Rectangular', 'Custom'
    [...                                                                % Matrix that represents the pilot pattern parameters
    Number_of_carriers,...                                                 % Number of subcarriers
    pilot_sep_freq; ...                                                                  % Pilot spacing in the frequency domain
    NumberOfSymbolsInTime,...                                                   % Number of OFDM Symbols
    pilot_sep_time ...                                                                 % Pilot spacing in the time domain
    ],...                                   
    'MovingBlockAverage' ,...                                           % Interpolation method: 'MovingBlockAverage' (takes the average of a few close pilots to estimate the channel at the data position), 'FullAverage','linear','nearest','natural'
    [pilot_sep_freq NumberOfSymbolsInTime] ...                                           % For 'MovingBlockAverage' defines the average region in frequency and time    
    );

% Number of samples
N_FBMC = FBMC.Nr.SamplesTotal;
N = N_FBMC;

% PAM and QAM Object
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
if strcmp(FBMC.Method(end-3),'O')
    % FBMC-OQAM transmission, only real valued data symbols
    PAMorQAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
else
    % FBMC-QAM transmission,  complex valued data symbols
    PAMorQAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');
end

% Pre-allocate Transmit Power
Ps_FBMC   = zeros(N_FBMC, 1);

% Pre-allocate Power Spectral Density
PSD_FBMC  = zeros(N_FBMC, 1);

%% Creating cell array of equalizers for each subcarrier
plot_sign = ["DFE_LMS_1"; "DFE_LMS_2"; "DFE_LMS_3"; "DFE_LMS_4"; "DFE_LMS_5"; "DFE_LMS_6"; "DFE_LMS_7"; "DFE_LMS_8"; "DFE_LMS_9"; "DFE_LMS_10"];
num_eq = length(plot_sign);
% eq_DFE_LMS_1 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.025, 'ReferenceTap', 1);
% eq_DFE_LMS_2 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.026, 'ReferenceTap', 1);
% eq_DFE_LMS_3 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.027, 'ReferenceTap', 1);
% eq_DFE_LMS_4 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.028, 'ReferenceTap', 1);
% eq_DFE_LMS_5 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.029, 'ReferenceTap', 1);
% eq_DFE_LMS_6 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.03, 'ReferenceTap', 1);
% eq_DFE_LMS_7 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.031, 'ReferenceTap', 1);
% eq_DFE_LMS_8 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.032, 'ReferenceTap', 1);
% eq_DFE_LMS_9 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.033, 'ReferenceTap', 1);
% eq_DFE_LMS_10 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.034, 'ReferenceTap', 1);
eq_DFE_LMS_1 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 1, 'NumFeedbackTaps', 1, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_2 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 2, 'NumFeedbackTaps', 2, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_3 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 3, 'NumFeedbackTaps', 3, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_4 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 4, 'NumFeedbackTaps', 4, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_5 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 5, 'NumFeedbackTaps', 5, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_6 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 6, 'NumFeedbackTaps', 6, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_7 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 7, 'NumFeedbackTaps', 7, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_8 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 8, 'NumFeedbackTaps', 8, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_9 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 9, 'NumFeedbackTaps', 9, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_DFE_LMS_10 = comm.DecisionFeedbackEqualizer('Algorithm', 'LMS', 'NumForwardTaps', 10, 'NumFeedbackTaps', 10, 'StepSize', 0.03, 'ReferenceTap', 1);
eq_array = {eq_DFE_LMS_1; eq_DFE_LMS_2; eq_DFE_LMS_3; eq_DFE_LMS_4; eq_DFE_LMS_5; eq_DFE_LMS_6; eq_DFE_LMS_7; eq_DFE_LMS_8; eq_DFE_LMS_9; eq_DFE_LMS_10};
BER_FBMC_Equalized = zeros(num_eq, length(Simulation_SNR_dB), Simulation_MonteCarloRepetitions);
for j=2:Number_of_carriers
    for k=1:num_eq
        eq_array(k, j) = {eval("eq_"+plot_sign(k))};
    end
end

%% Start Simulation
tic
for i_rep = 1:Simulation_MonteCarloRepetitions

    % Binary data
    BinaryDataStream_FBMC  = randi([0 1], FBMC.Nr.Subcarriers  * FBMC.Nr.MCSymbols  * log2(PAMorQAM.ModulationOrder),1);

    % Map bits to symbols
    x_FBMC  = reshape(PAMorQAM.Bit2Symbol(BinaryDataStream_FBMC), FBMC.Nr.Subcarriers, FBMC.Nr.MCSymbols);
    memory_X{i_rep} = x_FBMC;

    % Transmitted signal in the time domain
    s_FBMC  =  FBMC.Modulation( x_FBMC );

    % Свертка сигнала с нестационарным каналом в основной полосе 
    % Используется модифицированный интсрумет из Вотермарк replayfilter_bb
    r_FBMC_noNoise  = replayfilter_bb(s_FBMC, SamplingRate, h, fs_tau);
    r_FBMC_noNoise = r_FBMC_noNoise(1:length(s_FBMC));

    % Calculate the transmitted power over time
    Ps_FBMC  = Ps_FBMC  + abs(s_FBMC).^2;

    % Calculat the power spectral density
    PSD_FBMC  = PSD_FBMC  + abs(fft(s_FBMC)/sqrt(N_FBMC)).^2;

    for i_SNR = 1:length(Simulation_SNR_dB)
        % Add noise
        r_FBMC = awgn(r_FBMC_noNoise, Simulation_SNR_dB(i_SNR), 'measured');

        % Демодуляция сигналов
        y_FBMC  =  FBMC.Demodulation(r_FBMC);
        memory_Y{i_rep, i_SNR} = y_FBMC;

        for k=1:num_eq
            % Эквализация
%         myEqualizers;

        y_Equalized_FBMC = Simple_Equalazer(x_FBMC, y_FBMC, eq_array(k,:), true);

        
        memory_Y_eq{i_rep, i_SNR} = y_Equalized_FBMC;

%          y_Equalized_FBMC = y_FBMC;
        % Detect symbols (quantization and demapping to bits)
        DetectedBitStream_Equalized_FBMC  = PAMorQAM.Symbol2Bit(y_Equalized_FBMC(:));

        % Calculate the BER
        % with pilots
%         BER_FBMC_Equalized(i_SNR, i_rep)   = mean( BinaryDataStream_FBMC(ChannelEstimation.PilotMatrix==0)...
%             ~=DetectedBitStream_Equalized_FBMC(ChannelEstimation.PilotMatrix==0));
        % without pilots
        BER_FBMC_Equalized(k, i_SNR, i_rep)   = mean( BinaryDataStream_FBMC~=DetectedBitStream_Equalized_FBMC);
        
        end
        
        %BER_eq{k} = BER_FBMC_Equalized;
        
    end
        TimeNeededSoFar = toc;
        if mod(i_rep, 5)==0
            disp([int2str(i_rep/Simulation_MonteCarloRepetitions*100) '% Completed. Time Left, approx. ' int2str(TimeNeededSoFar/i_rep*(Simulation_MonteCarloRepetitions-i_rep)/60) 'min, corresponding to approx. '  int2str(TimeNeededSoFar/i_rep*(Simulation_MonteCarloRepetitions-i_rep)/3600) 'hour'])
        end
end

% constell = comm.ConstellationDiagram('NumInputPorts', 2);
% constell(x_OFDM((ChannelEstimation.PilotMatrix==0)), y_Equalized_OFDM((ChannelEstimation.PilotMatrix==0)))

% Take "average"
Ps_FBMC  = Ps_FBMC/Simulation_MonteCarloRepetitions;
PSD_FBMC = PSD_FBMC/Simulation_MonteCarloRepetitions;
   
%% Plot Stuff
% Plot BER
f = figure;
hold on;
semilogy(Simulation_SNR_dB, mean(BER_FBMC_Equalized, 3)'); 
% for k=1:num_eq
    % semilogy(Simulation_SNR_dB, mean(BER_eq{k},2)); 
% end
%set(f,{'DisplayName'}, plot_sign)
%legend show
legend(plot_sign, 'Interpreter', 'none');