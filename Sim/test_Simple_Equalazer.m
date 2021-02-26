function y_Equalized = test_Simple_Equalazer(tx_s, rx_s, eq_array, tmp, SNR)

y_Equalized  =  0 * rx_s;

% eq = comm.LinearEqualizer('Algorithm', 'RLS', 'NumTaps', 2, 'ForgettingFactor', 0.9, 'ReferenceTap', 1);


for k=1:length(rx_s(:,1))
    
    %% MMSE weights
%     chTaps = conj(rx_s(1,1:5))./conj(tx_s(1,1:5));
%     wgts = mmseweights(eq_array{1, k}, chTaps, SNR);
%     eq_array{1, k}.Channel = chTaps';
    
    if tmp
%         y = eq_array{1, k}(rx_s(k,:)'); % For MLSE and CMA
        y = eq_array{1, k}(rx_s(k,:)', tx_s(k,:)'); % For another EQ
    else
        y = eq_array{1, k}(rx_s(k,:)'); % For MLSE and CMA
    end
    
%     [y, ~, ~] = eq(rx_s(k,:)', tx_s(k,:)'); % For one EQ for all channels

    y_Equalized(k,:)=y'; 
%     release(eq_array{1, k})
end

end   
