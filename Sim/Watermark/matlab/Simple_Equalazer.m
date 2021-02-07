function y_Equalized = Simple_Equalazer(tx_s, rx_s, eq_array)

y_Equalized  =  0 * rx_s;

% eq = comm.LinearEqualizer('Algorithm', 'RLS', 'NumTaps', 2, 'ForgettingFactor', 0.9, 'ReferenceTap', 1);

for k=1:length(rx_s(:,1))
    [y, ~, ~] = eq_array{1, k}(rx_s(k,:)', tx_s(k,:)');
%     [y, ~, ~] = eq(rx_s(k,:)', tx_s(k,:)');
    y_Equalized(k,:)=y'; 
end

end   
