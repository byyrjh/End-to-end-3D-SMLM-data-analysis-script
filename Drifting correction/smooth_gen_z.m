function smoothed_trace = smooth_gen_z(data, vol_per_hyper, hyperpiecepos, smooth_f)
[num_vol num_FM] = size(data);
smoothed_trace = zeros(num_vol, num_FM);
for i = 1:num_FM
    if ~isempty(hyperpiecepos)
        num_piece = length(hyperpiecepos)-1;
        for j = 1:num_piece
            ini_idx = hyperpiecepos(j)-1;
            end_idx = hyperpiecepos(j+1)-1;
            smoothed_trace(1+ini_idx*vol_per_hyper:end_idx*vol_per_hyper,i) = medfilt1(data(1+ini_idx*vol_per_hyper:end_idx*vol_per_hyper,i),10);
            % smoothed_trace(1+ini_idx*vol_per_hyper:end_idx*vol_per_hyper,i)=smooth(data(1+ini_idx*vol_per_hyper:end_idx*vol_per_hyper,i),smooth_f);
        end
    else
        smoothed_trace(:,i) = medfilt1(data(:,i),10);
        % smoothed_trace(:,i) = smooth(data(:,i),smooth_f);
    end
end
end