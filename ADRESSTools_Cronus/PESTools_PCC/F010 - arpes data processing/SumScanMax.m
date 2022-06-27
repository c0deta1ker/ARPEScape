function OffsE = SumScanMax(Energy, Data)

%% - 1 - Extracting the maximum positions
OffsE = [];
figure(); hold on;
for i = 1:size(Data,3)
    % -- 
    dat     = squeeze(Data(:,:,i));
    dat     = sum(dat, 2);
    dat     = dat - min(dat(:));
    dat     = dat ./ max(dat(:));    
    [~, max_indx] = max(dat(:));
	OffsE(i) = Energy(max_indx);
    DataMax = 0;
    plot(Energy, dat, '-');
end

end
