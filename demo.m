% Evaluation of pSAX against SAX
%
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: October 2020
%

% clear
% close all
load('power.mat')
data = data(:)';
% data = generate_time_series(10000,1,'normal','walk',1);
% SOMEWHERE ABOVE THIS LINE YOU SHOULD HAVE SMH DEFINED THE 'data' OBJECT,
% which should be a single-row vector (higher dimensions are not supported
% at the moment).

% try
%     parpool; % Start parallel pool
% end
% try
%     pctRunOnAll javaaddpath ParforProgMon/java;
% end

evaluate_per_ab = 1; % 0 = evaluate for different lengths, 1 = for different alphabet sizes
plot_option = 0; % (boolean) plot SAX sequence.
normalize = 1; % (boolean) z-normalize data sequences. default = 1.
iter = 1000; % #Monte-Carlo iterations

dataset_len = length(data);
sequence_len_vector = [1920 480]; % Extract time-series subsequences with this length
nseg_vector = [32 48 64 80]; % Number of segments (dimensions) to project (compress) the above lengths onto.
alphabet_size_vector = 2.^[3 4 5]; % Cardinality of the alphabet size (codebook size)
if evaluate_per_ab
    TLB_pSAX = zeros(1,length(alphabet_size_vector));
    RMSE_pSAX = zeros(1,length(alphabet_size_vector));
    for_range_1 = 1:length(sequence_len_vector);
    for_range_2 = 1:length(alphabet_size_vector);
else
    TLB_pSAX = zeros(length(sequence_len_vector),length(nseg_vector));
    RMSE_pSAX = zeros(length(sequence_len_vector),length(nseg_vector));
    for_range_1 = 1:length(sequence_len_vector);
    for_range_2 = 1:length(nseg_vector);
end
TLB_SAX = TLB_pSAX;
RMSE_SAX = RMSE_pSAX;

for j = for_range_1
for jj = for_range_2

    if evaluate_per_ab
        alphabet_size = alphabet_size_vector(jj);
        sequence_len = sequence_len_vector(j);
        nseg = 12;
    else
        alphabet_size = 32;
        sequence_len = sequence_len_vector(j);
        nseg = nseg_vector(jj);
    end
    % Training set size is defined with some empirical heuristics. We
    % basically try to keep it at a minimum level, while attaining good
    % quality
    training_set_size = min(min(15,max(30,round(0.15*length(data)/sequence_len))),floor(length(data)/sequence_len));
    training_len = training_set_size*sequence_len;
    
    noise_power = 1;
    TLB_pSAX_i = 0; TLB_SAX_i = 0;
    RMSE_pSAX_i = 0; RMSE_SAX_i = 0;
    
%     % ppm will display the parfor progress (remove if parfor is not used)
%     ppm = ParforProgMon(['(' num2str((j-1)*length(for_range_2)+jj) '/' ...
%         num2str(length(for_range_2)*length(for_range_1)) ...
%         ') PARFOR Progress: '], iter,1, 395, 60);

    for i = 1:iter
%         close all
        
%         % Train from random sequences
%         training_set = zeros(sequence_len,training_set_size);
%         idx = randi(dataset_len-sequence_len,[1,training_set_size]);
%         for s = 1:training_set_size
%             training_set(:,s) = data(idx(s)+1:idx(s)+sequence_len)';
%         end
%         training_set = training_set(:)';
        
        % or...
        % Train from sequential data
        training_set = data(1:training_len);

        idx = randi(dataset_len-sequence_len);
        test_data = data(idx+1:idx+sequence_len); % this is the reference subsequence. See comment of query_data.
        idx = randi(dataset_len-sequence_len);
        query_data = data(idx+1:idx+sequence_len); % query_data is the sequence that is compared against the test_data
        if normalize
            training_set = znormalize(training_set);
            test_data = znormalize(test_data);
            query_data = znormalize(query_data);
        end

        % win_size is the number of data points on the raw time series that
        % will be mapped to a single symbol
        win_size = floor(sequence_len/nseg);

        % Create the training set from PAA attributes
        training_PAA = tsPAA(training_set,training_set_size*nseg);
        
        query_PAA = tsPAA(query_data,nseg);

        %% KDE + LLOYD-MAX (pSAX)
        
        % Estimate data distribution with KDE
        [f,x] = ksdensity(training_PAA','npoints',training_set_size*nseg,'Kernel','epanechnikov');
        
        % Quantize using Lloyd-Max algorithm
        [codewords,cutlines] = lloydmax(f,x,alphabet_size,training_PAA,[]);

        % Map the segments to string
        str = timeseries2symbol(test_data, sequence_len, nseg, cutlines, normalize);
        str = reshape(str,[1,nseg]);   

        % Assign values to string elements (SAX does not do this by
        % default. Cannot be used to lower-bound the raw Euclidiean distance)
        coded_data_rep = repmat(codewords(str)', 1, win_size)';
        coded_pSAX = reshape(coded_data_rep,1,numel(coded_data_rep));
        % Fix mis-matching resolution
        if length(coded_pSAX) ~= sequence_len
            coded_pSAX = changeres(coded_pSAX,sequence_len);
        end
        
        % Evaluate:
        % a. Tightness of Lower Bound
        compression_ratio = sequence_len/nseg;
        tlb = min_paa_dist(query_PAA, str, alphabet_size, compression_ratio, cutlines)/norm(query_data - test_data);
        if isinf(tlb), iter = iter-1; continue; end
        if isnan(tlb), tlb = 1; end
        TLB_pSAX_i = TLB_pSAX_i+tlb;
        % b. MSE from codeword
        RMSE_pSAX_i = RMSE_pSAX_i + sqrt(mse(test_data',coded_pSAX'));

        % Plot
        if plot_option==1 || plot_option==3
            data3_PAA = tsPAA(test_data,nseg);
            % in the line below, when the first argument is query_PAA, the
            % PAA points (mean values) are plotted. When it is codewords(str), the
            % segments have values equal to the respective codewords, instead.
            plot_SAX(data3_PAA,str,test_data,nseg,win_size,cutlines,f,x,[],codewords,alphabet_size);
            yl_temp = ylim;
            yl1 = yl_temp(1); yl2 = yl_temp(2);
        end

   
        %% ASSUMING NORMAL DISTRIBUTION (SAX)
        % set the breakpoints
        cutlines = normal_cutlines(alphabet_size);
        % Calculate only centroids using Lloyd-Max criterion (bounds are enforced)
        x = -3.5:0.01:3.5;
        f = pdf('norm',x,0,1);
        codewords = lloydmax(f,x,alphabet_size,[],cutlines,1);

        % Map the segments to string
        str = timeseries2symbol(test_data, sequence_len, nseg, cutlines, normalize);
        str = reshape(str,[1,nseg]);

        % Assign values to string elements (SAX does not do this by
        % default. Cannot be used to lower-bound the raw Euclidiean distance)
        coded_data_rep = repmat(codewords(str)', 1, win_size)';
        coded_SAX = reshape(coded_data_rep,1,numel(coded_data_rep));
        % Fix mis-matching resolution
        if length(coded_SAX) ~= sequence_len
            coded_SAX = changeres(coded_SAX,sequence_len);
        end
        
        % Evaluate:
        % a. Tightness of Lower Bound
        compression_ratio = sequence_len/nseg;
        tlb = min_paa_dist(query_PAA, str, alphabet_size, compression_ratio, cutlines)/norm(query_data - test_data);
        if isinf(tlb), iter = iter-1; continue; end
        if isnan(tlb), tlb = 1; end
        TLB_SAX_i = TLB_SAX_i+tlb;
        % b. MSE from codeword
        RMSE_SAX_i = RMSE_SAX_i + sqrt(mse(test_data',coded_SAX'));

        % Plot
        if plot_option==1 || plot_option==3
            % in the line below, when the first argument is data3_PAA, the
            % PAA points are plotted. When it is codewords(str), the
            % codewords are plotted instead.
            plot_SAX(data3_PAA,str,test_data,nseg,win_size,cutlines,f,x,[],codewords,alphabet_size);
            yl_temp = ylim;
            yl1 = min(yl1,yl_temp(1));
            yl2 = max(yl2,yl_temp(2));
        end
        
        % Scale all plots to match their y-axes
        h = findobj('type','figure'); h_count = length(h);
        for h_n = 1:h_count
            figure(h_n);  ylim([yl1 yl2]);
        end
        
%         ppm.increment();
        
    end

    RMSE_pSAX(j,jj) = RMSE_pSAX_i/iter;
    RMSE_SAX(j,jj) = RMSE_SAX_i/iter;
    TLB_pSAX(j,jj) = TLB_pSAX_i/iter;
    TLB_SAX(j,jj) = TLB_SAX_i/iter;
end
end

%% Plot metrics
if evaluate_per_ab
    xaxis_range = alphabet_size_vector;
    xlab = '\alpha';
else
    xaxis_range = nseg_vector;
    xlab = 'M';
end

for i = 1:size(RMSE_pSAX,1)
    figure;
    plot(xaxis_range,RMSE_pSAX(i,:),'k-.o');
    hold on
    plot(xaxis_range,RMSE_SAX(i,:),'k-*');
    legend('pSAX','SAX','Location','northeast');
    xlabel(xlab);
    ylabel('RMSE');
end

figure;
plot(xaxis_range,TLB_pSAX(1,:),'k-.o')
hold on
plot(xaxis_range,TLB_SAX(1,:),'k-*')
xlabel(xlab);
ylabel('TLB');
legend('pSAX','SAX','Location','southeast');
for i = 2:size(RMSE_pSAX,1)
    figure;
    plot(xaxis_range,TLB_pSAX(i,:),'k-.o')
    hold on
    plot(xaxis_range,TLB_SAX(i,:),'k-*')
    xlabel(xlab);
    ylabel('TLB');
    legend('pSAX','SAX','Location','southeast');
end

if ~evaluate_per_ab
    figure;
    h = bar3([TLB_pSAX ; TLB_SAX],0.5);
    set(gca,'XTickLabel',nseg_vector*ceil(log2(alphabet_size))/8);
    set(gca,'YTickLabel',sequence_len_vector);
    title('[pSAX, SAX]')
    [nBar, nGroup] = size([TLB_pSAX ; TLB_SAX]);
    nColors  = size(get(gcf, 'colormap'), 1);
    colorInd(1,1:4) = 1; colorInd(2,1:4) = 2; colorInd(3,1:4) = 3; colorInd(4,1:4) = 4;
    colorInd(5,1:4) = 15; colorInd(6,1:4) = 16; colorInd(7,1:4) = 17; colorInd(8,1:4) = 18;
    colorInd(9,1:4) = 25; colorInd(10,1:4) = 26; colorInd(11,1:4) = 27; colorInd(12,1:4) = 28;
    for i = 1:nGroup
       c     = get(h(i), 'CData');
       color = repelem(repmat(colorInd(:, i), 1, 4), 2, 1);
       set(h(i), 'CData', color);
    end
    zlabel('TLB');
    ylabel('sequence length');
    xlabel('bytes')
end