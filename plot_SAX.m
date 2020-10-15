% Plot time series with its SAX representation
% 
% Inputs:   PAA:            PAA representation of time series
%           str:            SAX     "   "   "
%           data:           time series
%           nseg:           number of segments (length of PAA,SAX)
%           win_size:       window length
%           cutlines:       SAX quantization intervals
%           f:              density function 
%           x:              sample space of f
%           xdata:          timestamps of data samples
%           codewords:      Numeric values for each of the symbols in the
%                           codebook of SAX
%           alphabet_size:  size of SAX's codebook
%                   
% Author: Konstantinos Bountrogiannis
% Contact: kbountrogiannis@gmail.com
% Date: July 2020
function plot_SAX(PAA,str,data,nseg,win_size,cutlines,f,x,xdata,codewords,alphabet_size)
    
    if isempty(xdata)
        xdata = 1:length(data);
    end
    [d,data_len] = size(data);
    data = data';
    [~,PAA_len] = size(PAA);
    % plot the PAA segments
    for i = 1:d
        PAA_plot_i = repmat(PAA(i,:)', 1, win_size);
        PAA_plot(:,i) = reshape(PAA_plot_i', 1, numel(PAA_plot_i))';
    end
    % Fix mis-matching resolution
    if length(PAA_plot) ~= data_len
        for i = 1:d
            PAA_plot_hat(:,i) = changeres(PAA_plot(:,i),data_len);
        end
        PAA_plot = PAA_plot_hat;
    end
    figure;
%     if d==1
%         plot(xdata,PAA_plot,'r');
%         xlim([1 data_len]);
%     else %d==2
%         plot3(PAA_plot(:,1)',1:data_len,PAA_plot(:,2)','r');
%     end
    hold on;

    plot(xdata,data,'.-','MarkerSize',8,'MarkerEdgeColor',[0.25 0.25 0.25]);
        
    % draw the gray guide lines in the background
    for i = 1:d
        guidelines = repmat(cutlines(i,:)', 1, data_len);    
        plot(xdata,guidelines', 'color', [0.4 0.4 0.4]);
    end

    if ~isempty(guidelines)
        color = {[1 0.5 0],[1 1 0],[0.5 1 0],[0 1 0],[0 1 0.5],...
            [0 1 1],[0 0.5 1],[0 0 1],[0.5 0 1],[1 0 1],[1 0 0.5],...
            [1 0 0],[1 0.5 0],[1 1 0]};
    else
        color = {'r'};
    end
%     symbols = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j','k','l',...
%                 'm','n','o','p','q','r','s','t','u','v','w','x','y','z'};
    symbols = cellstr(num2str(d2b(0:alphabet_size-1)));
    symbols = regexprep(symbols, '\s+', '');

    % high-light the segments and assign them to symbols
    for i = 1 : nseg

        % get the x coordinates for the segments
        x_start = (i-1) * win_size + 1;
        x_end   = x_start + win_size - 1;
        x_mid   = x_start + (x_end - x_start) / 2;

        % color-code each segment
        colorIndex = rem(str(i),length(color))+1;

        % draw the segments
        plot(xdata(x_start:x_end),PAA_plot(x_start:x_end), 'color', color{colorIndex}, 'linewidth', 2);

%         % show symbols
%         if alphabet_size>1
%             text(floor(x_start+x_mid-4)/2+1, PAA_plot(x_start)+1/floor(max(data)-min(data))/2,...
%                 symbols(str(i)), 'fontsize', 9, 'FontWeight', 'bold', 'Color', 'r');
%         end
    end

    if ~isempty(f)
        for i = 1:d
            plot_shaded(f(i,:)*data_len/4/max(abs(f(i,:))),x(i,:));
        end
    end

end

