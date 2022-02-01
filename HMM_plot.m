function [] = HMM_plot(estEmis, binsize)
%Takes in an emission matrix and bin size, and creates a plot of neural firing frequencies in Up states and Down states
    % estEmis is an output of fitHMM.m
    
    numChannel = length(estEmis);
    neuron_list = 1:numChannel;
    %nBoot = 10;
    %bootsam = unidrnd(numState,numState,nBoot);

    figure()
    plot(neuron_list, estEmis(1,:) * (1/binsize))
    hold on
    plot(neuron_list, estEmis(2,:) * (1/binsize))
    hold on
    xlim([0 numChannel+1])
    ylim([0 max(estEmis(2,:)) * (1/binsize) + 1])
    xlabel('Channel number')
    ylabel('Firing rate, Hz')


    
    %[bootstat, bootTrans, bootEmis, bootPi0] = bootstrapLL(emissionSeq, bootsam, numState, binSize, numIter);

    
    %hf = figure;
    %set(hf,'visible','off');
    %for iState = 1:numState
    %    errorbar(1:numChannel, data.estEmis(iState,:)/data.binSize, (data.estEmis(iState,:) - squeeze(data.estEmisBcb(1,iState,:))')/data.binSize, ...
    %        (squeeze(data.estEmisBcb(2,iState,:))'-data.estEmis(iState,:))/data.binSize, '-o','Color', colors(iState,:),'LineWidth',3,'MarkerFaceColor', ...
    %        colors(iState,:),'MarkerEdgeColor', colors(iState,:));
    %    hold on;
    %end;
    %box off;
    %ylim([0 1.05*max(max(data.estEmisBcb(2,:,:)))/data.binSize]);
    %xlim([0 numChannel+1])
    %set(gca,'fontsize',resultSave.axisFontSize);
    %set(gca, 'Color', 'none');
    %xlabel('Channel number','fontsize', resultSave.labelFontSize);
    %ylabel('Firing rate, Hz','fontsize', resultSave.labelFontSize);
    %tickMarks = 1:floor(size(data.estEmis,2)/4):size(data.estEmis,2);
    %set(gca,'XTick', tickMarks);
    %fileName = 'emissionMatrix.pdf';
    %export_fig(fileName, '-transparent');
    %close(hf); 


end

