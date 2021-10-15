function test_DTU

flist = {'low_az_15_render', 'low_2source_15_render',...
    'mid_az_15_render', 'mid_2source_15_render',...
    'high_az_15_render', 'high_2source_15_render'};

refSrcList = {-15, [-15 15], -15, [-15 15], -15, [-15 15]};

wavPath = 'DTU';

% Use durLen seconds for localisation
durLen = 1; % in seconds

for n = 1:length(flist)
    
    [sig, fsHz] = audioread(fullfile(wavPath, strcat(flist{n},'.wav')));
    
    % We use the middle durLen seconds for localisation
    mid = floor(length(sig)/2);
    len = floor(durLen * fsHz / 2);
    samples = mid-len:mid+len;
    sig = sig(samples,:);

    [probs,azimuths] = test_DNN_localisation(sig,fsHz);
    
    subplot(3, 2, n);
    hold off;
    plot(azimuths, probs, 'LineWidth', 2);
    hold on;
    refSrc = refSrcList{n};
    for m = 1:length(refSrc)
        plot(refSrc(m), 0.8, 'ro');
    end
    axis([-180 180 0 1]);
    xlabel('Azimuth', 'FontSize', 12);
    ylabel('Probability', 'FontSize', 12);
    set(gca, 'XTick', -150:30:150, 'FontSize', 12);
    set(gca, 'YTick', 0:0.2:1, 'FontSize', 12);
    title(sprintf('%s, %.1f sec', strrep(flist{n}, '_', '\_'), durLen),'FontSize', 12);
    grid on
end

fig = gcf;
fig.PaperPositionMode = 'auto';
fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
print(fig, fullfile(wavPath, 'dnnloc_DTU.pdf'), '-dpdf','-bestfit');
