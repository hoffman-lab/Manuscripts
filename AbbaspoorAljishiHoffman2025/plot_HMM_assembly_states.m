function plot_HMM_assembly_states(x, stateSeq)
% x: activation signal (T x 1)
% stateSeq: HMM states (T x 1), labeled 1..nStates

    x = x(:);
    stateSeq = stateSeq(:);

    nStates = max(stateSeq);

    % -------------------------------
    % Create figure with ratio 2:1
    % -------------------------------
    figure;
    subplot(3, 1, [1 2])

    % -------------------------------
    % 1) Plot activation trace
    % -------------------------------
    plot(x, 'k', 'LineWidth', 1.2);
    ylabel('Activation');
    title('Assembly Activation');
    xlim([1 numel(x)]);

    % -------------------------------
    % 2) Plot color-coded states
    % -------------------------------
    subplot(3, 1, 3)
    imagesc(stateSeq');       % row vector visualization
    
    cmap = [72 149 239; 255 255 255; 181 23 158]/255;
    
    colormap(cmap);   % any colormap you prefer
    caxis([1 nStates]);
    ylabel('State');
    xlabel('Time');
    title('HMM States (Color-coded)');
    xlim([1 numel(x)]);
    yticks([]);               % remove Y ticks for cleaner look

end
