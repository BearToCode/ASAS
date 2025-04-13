function save_figure(figure, filename)
    % SAVE_FIGURE Save the figure to a file with a specific filename.

    if ~exist('dist', 'dir')
        mkdir('dist')
    end

    exportgraphics(figure, strcat("dist/", filename));
end
