function save_figure(filename)
    % SAVE_FIGURE Save the figure to a file with a specific filename.

    if ~exist('dist', 'dir')
        mkdir('dist')
    end

    fig = gcf;

    try
        starting_title = copy(get(gca, 'title'));
        starting_pos = fig.Position;

        title('');
        fig.Position = [0 0 760 600];
        exportgraphics(fig, strcat("dist/", filename));

        set(gca, 'title', starting_title);
        fig.Position = starting_pos;

    catch
        starting_pos = fig.Position;

        fig.Position = [0 0 760 600];
        exportgraphics(fig, strcat("dist/", filename));

        fig.Position = starting_pos;

    end

end
