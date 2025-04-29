function save_figure(filename, options)
    % SAVE_FIGURE Save the figure to a file with a specific filename.
    arguments
        filename string
        options.keep_title (1, 1) {mustBeNumericOrLogical} = false
        options.aspect_ratio_multiplier (1, 1) double = 1
    end

    if ~exist('./dist', 'dir')
        mkdir('./dist')
    end

    fig = gcf;

    starting_pos = fig.Position;

    width = 760;
    height = 600 / options.aspect_ratio_multiplier;

    try
        starting_title = copy(get(gca, 'title'));

        if ~options.keep_title
            title('');
        end

        fig.Position = [0 0 width height];
        exportgraphics(fig, strcat("dist/", filename));

        set(gca, 'title', starting_title);
        fig.Position = starting_pos;

    catch
        fig.Position = [0 0 width height];
        exportgraphics(fig, strcat("dist/", filename));

        fig.Position = starting_pos;

    end

end
