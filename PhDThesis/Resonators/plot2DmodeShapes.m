%{
    Requires phi from twoDWave.m
%}

close all;
width = 0.28;
inc = 0.33;
start = 0.03
height = 0.4
figure('Position', [173 364 852 434])
titles = {"$(1, 1)$", "$(2, 1)$", "$(1, 2)$", "$(3, 1)$", "$(2, 2)$", "$(4, 1)$"};
for fig = 1:6
    subplot(2, 3, fig)
    data = reshape(phi(:,end-fig+1), Nyu, Nxu);
    data = data * sign(data(floor(Nyu/2), floor(Nxu/2)));
    [X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
    [X2,Y2] = meshgrid(1:0.01:size(data,2), 1:0.01:size(data,1));
    datainterp = interp2(X, Y, data, X2, Y2, 'cubic');
    imagesc(datainterp)
    xlabel('$x$', 'interpreter', 'latex', 'Fontsize', 16)
    ylabel('$y$', 'interpreter', 'latex', 'Fontsize', 16)
    title(titles{fig}, 'interpreter', 'latex', 'Fontsize', 16)
    colormap(gray);
    xticks([])
    yticks([])
    if fig < 4
        set(gca, 'Position', [start+(fig-1)*inc 0.5500 width height])
    else
        set(gca, 'Position', [start+(fig-4)*inc 0.0500 width height])
    end
end
