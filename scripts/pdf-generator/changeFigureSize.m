function changeFigureSize(value)
% This function takes in a two element vector that sets the paper size and
% paper position based on those two values.
% J. Fritzinger, updated 2/23/2021

if ~isempty(value)
    fig = gcf;
    fig.PaperSize = value;
    fig.PaperPosition = [0 0 value];
    fig.Units = 'inches';
    fig.Position(3:4) = value;
end

end