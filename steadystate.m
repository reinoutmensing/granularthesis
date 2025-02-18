
figure;
data=readMercuryCG('CSCRun.stat'); 
plot(data.t,data.StressYY,data.t,data.StressXY)

% Save the figure
saveas(gcf, 'steady_state.png'); % Saves the current figure as a PNG file
close figure;