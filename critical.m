% Read files
files = dir('stat_specific/*.stat');  % Only get files with the .stat extension
[~, idx] = sort({files.name});  % Sort based on file name
files = files(idx);  % Reorder files

cumulative_contour_values = [];
total_contour_values = [];
z_max_list = [];
top_filled_z_list = [];
avg_x_list = [];
z_fluct_list = [];
file_colors = []; 

% Loop through each file
for i = 1:length(files)
    disp(i)
    file = files(i);

    % Determine the label 
    file_label = file.name(1);

    % Get the color for the current label
    color = getColorForLabel(file_label);
    file_colors = [file_colors; color]; 

    % Get the full file path
    file_path = fullfile(file.folder, file.name);
    data = readMercuryCG(file_path);

    % Determine the top-filled height
    all_z_values = unique(data.z);
    particle_diameter = 1;
    zero_value_z_indices = find_zero_value_layers(data, all_z_values);

    if ~isempty(zero_value_z_indices)
        first_zero_z = all_z_values(zero_value_z_indices(1));
        effective_top_z = first_zero_z - (5 * particle_diameter);
        [~, idx_closest_top_z] = min(abs(all_z_values - effective_top_z));
        top_filled_z = all_z_values(idx_closest_top_z);
    else
        top_filled_z = max(all_z_values);
    end

    top_filled_z_list = [top_filled_z_list; top_filled_z];
    z_max_list = [z_max_list; top_filled_z];

    % Process zero crossings
    z_start = 0.25;
    z_step = 1;
    x_values = unique(data.x(:));
    zero_crossing_x = zeros(1, floor(top_filled_z) + 1);
    z_val_list = zeros(1, floor(top_filled_z) + 1);
    tolerance = 1e-3;
    step = 0;

    for z_val = z_start:z_step:top_filled_z
        step = step + 1;
        indices = abs(data.z - z_val) < tolerance;
        velocityY_values = data.VelocityY(indices);
        x_values_filtered = data.x(indices);

        if isempty(velocityY_values) || length(velocityY_values) < 2
            continue;
        end
        
        % Check for zero crossing
        zero_crossing_found = false;
        for k = 1:(length(velocityY_values) - 1)
            
            if velocityY_values(k) < 0 && velocityY_values(k + 1) > 0
                x1 = x_values_filtered(k);
                x2 = x_values_filtered(k + 1);
                y1 = velocityY_values(k);
                y2 = velocityY_values(k + 1);

                zero_crossing_x(step) = x1 - y1 * (x2 - x1) / (y2 - y1);
                z_val_list(step) = z_val;
                zero_crossing_found = true;
                break;
            end

            % Additional check
            if velocityY_values(k) == 0 && velocityY_values(k + 1) > 0
                zero_crossing_x(step) = x_values_filtered(k);
                z_val_list(step) = z_val;
                zero_crossing_found = true;
                break;
            end
        end

        if ~zero_crossing_found
            break;
        end
    end

    non_zero_indices = find(z_val_list);
    zero_crossing_x = zero_crossing_x(1:length(non_zero_indices));
    z_val_list = z_val_list(1:length(non_zero_indices));

    % Define the threshold for avg_x based on the label
    if ismember(file_label, 'A':'M')
        threshold = 70;
    else ismember(file_label, 'N':'S')
        threshold = 140;
    end

    % Average x-value calculation, accounting for points below the z-threshold
    high_z_indices = z_val_list >= threshold;
    
    if any(high_z_indices)
        % Calculate the average for values above the z-threshold
        avg_x = mean(zero_crossing_x(high_z_indices));
    else
        % If no values are above the threshold, use the average of all x-positions
        avg_x = mean(zero_crossing_x);
    end
    
    % Append the avg_x value to the list
    avg_x_list = [avg_x_list; avg_x];
    
    % First height where x exceeds avg_x
    if avg_x < 0
        idx_first_exceeds_avg_x = find(zero_crossing_x < avg_x, 1);
    else
        idx_first_exceeds_avg_x = find(zero_crossing_x > avg_x, 1);
    end
    
    if isempty(idx_first_exceeds_avg_x)
        height_at_avg_x = NaN;
    else
        height_at_avg_x = z_val_list(idx_first_exceeds_avg_x);
    end
    
    % Append height fluctuation to the list
    z_fluct_list = [z_fluct_list; height_at_avg_x];

end

% Define colors explicitly for legend and points
blue_color = [0, 0.4470, 0.7410];      % Blue (Width 30)
green_color = [0.4660, 0.6740, 0.1880]; % Green (Width 60)

% Add dummy points specifically for the legend
h1 = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', blue_color, 'MarkerFaceColor', blue_color); % Blue dummy
h2 = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', green_color, 'MarkerFaceColor', green_color); % Green dummy

% Plot: Average x-position when z > threshold
valid_indices = ~isnan(avg_x_list); % Include all valid entries
filtered_avg_x_list = abs(avg_x_list(valid_indices)); % Filter NaN values
filtered_z_max_list = z_max_list(valid_indices); % Corresponding z_max values
filtered_colors = file_colors(valid_indices, :); % Corresponding colors

% Define y-threshold
y_threshold_min = 9;

% Separate points above and below the threshold
above_threshold_indices = filtered_avg_x_list >= y_threshold_min;
below_threshold_indices = ~above_threshold_indices; % Opposite of above_threshold_indices

% Separate data for above and below threshold
avg_x_above = filtered_avg_x_list(above_threshold_indices);
z_max_above = filtered_z_max_list(above_threshold_indices);
colors_above = filtered_colors(above_threshold_indices, :);

avg_x_below = filtered_avg_x_list(below_threshold_indices);
z_max_below = filtered_z_max_list(below_threshold_indices);
colors_below = filtered_colors(below_threshold_indices, :);

figure;
hold on;
 
%Plot points above the threshold (filled circles)
for i = 1:length(z_max_above)
    plot(z_max_above(i), avg_x_above(i), 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', colors_above(i, :), ...
    'MarkerFaceColor', colors_above(i, :)); % Filled circle
end

% Plot points below the threshold (hollow circles)
for i = 1:length(z_max_below)
    plot(z_max_below(i), avg_x_below(i), 'o', 'MarkerSize', 8, ...
         'MarkerEdgeColor', colors_below(i, :), ...
         'MarkerFaceColor', 'none'); % Hollow circle
end

% Calculate and plot the average line for points above the threshold
mean_avg_x = mean(avg_x_above);
yline(mean_avg_x, 'r--', 'LineWidth', 1.5); % Red dotted line

% Add a text label for the average line
text(mean(filtered_z_max_list) - 50, mean_avg_x + 0.2, sprintf('$x_r = %.2f$', mean_avg_x), ...
    'Interpreter', 'latex', ...
    'FontSize', 12, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');


% Dummy plots for legend
dummy_filled = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [0, 0.4470, 0.7410], ...
    'MarkerFaceColor', [0, 0.4470, 0.7410]); % Filled blue
dummy_hollow = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [0, 0.4470, 0.7410], ...
    'MarkerFaceColor', 'none'); % Hollow blue

legend([dummy_filled, dummy_hollow], {'$<x>$ above $z_r$', '$x$ below $z_r$'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 10, ...
    'Location', 'northeast');

% Add labels and title
xlabel('$H(d)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x/W$', 'Interpreter', 'latex', 'FontSize', 12);
%title('Average x-position vs Height', 'Interpreter', 'latex', 'FontSize', 14);

grid on;

% Save the figure
saveas(gcf, 'figures/avg_x_plot_all.png');

% Plot: z-position when average x is reached
valid_indices = ~isnan(z_fluct_list);
filtered_z_fluct_list = z_fluct_list(valid_indices);
filtered_z_max_list_for_fluct = z_max_list(valid_indices);
filtered_colors_for_fluct = file_colors(valid_indices, :);

% Sort the points by x-values for line connectivity
[sorted_z_max, sort_idx] = sort(filtered_z_max_list_for_fluct);
sorted_z_fluct = filtered_z_fluct_list(sort_idx);

% Identify indices where sorted_z_max exceeds 80
valid_indices = sorted_z_max > 80;

% Plot the connected line only if at least two points exceed the threshold
figure;
hold on;
if sum(valid_indices) > 1
    plot(sorted_z_max(valid_indices), sorted_z_fluct(valid_indices), ':', ...
         'LineWidth', 1.5, 'Color', [0 0 0]); % Black dotted line
end

% Plot the individual points on top (only those exceeding 80)
for i = 1:length(sorted_z_fluct)
    if sorted_z_max(i) > 80
        plot(sorted_z_max(i), sorted_z_fluct(i), 'o', 'MarkerSize', 8, ...
             'MarkerEdgeColor', filtered_colors_for_fluct(sort_idx(i), :), ...
             'MarkerFaceColor', filtered_colors_for_fluct(sort_idx(i), :));
    end
end
hold off;


% Add labels, title, and grid
xlabel('$H(d)$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$z_r(d)$','Interpreter', 'latex', 'FontSize', 12);
grid on;

% Save the figure
saveas(gcf, 'figures/avgx_z_connected_plot.png');

% Helper function: Assign colors based on file label
function color = getColorForLabel(label)
    if ismember(label, 'A':'M') % Labels A-M -> Width 30 (blue)
        color = [0, 0.4470, 0.7410];
    else ismember(label, 'N':'S') % Labels R-W -> Width 60 (green)
        color = [0.4660, 0.6740, 0.1880];
    end
end

%Find zero value layers
function zero_value_z_indices = find_zero_value_layers(data, all_z_values)
    zero_value_z_indices = [];
    for i = 1:length(all_z_values)
        z_val = all_z_values(i);
        indices = data.z == z_val;
        if all(data.ParticleSize1(indices) == 0)
            zero_value_z_indices(end+1) = i;
        end
    end
end








