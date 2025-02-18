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
file_colors = []; % Store color corresponding to each file
oscillation_contribution_list = []; % Store oscillation contributions

% Loop through each file
for i = 1:length(files)
    file = files(i);

    % Determine the label (e.g., A, B, ..., W)
    file_label = file.name(1); % Assuming the label is the first character of the file name

    % Get the color for the current label
    color = getColorForLabel(file_label);
    file_colors = [file_colors; color]; % Save color for this file

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

        zero_crossing_found = false;
        for k = 1:(length(velocityY_values) - 1)
            % Check for zero crossing (negative to positive)
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

            % Additional check for zero transition
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

    % Average x-value at z > threshold
    high_z_indices = z_val_list >= threshold;
    avg_x = mean(zero_crossing_x(high_z_indices));
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
    z_fluct_list = [z_fluct_list; height_at_avg_x];

    % Adjusted contour value calculation
    if ~isempty(zero_crossing_x) && ~isempty(z_val_list)
        first_contour_value = sqrt((zero_crossing_x(1))^2 + (z_val_list(1))^2) - z_start;
        cumulative_contour_value = first_contour_value;
        total_contour_value = first_contour_value;
        oscillation_contribution = 0; % Initialize contribution value

        for j = 1:(length(zero_crossing_x) - 1)
            contour_value = sqrt((zero_crossing_x(j) - zero_crossing_x(j + 1))^2 + (z_step)^2);
            total_contour_value = total_contour_value + contour_value;
            contour_value_adjusted = contour_value - z_step;
            cumulative_contour_value = cumulative_contour_value + contour_value_adjusted;
        
        
            % Check if within the oscillation region (z > height_at_avg_x)
            if ~isnan(height_at_avg_x) && z_val_list(j) > height_at_avg_x
                oscillation_contribution = oscillation_contribution + contour_value_adjusted;
            end
        end


        total_contour_values = [total_contour_values; total_contour_value];
        cumulative_contour_values = [cumulative_contour_values; cumulative_contour_value];
        oscillation_contribution_list = [oscillation_contribution_list; oscillation_contribution];
    else
        total_contour_values = [total_contour_values; NaN];
        cumulative_contour_values = [cumulative_contour_values; NaN];
        oscillation_contribution_list = [oscillation_contribution_list; 0];
    end
end


disp(size(oscillation_contribution_list));
% Loop through each file and calculate oscillation contribution
for i = 1:length(files)
    file_label = files(i).name(1); % Extract first letter of filename
    
   
    
    % Apply division rule
    if ismember(file_label, 'A':'M')
        oscillation_contribution_list(i) = oscillation_contribution_list(i) / 30;
    elseif ismember(file_label, 'N':'S')
        oscillation_contribution_list(i) = oscillation_contribution_list(i) / 60;
    else
        error('Unexpected file label: %s', file_label);
    end
end

disp('oscillation')
disp(oscillation_contribution_list)

% Plot: Adjusted contour values
figure;
hold on;

% Define colors explicitly for legend and points
blue_color = [0, 0.4470, 0.7410];      % Blue (Width 30)
green_color = [0.4660, 0.6740, 0.1880]; % Green (Width 60)

% Plot actual data points
for i = 1:length(z_max_list)
    if i < 11
        cumulative_contour_values(i) = cumulative_contour_values(i)/30;
    else
        cumulative_contour_values(i) = cumulative_contour_values(i)/60;
    end
    plot(z_max_list(i), cumulative_contour_values(i), 'o', 'MarkerSize', 8, ...
         'MarkerEdgeColor', file_colors(i, :), ...
         'MarkerFaceColor', file_colors(i, :));

    % Plot oscillation contribution points and connect with a vertical line
    if ~isnan(cumulative_contour_values(i)) && oscillation_contribution_list(i) > 0
        % Calculate the lower point (empty marker)
        lower_point = cumulative_contour_values(i) - oscillation_contribution_list(i);

        % Draw a vertical dotted line
        plot([z_max_list(i), z_max_list(i)], [lower_point, cumulative_contour_values(i)], ':', 'Color', file_colors(i, :), 'LineWidth', 1.5);

        % Plot the lower point as an empty marker
        plot(z_max_list(i), lower_point, 'o', 'MarkerSize', 8, ...
             'MarkerEdgeColor', file_colors(i, :), ...
             'MarkerFaceColor', 'none'); % Empty marker with blue edge
    end
end

% Add dummy points specifically for the legend
h1 = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', blue_color, 'MarkerFaceColor', blue_color); % Blue dummy
h2 = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', green_color, 'MarkerFaceColor', green_color); % Green dummy

% Add legend explicitly with handles
legend([h1, h2], {'Width 30', 'Width 60'}, 'Location', 'NorthWest', 'FontSize', 10);

% Add labels and title
xlabel('$H(d)$', 'Interpreter', 'latex');
ylabel('$C_{a}(d)/W(d)$', 'Interpreter', 'latex');
grid on;

% Save the figure
saveas(gcf, 'figures/cumulative_contour_adjusted.png');

% Helper function: Assign colors based on file label
function color = getColorForLabel(label)
    if ismember(label, 'A':'M') % Labels A-M -> Width 30 (blue)
        color = [0, 0.4470, 0.7410];
    else ismember(label, 'N':'S') % Labels R-W -> Width 60 (green)
        color = [0.4660, 0.6740, 0.1880];
    end
end

% Helper function: Find zero value layers
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

% The oscillation_contribution_list now contains the additional contour length contribution for each dataset.
