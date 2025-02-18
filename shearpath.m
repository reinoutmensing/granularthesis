% Retrieve all files with the .stat extension and sort them by name
files = dir('stat_specific/*.stat');
[~, idx] = sort({files.name});
files = files(idx);

% Initialize storage arrays for computed values
cumulative_contour_values = [];
normalized_contour_values = [];
total_contour_values = [];
z_files = [];
top_filled_z_list = [];

avg_x_list = [];
z_fluct_list = [];

% Define subplot grid dimensions
num_rows = 2;
num_cols = 5;

% Create a figure for visualization
figure;
z_max_list = [];

% Process each file in the dataset
for i = 1:num_files
    file = files(i);
    file_path = fullfile(file.folder, file.name);
    data = readMercuryCG(file_path);
    
    % Identify the top-filled height within the system
    all_z_values = unique(data.z);
    particle_diameter = 1;
    zero_value_z_indices = find_zero_value_layers(data, all_z_values);
    
    % Determine the effective top z-value for valid data
    if ~isempty(zero_value_z_indices)
        first_zero_z = all_z_values(zero_value_z_indices(1));
        effective_top_z = first_zero_z - (3 * particle_diameter);
        [~, idx_closest_top_z] = min(abs(all_z_values - effective_top_z));
        top_filled_z = all_z_values(idx_closest_top_z);
    else
        top_filled_z = max(all_z_values); % Default to the highest z-value if no zero layers exist
    end

    % Store computed top-filled height values
    top_filled_z_list = [top_filled_z_list; top_filled_z];
    z_max_list = [z_max_list; top_filled_z];

    % Define z-scan parameters
    z_start = 0.25;
    z_step = 1;
    
    % Extract unique x and z values from the dataset
    x_values = unique(data.x(:));
    z_values = unique(data.z(:));

    % Initialize arrays for zero-crossing computations
    zero_crossing_x = zeros(1, floor(top_filled_z) + 1);
    z_val_list = zeros(1, floor(top_filled_z) + 1);
    tolerance = 1e-3;
    
    step = 0; % Counter for z-value iterations

    % Scan through z-values to detect zero crossings
    for z_val = z_start:z_step:top_filled_z
        step = step + 1;
        indices = abs(data.z - z_val) < tolerance;
        velocityY_values = data.VelocityY(indices);

        zero_crossing_found = false; % Flag for detecting zero crossings

        % Identify zero crossing in velocity field
        for k = 1:(length(x_values) - 1)
            if velocityY_values(k) < 0 && velocityY_values(k + 1) > 0
                x1 = x_values(k);
                x2 = x_values(k + 1);
                y1 = velocityY_values(k);
                y2 = velocityY_values(k + 1);

                % Compute interpolated x-location of zero crossing
                zero_crossing_x(step) = x1 - y1 * (x2 - x1) / (y2 - y1);
                z_val_list(step) = z_val;
                zero_crossing_found = true;
                break;
            end
        end

        % Stop scanning if no zero crossing is found
        if ~zero_crossing_found
            break;
        end
    end

    % Trim zero values from computed arrays
    non_zero_indices = find(z_val_list);
    zero_crossing_x = zero_crossing_x(1:length(non_zero_indices));
    z_val_list = z_val_list(1:length(non_zero_indices));

    % Compute average x-value at high z levels
    high_z_indices = z_val_list >= 75;
    avg_x = mean(zero_crossing_x(high_z_indices));
    avg_x_list = [avg_x_list; avg_x];

    % Identify first height where x exceeds avg_x
    if avg_x < 0
        idx_first_exceeds_avg_x = find(zero_crossing_x < avg_x, 1);
    else
        idx_first_exceeds_avg_x = find(zero_crossing_x > avg_x, 1);
    end
    height_at_avg_x = z_val_list(idx_first_exceeds_avg_x);
    z_fluct_list = [z_fluct_list, height_at_avg_x];

    % Display computed results for the current file
    fprintf('For file %d:\n', i);
    fprintf('Average x position above height 70: %.2f\n', avg_x);
    fprintf('Height at which x first exceeds average: %.2f\n\n', height_at_avg_x);

    % Plot zero-crossing locations
    subplot(num_rows, num_cols, i);
    plot(zero_crossing_x, z_val_list, 'b-');
    hold on;
    plot(avg_x * ones(size(z_val_list)), z_val_list, 'r--');
    plot(zero_crossing_x(idx_first_exceeds_avg_x), height_at_avg_x, 'kx', 'MarkerSize', 8, 'LineWidth', 2);
    hold off;

    % Label axes conditionally for clarity
    if i == 1 || i == 6
        ylabel('$z(d)$', 'Interpreter', 'latex', 'FontSize', 12);
    end
    xlabel('$x(d)$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('$H = %.2f$', top_filled_z_list(i)), 'Interpreter', 'latex', 'FontSize', 14);

    xlim([-15, 15]);
    grid on;
    saveas(gcf, 'figures/shearband_path.png');
end


% Function to identify zero-value layers in the dataset
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
