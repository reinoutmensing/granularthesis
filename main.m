%Main script

data = read_and_preprocess('stat_specific/B_CSCRun_30_10_50.stat'); % Read and preprocess the data

z_values = get_z_value(data);

plot_velocity_contour(data); % Plot the contour for VelocityY

plot_density_contour(data);

V = 1/40; % Velocity at boundaries

[best_fits1, best_fits2] = shear_rate_velocity_analysis(data, V); % Analyze shear rates and fits

plot_width_vs_height(z_values, best_fits2); % Plot Width vs Height

plot_log_log_width_vs_height(z_values, best_fits2); % Log-log plot of Width vs Height

plot_center_and_width(z_values, best_fits2); % Plot the center and width of the shear band

% read_and_preprocess.m
function data = read_and_preprocess(filename)
    data = readMercuryCG(filename);
    % Any other preprocessing steps if needed can be added here
end


% plot_velocity_contour.m
function plot_velocity_contour(data)
    figure;
    contourf(data.x, data.z, data.VelocityY, linspace(- 0.0250, 0.0250), 'EdgeColor', 'none');
    colormap(jet);
    colorbar;
    title('VelocityY Contour');
    saveas(gcf, "figures/Velocity_profile.png");
end

% plot_velocity_contour.m
function plot_density_contour(data)
    figure;
    contourf(data.x, data.z, data.Density, linspace(0, 1),'EdgeColor', 'none');
    colormap(jet);
    colorbar;
    title('Density Contour');
    saveas(gcf, "figures/density_profile.png");
end

function [best_fits1, best_fits2] = shear_rate_velocity_analysis(data, V)
    initial_w = linspace(0, 5, 20);
    initial_x0 = linspace(min(data.x(:)), max(data.x(:)), 20);
    
    ContactStressZZ = get_rid_of_nan(data.ContactStressZZ);
    ContactStressXZ = get_rid_of_nan(data.ContactStressXZ);
    ContactStressXY = get_rid_of_nan(data.ContactStressXY);

    z_values = get_z_value(data);
    particle_diameter = 1;
    
    % Initialize storage for fit parameters
    best_fits1 = zeros(length(z_values), 4);
    best_fits2 = zeros(length(z_values), 3);

    % Initialize figures for 3x3 subplots
    figure('Name', 'Shear Rate in YX Direction'); % Shear Rate plot
    figure('Name', 'Shear Stress in YX Direction'); % Shear Stress plot
    figure('Name', 'Proportionality Between Shear Rate and Shear Stress'); % Proportionality plot
    figure('Name', 'VelocityY Plot'); % VelocityY plot
    figure('Name','Mu Plot');
    figure('Name','Inertial Plot');
    figure('Name', 'Shear Velocity vs Shear Rate'); % Shear Velocity vs Shear Rate plot
    figure('Name', 'Mu_together');
    figure('Name', 'Inertial together');
    figure('Name','Mu versus Inertial');
    figure("Name",'Mu VS log inertial');
    figure("Name",'Zoomed in muI');
    % Initialize cell arrays to store values for each height
    all_mu = {};
    all_inertial = {};
    color_map = lines(length(z_values)); % Define colors for each height

    for i = 1:length(z_values)
        z_val = z_values(i);
        indices = data.z == z_val;
        
        x_values = data.x(indices);
        velocityY_values = data.VelocityY(indices);
        shear_rate_yx = diff(velocityY_values) ./ diff(x_values);

         % Calculate the inertial number I
        pressure = mean(data.ContactStressZZ(indices)); % Assuming constant pressure in z layer
        density = mean(data.Density(indices)); % Assuming constant density in z layer
        inertial_number = mean(shear_rate_yx) * particle_diameter ./ sqrt(pressure / density);
        
        % Inertial number per point
        % Extrapolate the last value based on the trend of the last two values
        last_value = shear_rate_yx(end) + (shear_rate_yx(end) - shear_rate_yx(end-1));

        % Append the extrapolated value to the original array
        shear_rate_yx_extrapolated = [shear_rate_yx; last_value];

        % Now use the extrapolated array for further calculations
        inertial_per_point = shear_rate_yx_extrapolated .* particle_diameter ./ sqrt(data.ContactStressZZ(indices) ./ data.Density(indices));
        
        % Store mu and inertial number for each height
        mu = data.ContactStressXY(indices)./data.ContactStressZZ(indices); % Calculate mu
        all_mu{i} = mu;
        all_inertial{i} = inertial_per_point;

        % Fit error function model for velocityY
        [bestParams1, maxRsq1] = fit_error_function(x_values, velocityY_values, V, initial_w, initial_x0);
        width = bestParams1(1);
        best_fits1(i, :) = [bestParams1, maxRsq1, width];
        
        % Fit Gaussian model for shear rate
        [bestParams2, maxRsq2] = fit_gaussian_model(x_values(1:end-1), shear_rate_yx, initial_w);
        bestParams2(3) = abs(bestParams2(3)); % Ensure width is positive
        best_fits2(i, :) = [bestParams2(3), maxRsq2, bestParams2(2)];
        
        % Plot shear rate in the yx direction with Gaussian fit
        figure(1);
        subplot(3, 3, i);
        scatter(x_values(1:end-1), shear_rate_yx, 'bo'); % Shear rate data points
        hold on;
        fittedX = linspace(min(x_values), max(x_values), 100);
        fittedRateY = bestParams2(1) * exp(-((fittedX - bestParams2(2)).^2) / (2 * bestParams2(3)^2));
        plot(fittedX, fittedRateY, 'r-'); % Fitted Gaussian curve
        xlabel('X Position');
        ylabel('Shear Rate (\gamma_{yx})');
        title(sprintf('z = %.3f', z_val));
        grid on;
        hold off;
        
        % Plot shear stress in the yx direction
        figure(2);
        subplot(3, 3, i);
        scatter(x_values, data.ContactStressXY(indices)./data.ContactStressZZ(indices), 'go'); % Shear stress data points
        xlabel('X Position');
        ylabel('Shear Stress (\tau_{yx})');
        title(sprintf('z = %.3f', z_val));
        grid on;
        
        % Ensure shear rate and stress arrays match in size
        shear_rate_yx = shear_rate_yx';
        shear_stress_xy = data.StressXY(indices);
        shear_stress_xy = shear_stress_xy(1:end-1); % Reduce to match shear_rate_yx size
        
        % Linear fit for proportionality and calculation of slope and R^2
        coeffs = polyfit(shear_rate_yx, shear_stress_xy, 1);
        fitted_values = polyval(coeffs, shear_rate_yx);
        residuals = shear_stress_xy - fitted_values;
        SSresid = sum(residuals.^2);
        SStotal = (length(shear_stress_xy) - 1) * var(shear_stress_xy);
        Rsq = 1 - SSresid / SStotal;
        slope = coeffs(1);
        
        % Plot shear rate vs shear stress for proportionality
        figure(3);
        subplot(3, 3, i);
        scatter(shear_rate_yx, shear_stress_xy, 'b*'); % Shear rate vs. shear stress
        hold on;
        plot(shear_rate_yx, fitted_values, 'r-'); % Linear fit line
        xlabel('Shear Rate (\gamma_{yx})');
        ylabel('Shear Stress (\tau_{yx})');
        grid on;
        hold off;
        
        % Plot VelocityY as a regular line plot for each z value
        figure(4);
        subplot(3, 3, i);
        plot(x_values, velocityY_values, 'b-'); % Plot X vs VelocityY
        xlabel('X Position');
        ylabel('VelocityY');
        title(sprintf('VelocityY at z = %.3f', z_val));
        grid on;
        
        % Plot VelocityY as a regular line plot for each z value
        figure(5);
        subplot(3, 3, i);
        plot(x_values, mu, 'b-'); % Plot X vs VelocityY
        xlabel('X Position');
        ylabel('mu');
        title(sprintf('Mu at z = %.3f', z_val));
        grid on;

        % Plot VelocityY as a regular line plot for each z value
        figure(6);
        subplot(3, 3, i);
        plot(x_values, inertial_per_point, 'b-'); % Plot X vs VelocityY
        xlabel('X Position');
        ylabel('Inertial number');
        title(sprintf('Inertial number at z = %.3f', z_val));
        grid on;
        % Plot Shear Velocity vs Shear Rate
     
    shear_velocity_yx = velocityY_values(2:end) - velocityY_values(1:end-1);
    figure(7);
    subplot(3, 3, i);
    scatter(shear_velocity_yx, shear_rate_yx, 'b*'); % Shear velocity vs. shear rate
    hold on;
    
    % Perform linear fit and calculate the slope
    fit_sv_sr = polyfit(shear_velocity_yx, shear_rate_yx, 1);
    slope_sv_sr = fit_sv_sr(1); % The slope of the fit
    
    % Plot the linear fit line
    plot(shear_velocity_yx, polyval(fit_sv_sr, shear_velocity_yx), 'r-'); % Linear fit
    xlabel('Shear Velocity Difference');
    ylabel('Shear Rate (\gamma_{yx})');
    title(sprintf('z = %.3f', z_val));
    sgtitle(sprintf('Slope: %.4f', slope_sv_sr))
    grid on;
    hold off;

    end

    markers = {'o', 's', 'd', '^', 'v', '>', '<', '*', 'x', '+'};  
        % Plot all mu values together in one plot
    figure(8)
    hold on;
    for j = 1:length(z_values)
        plot(x_values, all_mu{j}, 'Color', color_map(j, :), 'DisplayName', sprintf('z = %.3f', z_values(j)));
    end
    xlabel('Width');
    ylabel('Friction Coefficient (\mu)');
    legend('show');
    title('Friction Coefficient (\mu) for All Heights');
    grid on;
    hold off;
    

    % Plot all I values together in one plot
    figure(9)
    hold on;
    for j = 1:length(z_values)
        plot(x_values, all_inertial{j}, 'Color', color_map(j, :), 'DisplayName', sprintf('z = %.3f', z_values(j)));
    end
    xlabel('Width');
    ylabel('Inertial Number (I)');
    legend('show');
    title('Inertial Number (I) for All Heights');
    grid on;
    hold off;
    
    % Plot mu vs. Inertial number (I) for each height in a new figure
    figure(10);
    hold on;
    for j = 1:length(z_values)
        % Plot mu vs. I
        plot(all_inertial{j}, abs(all_mu{j}), 'o-', 'Color', color_map(j, :), 'DisplayName', sprintf('z = %.3f', z_values(j)));
    end
    xlabel('$I$', 'Interpreter', 'latex', 'FontSize', 16);  % Ensure correct LaTeX syntax for log
    ylabel('$\mu_{macro}$', 'Interpreter', 'latex', 'FontSize', 16);
    legend('show');
    title('$\mathbf{CG\ width\ 3}$', 'Interpreter', 'latex', 'FontSize', 16); 
    hold off;
    
    

    % Plot mu vs. Inertial number (I) for each height in a new figure
    figure(11);
    clf;
    hold on;
    for j = 1:length(z_values)
        % Use both distinct colors and markers
        plot(all_inertial{j}, abs(all_mu{j}), 'Color', color_map(j, :), 'Marker', markers{mod(j-1, length(markers))+1}, 'LineStyle', 'none', 'DisplayName', sprintf('z = %.3f', z_values(j)), 'MarkerSize', 6);
    end
    set(gca, 'XScale', 'log');  % Set the x-axis to a logarithmic scale
    
    xlabel('$I$', 'Interpreter', 'latex', 'FontSize', 16);  % Ensure correct LaTeX syntax for log
    ylabel('$\mu_{macro}$', 'Interpreter', 'latex', 'FontSize', 16);  % Correct LaTeX syntax for mu
    legend('show', 'Interpreter', 'latex', 'Location', 'northwest','FontSize', 12);  % Set legend to LaTeX mode and position it top left
    title('$\mathbf{CG\ width\ 3}$', 'Interpreter', 'latex', 'FontSize', 16);  % Ensure consistent font size for title
    grid on;
    hold off;



    % Define different marker shapes and a colormap
    markers = {'o', 's', 'd', '^', 'v', '>', '<', '*', 'x', '+'};  % Circle, square, diamond, triangle up/down/left/right, star, cross, plus
    
    figure(12);
    clf;
    hold on;
    for j = 1:length(z_values)
        % Use both distinct colors and markers
        plot(all_inertial{j}, abs(all_mu{j}), 'Color', color_map(j, :), 'Marker', markers{mod(j-1, length(markers))+1}, 'LineStyle', 'none', 'DisplayName', sprintf('z = %.3f', z_values(j)), 'MarkerSize', 6);
    end
    %set(gca, 'XScale', 'log');  % Set the x-axis to a logarithmic scale
    xlim([0, 0.08]);  % Adjust x-axis limits to focus on relevant region
    ylim([0.25,0.45]);
    xlabel('Inertial Number (I)');
    ylabel('Friction Coefficient (\mu)');
    legend('show', 'Location', 'northwest');  % Set legend location to top left
    title('Friction Coefficient (\mu) vs. Inertial Number (I) for All Heights');
    grid on;
    hold off;
        
    % Save figures
    saveas(figure(1), 'figures/shear_rate_yx.png');
    saveas(figure(2), 'figures/shear_stress_yx.png');
    saveas(figure(3), 'figures/proportionality_yx.png');
    saveas(figure(4), 'figures/velocity_y.png');
    saveas(figure(5), 'figures/mu.png');
    saveas(figure(6), 'figures/Inertial.png');
    saveas(figure(7), 'figures/shear_velocity_vs_shear_rate.png');
    saveas(figure(8), 'figures/all_mu.png'); % Save the mu plot
    saveas(figure(9), 'figures/all_inertial.png'); % Save the inertial number plot
    saveas(figure(10), 'figures/inertialvsmunieuwste.png');
    saveas(figure(11), 'figures/loginertialvsmunieuwste.png');
    saveas(figure(12), 'figures/zoominMuI.png');
    % Display results
    disp('BESTFITS1');
    disp(array2table(best_fits1, 'VariableNames', {'BestW', 'BestX0', 'BestRsq', 'Width'}, 'RowNames', string(z_values)));
    disp('BESTFITS2');
    disp(array2table(best_fits2, 'VariableNames', {'BestW', 'BestRsq', 'Center'}, 'RowNames', string(z_values)));
end

% fit_error_function.m
function [bestParams, maxRsq] = fit_error_function(x, y, V, initial_w, initial_x0)
    model = @(params, x) (V/2) * erf((x - params(2)) / (sqrt(2) * params(1)));
    maxRsq = -inf;
    bestParams = [0, 0];
    
    for w = initial_w
        for x0 = initial_x0
            params = [w, x0];
            [beta, R] = nlinfit(x, y, model, params);
            rsq = calculate_rsq(R, y);
            if rsq > maxRsq
                maxRsq = rsq;
                bestParams = beta;
            end
        end
    end
end

% fit_gaussian_model.m
function [bestParams, maxRsq] = fit_gaussian_model(x, y, initial_w)
    model = @(params, x) params(1) * exp(-((x - params(2)).^2) / (2 * params(3)^2));
    maxRsq = -inf;
    bestParams = [0, 0, 0];
    
    for w = initial_w
        params = [max(y), mean(x), w];
        [beta, R] = nlinfit(x, y, model, params);
        rsq = calculate_rsq(R, y);
        if rsq > maxRsq
            maxRsq = rsq;
            bestParams = beta;
        end
    end
end

% plot_width_vs_height.m
function plot_width_vs_height(z_values, best_fits1)
    figure;
    plot(z_values,best_fits1(:,1), 'bo-');
    ylabel('Width (W)');
    xlabel('Height');
    title('Width vs. Height');
    grid on;
    saveas(gcf, 'figures/shear_band_width_vs_height.png');
end

% plot_log_log_width_vs_height.m
function plot_log_log_width_vs_height(z_values, best_fits1)
    figure;
    loglog(z_values, best_fits1(:,1), '-o');
    xlabel('log(Height)');
    ylabel('log(Width)');
    title('Log-Log Plot of Width vs. Height');
    grid on;
    saveas(gcf, 'figures/log_log_plot.png');
end

% plot_center_and_width.m
function plot_center_and_width(z_values, best_fits1)
    figure;
    hold on;
    centers = best_fits1(:,2);
    widths = best_fits1(:,1);

    for i = 1:length(z_values)
        center = centers(i);
        width = widths(i);
        plot([center - width/2, center + width/2], [z_values(i), z_values(i)], 'b-o');
    end

    xlabel('X Position');
    ylabel('Height');
    title('Shear Band Center and Width at Different Heights');
    xlim([-15, 15]);
    grid on;
    hold off;
    saveas(gcf, 'figures/shear_band_center_width.png');
end

% calculate_rsq.m
function rsq = calculate_rsq(residuals, actual)
    % Calculate Sum of Squares of Residuals
    SSresid = sum(residuals .^ 2);
    
    % Calculate Total Sum of Squares
    SStotal = (length(actual) - 1) * var(actual);
    
    % Compute R-squared
    rsq = 1 - SSresid / SStotal;
end

function z_values = get_z_value(data)

    % Get all unique z values from data.z
    all_z_values = unique(data.z);
    particle_diameter = 1; % Define the particle diameter, adjust this based on your system
    
    
    % Find the first zero-value z layer and deduct three particle diameters
    zero_value_z_indices = find_zero_value_layers(data, all_z_values);
    
    if ~isempty(zero_value_z_indices)
        first_zero_z = all_z_values(zero_value_z_indices(1));
        effective_top_z = first_zero_z - (3 * particle_diameter);
        
        % Find the closest actual z value in data.z to the effective top z
        [~, idx_closest_top_z] = min(abs(all_z_values - effective_top_z));
        top_filled_z = all_z_values(idx_closest_top_z);
        
        % Generate 9 equally spaced points from the nearest to 1 to top_filled_z
        nearest_to_one_idx = find(all_z_values >= 1, 1);  % Nearest point to 1 or above
        %REPLACED TOP_FILLED_Z WITH 41.5 SO EVERY COARSE GRAINING WIDTH HAS
        %THE SAME VALUE. 
        equally_spaced_z = linspace(all_z_values(nearest_to_one_idx), top_filled_z, 9);
        
        % Adjust to closest actual values in data.z
        z_values = arrayfun(@(z) all_z_values(find(abs(all_z_values - z) == min(abs(all_z_values - z)), 1)), equally_spaced_z);
        disp('Selected z_values for analysis:');
        disp(z_values);
    else
        disp('No layers with zero values were found.')
        return;
    end

end

function zero_value_z_indices = find_zero_value_layers(data, all_z_values)
    zero_value_z_indices = [];
    for i = 1:length(all_z_values)
        z_val = all_z_values(i);
        indices = data.z == z_val; % Logical indices for current z layer
        % Check if all data points in this z layer are zero
        if all(data.ParticleSize1(indices) == 0) % Replace 'ParticleSize1' with your field name
            zero_value_z_indices(end+1) = i;
        end
    end
end

function nonan = get_rid_of_nan(contact_stress)
% Getting rid of the Inf and the NaN values in ContactStress
    for k = 1:5  % Run this loop multiple times to propagate values across consecutive NaNs/Infs
        for i = 1:size(contact_stress, 1)
            for j = 2:size(contact_stress, 2)-1  % Ignore first and last columns initially
                if isnan(contact_stress(i, j)) || isinf(contact_stress(i, j))
                    % Check if previous and next elements are valid numbers
                    if ~isnan(contact_stress(i, j-1)) && ~isinf(contact_stress(i, j-1)) && ...
                       ~isnan(contact_stress(i, j+1)) && ~isinf(contact_stress(i, j+1))
                        % Calculate the average if both neighbors are valid
                        contact_stress(i, j) = (contact_stress(i, j-1) + contact_stress(i, j+1)) / 2;
                    elseif ~isnan(contact_stress(i, j-1)) && ~isinf(contact_stress(i, j-1))
                        % If only the previous value is valid, use it
                        contact_stress(i, j) = contact_stress(i, j-1);
                    elseif ~isnan(contact_stress(i, j+1)) && ~isinf(contact_stress(i, j+1))
                        % If only the next value is valid, use it
                        contact_stress(i, j) = contact_stress(i, j+1);
                    end
                end
            end
            
            % Handle NaNs or Infs at the edges (first and last columns)
            if isnan(contact_stress(i, 1)) || isinf(contact_stress(i, 1))
                contact_stress(i, 1) = contact_stress(i, 2);  % Replace with the next value
            end
            if isnan(contact_stress(i, end)) || isinf(contact_stress(i, end))
                contact_stress(i, end) = contact_stress(i, end-1);  % Replace with the previous value
            end
        end
    end
    nonan = contact_stress;
    return;
end
