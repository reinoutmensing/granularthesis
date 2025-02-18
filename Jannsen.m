% Read the data from MercuryCG
data = readMercuryCG('CSCRun.stat');

% Apply NaN handling to ContactStress
data.ContactStressZZ = get_rid_of_nan(data.ContactStressZZ);
data.ContactStressXX = get_rid_of_nan(data.ContactStressXX);
data.ContactStressYY = get_rid_of_nan(data.ContactStressYY);

figure;
contourf(data.x, data.z, data.ContactStressXX, linspace(0, 400), 'EdgeColor', 'none');
colormap(jet);
xlabel('x [pd]', 'FontName', 'Arial', 'FontSize', 12);
ylabel('z [pd]', 'FontName', 'Arial', 'FontSize', 12);
title('$\sigma_{xx}$', 'Interpreter', 'latex', 'FontSize', 14); % Correct LaTeX title
colorbar; % Adds a colorbar to the plot
saveas(gcf, 'figures/Velocity_profile_xx_400_vhigh.png');

% Plot the contour of ContactStressYY
figure;
contourf(data.x, data.z, data.ContactStressYY, linspace(0, 400), 'EdgeColor', 'none');
colormap(jet);
xlabel('x [pd]', 'FontName', 'Arial', 'FontSize', 12);
ylabel('z [pd]', 'FontName', 'Arial', 'FontSize', 12);
title('$\sigma_{yy}$', 'Interpreter', 'latex', 'FontSize', 14); % Correct LaTeX title
colorbar; % Adds a colorbar to the plot
saveas(gcf, 'figures/Velocity_profile_yy_400_vhighplop.png');

% Plot the contour of ContactStressZZ
figure;
contourf(data.x, data.z, data.ContactStressZZ, linspace(0, 400), 'EdgeColor', 'none');
colormap(jet);
xlabel('$x(d)$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 16);
ylabel('$z(d)$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 16);
%title('$\sigma_{zz}$', 'Interpreter', 'latex', 'FontSize', 14); % Use LaTeX for rendering
ylim([0,285]);
% Add colorbar with \sigma_{zz} as the description
c = colorbar;
c.Label.String = '$\sigma_{zz}(\rho g d)$'; % Set the colorbar description
c.Label.Interpreter = 'latex'; % Use LaTeX for proper rendering
c.Label.FontSize = 16; % Make the label slightly larger
c.Label.Rotation = 0; % Ensure the text is horizontal
c.Label.Position = [c.Label.Position(1) + 1, c.Label.Position(2)]; % Adjust label position

% Save the plot
saveas(gcf, 'figures/Velocity_profile_zz_400_slow.png');

% Unique z-values (height levels)
z_values = unique(data.z(:)); 

% Initialize array for averaged contact stress at each height
average_ContactStressZZ = zeros(length(z_values), 1);
average_ShearStressXY = zeros(length(z_values), 1);

% Loop over each unique z-value and compute the average ContactStressZZ
for i = 1:length(z_values)
    % Find all particles at the current height
    indices_at_height = abs(data.z - z_values(i)) < 1e-6; % Tolerance for floating-point comparison
    
    % Compute the average ContactStressZZ for this height
    average_ContactStressZZ(i) = mean(data.StressZZ(indices_at_height));
    average_ShearStressXY(i) = mean(data.StressXY(indices_at_height));
end

% Compute the average density at each height level
average_Density = zeros(length(z_values), 1);
for i = 1:length(z_values)
    % Find all particles at the current height
    indices_at_height = abs(data.z - z_values(i)) < 1e-6; % Tolerance for floating-point comparison
    % Compute the average Density for this height
    average_Density(i) = mean(data.Density(indices_at_height));
end

% Gravitational acceleration
g_z = -1;

% Compute sigma_L using cumulative integration
sigma_L = cumtrapz(z_values, average_Density * g_z); % Gravitational stress

% Display the minimum value of sigma_L
disp(['Minimum value of sigma_L (after offset): ', num2str(min(sigma_L))]);

% Compute wall traction (T_z)
gradient_sigma_zz = gradient(average_ContactStressZZ, z_values); % Gradient of ContactStressZZ

% Fill NaN values in average_ContactStressZZ with linear interpolation
average_ContactStressZZ = fillmissing(average_ContactStressZZ, 'linear', 'EndValues', 'nearest');

% Replace Inf values with interpolated values
average_ContactStressZZ(inf_indices) = interp1(z_values(~inf_indices), average_ContactStressZZ(~inf_indices), z_values(inf_indices), 'linear');

% Recompute wall traction with the filled average_ContactStressZZ
wall_traction = sigma_L - average_ContactStressZZ;

disp(['NaN in sigma_L: ', num2str(any(isnan(sigma_L)))]);
disp(['NaN in average_ContactStressZZ: ', num2str(any(isnan(average_ContactStressZZ)))]);
disp(['NaN in wall_traction: ', num2str(any(isnan(wall_traction)))]);

% Check for discontinuities in z_values
z_diff = diff(z_values);
disp(['Max difference in z_values: ', num2str(max(z_diff))]);
disp(['Min difference in z_values: ', num2str(min(z_diff))]);

% Check for discontinuities in average_ContactStressZZ
stress_diff = diff(average_ContactStressZZ);
disp(['Max difference in average_ContactStressZZ: ', num2str(max(stress_diff))]);
disp(['Min difference in average_ContactStressZZ: ', num2str(min(stress_diff))]);

% Plot sigma_L (gravitational stress), sigma_ZZ (measured stress), and T_z
figure;
hold on;
plot(z_values, sigma_L, 'r-', 'LineWidth', 2); % Plot sigma_L
plot(z_values, average_ContactStressZZ, 'b-', 'LineWidth', 2); % Plot average sigma_ZZ
plot(z_values, wall_traction, 'g-', 'LineWidth', 2); % Plot T_z
hold off;

% Add labels, legend, and title
xlabel('$z(d)$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 16);
ylabel('$\langle \sigma_{zz} \rangle_x \, (\rho g d)$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 16);
legend('$\sigma_L$', '$\sigma_{zz}$', '$T_z$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
grid on;

% Save the figure
saveas(gcf, 'figures/Stress_Comparison_slow.png');

disp(['Size of z_values: ', num2str(length(z_values))]);
disp(['Size of sigma_L: ', num2str(length(sigma_L))]);
disp(['Size of average_ContactStressZZ: ', num2str(length(average_ContactStressZZ))]);
disp(['Size of wall_traction: ', num2str(length(wall_traction))]);

% Plot the averaged contact stress vs height
figure;
plot(z_values, abs(average_ContactStressZZ), 'o-', 'LineWidth', 1.5);
xlabel('$z(d)$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 12);
ylabel('$\langle \sigma_{zz} \rangle_x$', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 12);
title('$\langle \sigma_{zz} \rangle$', 'Interpreter', 'latex', 'FontSize', 14); % Updated title with averaging notation
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.5);

% Save the figure
saveas(gcf, 'figures/averaged_sigma_zz_NEW.png');

% Function to handle NaN and Inf values
function nonan = get_rid_of_nan(contact_stress)
    % Getting rid of the Inf and the NaN values in ContactStress
    for k = 1:5 % Run this loop multiple times to propagate values across consecutive NaNs/Infs
        for i = 1:size(contact_stress, 1)
            for j = 2:size(contact_stress, 2)-1 % Ignore first and last columns initially
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
                contact_stress(i, 1) = contact_stress(i, 2); % Replace with the next value
            end
            if isnan(contact_stress(i, end)) || isinf(contact_stress(i, end))
                contact_stress(i, end) = contact_stress(i, end-1); % Replace with the previous value
            end
        end
    end
    nonan = contact_stress;
    return;
end



