%-------------------------------------------------------------------------
%               Easy airfoil generator using Joukowsky transform
%
%                              Zexian Xiao
%-------------------------------------------------------------------------

% This project is inspired by Jens Nørkær Sørensen's course 46110 Fundamentals of Aerodynamics 
% Using the Joukowsky transformation to quickly generate airfoils

% Parameters for the airfoil.
% For example, NACA 4512 airfoil has a maximum camber of 4% located 50% (0.5 chords) from the leading edge with a maximum thickness of 12% of the chord.
Maximum_camber = 4;       % Target maximum camber (in percentage)
Relative_thickness = 12;  % Target relative thickness (in percentage)
Total_points_of_the_airfoil = 1001;  % Total points of the airfoil

% The Joukowsky transform cannot directly determine the location of the maximum camber.
% By scaling the x-coordinates of the airfoil, the location of the maximum camber can be aligned to the target location.
% If the airfoil is severely deformed, consider disabling it.
% Set to true to enable transformation, false to disable
enable_transformation = true;
maximum_camber_target_location = 50;  % Target location for maximum camber (in percentage)

s1 = Relative_thickness / 100;
s2 = Maximum_camber / 100;
maximum_camber_location = maximum_camber_target_location / 100;
%%
tolerance = 1e-6;
max_iterations = 1000;

for outer_iteration = 1:max_iterations
    for inner_iteration = 1:max_iterations
        c = 1;
        a = sqrt((s1 + c)^2 + s2^2);
        angle = linspace(0, 2 * pi, Total_points_of_the_airfoil);
        z_circle = a * (cos(angle) + 1i * sin(angle)) + (-1 * s1 + 1i * s2);
        z_ellipse = z_circle + c^2 ./ z_circle;

        % Normalize the airfoil coordinates
        Airfoil_x_coordinate = real(z_ellipse);
        Airfoil_y_coordinate = imag(z_ellipse);
        max_x_value = max(Airfoil_x_coordinate(:));
        min_x_value = min(Airfoil_x_coordinate(:));
        zoom = max_x_value - min_x_value;

        normalized_x = (Airfoil_x_coordinate - min_x_value) / zoom;
        normalized_y_1 = Airfoil_y_coordinate / zoom;
        [~, idx_x0] = min(abs(normalized_x));
        y_at_x0 = normalized_y_1(idx_x0);
        normalized_y = normalized_y_1 + (y_at_x0 * normalized_x - y_at_x0);

        % Store linear functions
        linear_segments = cell(length(normalized_x)-1, 1);
        for i = 1:length(normalized_x)-1
            x1 = normalized_x(i);
            y1 = normalized_y(i);
            x2 = normalized_x(i+1);
            y2 = normalized_y(i+1);
            m = (y2 - y1) / (x2 - x1);
            b = y1 - m * x1;
            linear_segments{i} = @(x) m * x + b;
        end

        % Lines for chord line divisions
        n = ceil(Total_points_of_the_airfoil / 2);
        camber_x_coordinate = linspace(min_x_value, max_x_value, n + 2);
        camber_x_coordinate = (camber_x_coordinate - min_x_value) / zoom;

        x_targets = camber_x_coordinate;
        y_targets = cell(length(x_targets), 1);
        camber_points = [];
        Actual_relative_thickness = 0;
        x_at_max_thickness = 0;

        % Iterate through all x_targets to find the corresponding y-coordinates and midpoints
        for j = 1:length(x_targets)
            x_target = x_targets(j);
            y_coords = [];

            for i = 1:length(normalized_x)-1
                x1 = normalized_x(i);
                x2 = normalized_x(i+1);

                if (x1 <= x_target && x_target <= x2) || (x2 <= x_target && x_target <= x1)
                    y_target = linear_segments{i}(x_target);
                    y_coords = [y_coords, y_target];
                end
            end

            y_targets{j} = y_coords;

            % Calculate the midpoint
            if length(y_coords) == 2
                midpoint_x = x_target;
                midpoint_y = mean(y_coords);
                camber_points = [camber_points; midpoint_x, midpoint_y];

                % Calculate the y difference and update max_y_diff and x_at_max_y_diff if necessary
                y_diff = abs(y_coords(1) - y_coords(2));
                if y_diff > Actual_relative_thickness
                    Actual_relative_thickness = y_diff;
                    x_at_max_thickness = x_target;
                end
            end
        end

        % Find the maximum value of y in the camber line and its corresponding x value
        [Actual_Maximum_camber, index] = max(camber_points(:, 2));
        x_at_max_camber = camber_points(index, 1);

        if abs(Actual_relative_thickness * 100 - Relative_thickness) < tolerance
            break;
        end

        s1 = s1 * Relative_thickness / (Actual_relative_thickness * 100);
    end

    if abs(Actual_Maximum_camber * 100 - Maximum_camber) < tolerance
        break;
    end
    
    s2 = s2 * Maximum_camber / (Actual_Maximum_camber * 100);
end

% Adjust airfoil x-coordinates if transformation is enabled
if enable_transformation
    x_at_max_camber = min(max(x_at_max_camber, 0), 1);
    transformed_x = normalized_x;
    transformed_y = normalized_y;
    transformed_x(normalized_x <= x_at_max_camber) = normalized_x(normalized_x <= x_at_max_camber) * (maximum_camber_location / x_at_max_camber);
    transformed_x(normalized_x > x_at_max_camber) = (normalized_x(normalized_x > x_at_max_camber) - x_at_max_camber) * ((1 - maximum_camber_location) / (1 - x_at_max_camber)) + maximum_camber_location;

    transformed_camber_x = camber_points(:, 1);
    transformed_camber_y = camber_points(:, 2);
    transformed_camber_x(camber_points(:, 1) <= x_at_max_camber) = camber_points(camber_points(:, 1) <= x_at_max_camber, 1) * (maximum_camber_location / x_at_max_camber);
    transformed_camber_x(camber_points(:, 1) > x_at_max_camber) = (camber_points(camber_points(:, 1) > x_at_max_camber, 1) - x_at_max_camber) * ((1 - maximum_camber_location) / (1 - x_at_max_camber)) + maximum_camber_location;
else
    transformed_x = normalized_x;
    transformed_y = normalized_y;
    transformed_camber_x = camber_points(:, 1);
    transformed_camber_y = camber_points(:, 2);
end

% Plot
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on;
plot([0, 1], [0, 0], '-r', 'LineWidth', 1);
plot(transformed_camber_x, transformed_camber_y, '-g', 'LineWidth', 1);
plot(transformed_x, transformed_y, '-b', 'LineWidth', 3, 'MarkerSize', 2);

if enable_transformation
    x_at_max_camber_transformed = x_at_max_camber;
    x_at_max_camber_transformed(x_at_max_camber <= x_at_max_camber) = x_at_max_camber * (maximum_camber_location / x_at_max_camber);
    x_at_max_camber_transformed(x_at_max_camber > x_at_max_camber) = (x_at_max_camber - x_at_max_camber) * ((1 - maximum_camber_location) / (1 - x_at_max_camber)) + maximum_camber_location;
    line([x_at_max_camber_transformed x_at_max_camber_transformed], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
    fprintf('Maximum camber of %.2f%% located at %.2f%% chords.\n', Maximum_camber, x_at_max_camber_transformed * 100);
    title_text = sprintf('Airfoil - %.0f %.0f %.0f', Maximum_camber, x_at_max_camber_transformed * 10, Relative_thickness);
    title(title_text);
else
    line([x_at_max_camber x_at_max_camber], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
    fprintf('Maximum camber of %.2f%% located at %.2f%% chords.\n', Maximum_camber, x_at_max_camber * 100);
    title_text = sprintf('Airfoil - %.0f %.0f %.0f', Maximum_camber, x_at_max_camber * 10, Relative_thickness);
    title(title_text);
end

line([x_at_max_thickness x_at_max_thickness], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1);
fprintf('Maximum thickness of %.2f%% located at %.2f%% chords.\n', Relative_thickness, x_at_max_thickness * 100);
xlabel('Normalized X');
ylabel('Normalized Y');
grid on;
axis equal;
xlim([0, 1]);
legend('Chord Line', 'Camber Line', 'Airfoil', 'Maximum Camber', 'Maximum Thickness', 'Location', 'southeast');
hold off;

%% Export
If the DXF file cannot be exported, move this code to the DXFLib_v0.9.1 folder and run it again.

filename = sprintf('%s.dxf', title_text);
FID = dxf_open(filename);

% Export chord line
dxf_polyline(FID, [0; 1], [0; 0], [0; 0]);

% Export camber line
for i = 1:numel(transformed_camber_x)-1
    dxf_polyline(FID, [transformed_camber_x(i); transformed_camber_x(i+1)], [transformed_camber_y(i); transformed_camber_y(i+1)], [0; 0]);
end

% Export airfoil
for i = 1:numel(transformed_x)-1
    dxf_polyline(FID, [transformed_x(i); transformed_x(i+1)], [transformed_y(i); transformed_y(i+1)], [0; 0]);
end

dxf_close(FID);
disp('DXF Created');
fprintf('File saved at: %s\n', fullfile(pwd, filename));