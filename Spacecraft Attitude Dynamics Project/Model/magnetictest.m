clear
clc

[g, h] = IGRF13;
r=500+6371;

for i=1:999
    theta=2*pi*i/1000;
    for j=1:999
        phi=-pi/2+pi*j/1000;
        r_N=r.*[cos(theta)*cos(phi);sin(theta)*cos(phi);sin(phi)];
        [B_N]  = earthmagfield13(r_N, 3600*9, g, h, 0, 13);
        B(j,i)=norm(B_N);
    end
end
X=0:2*pi/1000:2*pi-4*pi/1000;
Y=-pi/2:2*pi/2000:pi/2-4*pi/2000;
%%



% Load the background image of Earth
earth = imread('earth.jpeg'); % Ensure the image is in the current directory or provide the full path

% Display the background image
figure;
imagesc([0, 360], [-90, 90], flipud(earth)); % Adjust to match degrees
set(gca, 'YDir', 'normal'); % Correct orientation
hold on;

% Convert X and Y to degrees for the contour plot
contour(X * 180/pi, Y * 180/pi, B, 30, 'LineColor', 'w', 'LineWidth', 1); % White contour lines for contrast

% Customize axes
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
%title('Contour Plot with Earth Background');
axis tight equal; % Ensure axes match
hold off;
