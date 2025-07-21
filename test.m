angles = -90:0.5:90;
P = sin(deg2rad(angles)).^2;
figure;
plot(angles, P, 'r'); grid on;
xlabel('Angle'); ylabel('Power');
title('Test Plot');