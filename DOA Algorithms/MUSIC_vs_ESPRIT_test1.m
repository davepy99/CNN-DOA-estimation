% === Plot MUSIC Spectrum with ESPRIT DOAs ===
figure;
plot(angles_scan, P_music_dB, 'b', 'LineWidth', 2); hold on;
ylim([-50 0]); grid on;

% Overlay ESPRIT DOAs as vertical lines
for i = 1:length(doa_esprit)
    xline(real(doa_esprit(i)), 'r--', 'LineWidth', 2, ...
        'DisplayName', sprintf('ESPRIT DOA: %.2f°', real(doa_esprit(i))));
end

% Overlay True DOAs
for i = 1:length(trueAngles)
    xline(trueAngles(i), 'k:', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('True DOA: %d°', trueAngles(i)));
end

xlabel('Angle (degrees)');
ylabel('Spatial Spectrum (dB)');
title('MUSIC Spectrum with ESPRIT DOA Estimates');
legend('MUSIC Spectrum', 'Location', 'Best');
xlim([-90 90]);
