close all;
clear;
clc;

% The antenna's final version is a result of *4* iterations; the final
% iteration will undergo parametric analysis for Wslot-1 and Wslot-3 (*2*
% parameters) and those parameters will each be tested for *3* values
numberOfAntennas = 4;
numberOfParams = 2;
numberOfVariations = 3;

%% Dimensions
% Antenna
Wsub = 15e-3;
Lsub = 15e-3;
% Radiator patch
Wr = 7.5e-3;
Lr = 7.5e-3;
% Line feed
Wf = 2e-3;
Lf = 7e-3;
% First slot
Wslot_1 = [7.5e-3, 8.5e-3, 9.5e-3];
Lslot_1 = 2e-3;
% Second slot
Wslot_2 = 7e-3;
Lslot_2 = 1e-3;
% Third slot
Wslot_3 = [12.5e-3, 13.5e-3, 14.5e-3];
Lslot_3 = 1e-3;
% Spacers
L1 = 1.6e-3;
L2 = 1e-3;
% Ground and conductor-backed plane
Lgnd = 3e-3;
Lg = 1.8e-3;
Lremaining = Lsub - Lgnd - Lg - L1 - Lslot_1 - L2 - Lslot_2 - L2 - Lslot_3;
Ltop = 0.54*Lremaining;
% Feed
feedDiameter = Wf/50;
feedEntry = -Lsub/2;

%% Create Top Layer Shape
topPlane = antenna.Rectangle('Length', Wsub, 'Width', Lsub, 'NumPoints', 60);
feedRemoval1 = antenna.Rectangle('Length', Wsub/2 - Wf/2, 'Width', Lf, ...
    'Center', [(Wf/2 + (Wsub/2 - Wf/2)/2) (-Lsub/2 + Lf/2)]);
feedRemoval2 = antenna.Rectangle('Length', Wsub/2 - Wf/2, 'Width', Lf, ...
    'Center', [(-Wf/2 - (Wsub/2 - Wf/2)/2) (-Lsub/2 + Lf/2)]);
radRemoval1 = antenna.Rectangle('Length', Wsub/2 - Wr/2, 'Width', Lr, ...
    'Center', [(Wr/2 + (Wsub/2 - Wr/2)/2) (-Lsub/2 + Lf + Lr/2)]);
radRemoval2 = antenna.Rectangle('Length', Wsub/2 - Wr/2, 'Width', Lr, ...
    'Center', [(-Wr/2 - (Wsub/2 - Wr/2)/2) (-Lsub/2 + Lf + Lr/2)]);
topRemoval = antenna.Rectangle('Length', Wsub, 'Width', Lsub - Lf - Lr, ...
    'Center', [0 (Lsub/2 - (Lsub - Lf - Lr)/2)]);
topLayer = topPlane - feedRemoval1 - feedRemoval2 - radRemoval1 - radRemoval2 - topRemoval;

%% Create Bottom Layer Shapes
bottomPlane = antenna.Rectangle('Length', Wsub, 'Width', Lsub, 'NumPoints', 60);
% The final/proposed antenna will be analyzed for different slot sizes
thirdSlot = cell(1, numberOfVariations);
firstSlot = cell(1, numberOfVariations);
for n=1:1:numberOfVariations
    thirdSlot{n} = antenna.Rectangle('Length', Wslot_3(n), 'Width', Lslot_3, ...
        'Center', [(-Wsub/2 + Wslot_3(n)/2) (Lsub/2 - Lremaining - Lslot_3/2)]);
    firstSlot{n} = antenna.Rectangle('Length', Wslot_1(n), 'Width', Lslot_1, ...
        'Center', [(-Wsub/2 + Wslot_1(n)/2) (-Lsub/2 + Lgnd + Lg + L1 + Lslot_1/2)]);
end
topSlot = antenna.Rectangle('Length', Wsub, 'Width', Ltop, 'Center', ...
    [0 (Lsub/2 - Ltop/2)]);
secondSlot = antenna.Rectangle('Length', Wslot_2, 'Width', Lslot_2, ...
    'Center', [(Wsub/2 - Wslot_2/2) (Lsub/2 - Lremaining - Lslot_3 - L2 - Lslot_2/2)]);
gSlot = antenna.Rectangle('Length', Wsub, 'Width', Lg, ...
    'Center', [0 (-Lsub/2 + Lgnd + Lg/2)]);
% The bottom layers of all 4 versions of the antenna
bottomLayers = cell(1, numberOfAntennas);
bottomLayers{1} = bottomPlane - topSlot - gSlot;
bottomLayers{2} = bottomPlane - topSlot - firstSlot{numberOfVariations} - gSlot;
bottomLayers{3} = bottomPlane - topSlot - secondSlot - firstSlot{numberOfVariations} - gSlot;
bottomLayers{4} = bottomPlane - topSlot - thirdSlot{numberOfVariations} - secondSlot - firstSlot{numberOfVariations} - gSlot;

% All slot-size variants of the final bottom layer
bottom4Variations = cell(numberOfParams, numberOfVariations-1);
for m=1:1:numberOfParams
    for n=1:1:numberOfVariations-1
        if m==1
            bottom4Variations{m,n} = bottomPlane - topSlot - thirdSlot{numberOfVariations} - secondSlot - firstSlot{n} - gSlot;
        else
            bottom4Variations{m,n} = bottomPlane - topSlot - thirdSlot{n} - secondSlot - firstSlot{numberOfVariations} - gSlot;
        end
    end
end

%% Define Dielectric
d = dielectric('FR4');
d.EpsilonR = 4.4;
d.Thickness = 1.6e-3;
d.LossTangent = 0.024;

%% Combine the above to create and show the antennas
% All 4 versions of the antenna
antennas = cell(1, numberOfAntennas);
for n=1:1:numberOfAntennas
    antennas{n} = pcbStack;
    antennas{n}.Name = sprintf('Microstrip-fed slotted dual-band monopole antenna v.%d', n);
    antennas{n}.BoardShape = bottomPlane;
    antennas{n}.BoardThickness = d.Thickness;
    antennas{n}.Layers = {topLayer, d, bottomLayers{n}};
    antennas{n}.FeedLocations = [0, feedEntry, 1, 3];
    antennas{n}.FeedDiameter = feedDiameter;
    figure;
    show(antennas{n});
    % Create and show mesh
    figure;
    mesh(antennas{n}, 'MaxEdgeLength', .001, 'GrowthRate', 0.7);
end

% All slot-size variations of the proposed antenna
ant4Variations = cell(numberOfParams, numberOfVariations-1);
for m=1:1:numberOfParams
    for n=1:1:numberOfVariations-1
        ant4Variations{m,n} = pcbStack;
        ant4Variations{m,n}.Name = sprintf('Microstrip-fed slotted dual-band monopole antenna, scenario %d&%d', m, n);
        ant4Variations{m,n}.BoardShape = bottomPlane;
        ant4Variations{m,n}.BoardThickness = d.Thickness;
        ant4Variations{m,n}.Layers = {topLayer, d, bottom4Variations{m,n}};
        ant4Variations{m,n}.FeedLocations = [0, feedEntry, 1, 3];
        ant4Variations{m,n}.FeedDiameter = feedDiameter;
        figure;
        show(ant4Variations{m,n});
        % Create and show mesh
        figure;
        mesh(ant4Variations{m,n}, 'MaxEdgeLength', .001, 'GrowthRate', 0.7);
    end
end

%% Compute and plot S11s
fmin = 2e9;
fmax = 10e9;
Z0 = 50;
N = 500;
freq = linspace(fmin, fmax, N);
% All 4 versions of the antenna
lineStyles = {'--', '-.', ':', '-'};
lineColors = {'k', 'g', 'b', 'r'};
lineMarkers = {'none', '*', 'o', 'none'};
figure;
ax = gca;
for n=1:1:numberOfAntennas
    s = sparameters(antennas{n}, freq, Z0);
    lnBefore = findobj(ax, 'Type', 'line');
    rfplot(s, 1, 1);
    hold on;
    lnAfter = findobj(ax, 'Type', 'line');
    newLn = setdiff(lnAfter, lnBefore);
    set(newLn, 'Color', lineColors{n}, 'LineStyle', lineStyles{n}, 'DisplayName', sprintf('Ant. %d', n), 'Marker', lineMarkers{n});
end
hold off;

% All slot-size variations of the proposed antenna
labelNums = [1 3];
labelValues = [Wslot_1*10^3; Wslot_3*10^3];
lineStyles = {'--', '-.', '-'};
lineColors = {'k', 'g', 'r'};
lineMarkers = {'none', '*', 'none'};
for m=1:1:numberOfParams
    figure;
    ax = gca;
    for n=1:1:numberOfVariations-1
        sp = sparameters(ant4Variations{m,n}, freq, Z0);
        lnBefore = findobj(ax, 'Type', 'line');
        rfplot(sp, 1, 1);
        hold on;
        lnAfter = findobj(ax, 'Type', 'line');
        newLn = setdiff(lnAfter, lnBefore);
        set(newLn, 'Color', lineColors{n}, 'LineStyle', lineStyles{n}, 'DisplayName', sprintf('W slot-%d = %g mm', labelNums(m), labelValues(m,n)), 'Marker', lineMarkers{n});
    end
    lnBefore = findobj(ax, 'Type', 'line');
    rfplot(s, 1, 1);
    lnAfter = findobj(ax, 'Type', 'line');
    newLn = setdiff(lnAfter, lnBefore);
    set(newLn, 'Color', lineColors{n+1}, 'LineStyle', lineStyles{n+1}, 'DisplayName', sprintf('W slot-%d = %g mm', labelNums(m), labelValues(m,n+1)), 'Marker', lineMarkers{n+1});
    hold off;
end

% Compute and plot current distribution
freqs = [2.4e9, 4e9, 5e9, 7.14e9, 9.05e9];
for n=1:1:length(freqs)
    figure;
    current(antennas{numberOfAntennas}, freqs(n), "metal", Direction="on", Scale="log");
end

%% Compute and plot radiation patterns
freqs = [2.4e9, 4e9; 5.2e9, 5.8e9; 7.14e9, 9.05e9];
planes = ["E-Plane (y-z)", "H-Plane (x-z)"];
pol = ["V", "H"; "H", "V"];
azimuthAngle = [90, 0];
for m=1:1:size(freqs, 1)
    for n=1:1:length(planes)
        figure;
        pattern(antennas{numberOfAntennas}, freqs(m,:), azimuthAngle(n), 0:1:360, Polarization=pol(n, 1), PlotStyle="overlay");
        currentPattern = polarpattern("gco");
        currentPattern.TitleTop = sprintf("%s, Co-polarization", planes(n));
        currentPattern.LineWidth = 3;

        figure;
        pattern(antennas{numberOfAntennas}, freqs(m,:), azimuthAngle(n), 0:1:360, Polarization=pol(n, 2), PlotStyle="overlay");
        currentPattern = polarpattern("gco");
        currentPattern.TitleTop = sprintf("%s, Cross-polarization", planes(n));
        currentPattern.LineStyle = '--';
        currentPattern.LineWidth = 3;
    end
end
