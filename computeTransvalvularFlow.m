function [mv,tv,av,pv] = computeTransvalvularFlow(flowDir, laDir)
    [mv,tv,av,pv] = RetrospectiveValveTracker([], flowDir, laDir, [], [], 0, 'valveLocations.mat');
end

