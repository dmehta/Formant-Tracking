function [] = plotCoast(E)
% If coasting was used, function indicates at what samples coasting was
% implemented
% Created : 3/21/2007, Dan S.
% Last Updated : 3/21/2007

for i = 1:size(E.formantInds,2)
    % Find non-zero indices of mask
    inds = find(E.formantInds(:,i) == 0);

    % Overlay data onto plots
    errorbar(inds,E.estTracks(i,inds,1)',sqrt(E.estVar(i,i,inds)),'k.');
end