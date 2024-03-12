function st = localSTfun(ants,w)
% Calculates for each focal point the local straightness as the distance
% between the 2 points ±w away from the focal point

if ~mod(w,2); w = w+1; end % If w is even it's made odd 2B divisible by 2
pad = nan(floor(w/2),1); % NaN padding start & end

% % Make table of x,y of pre point, focal point, post point
% ids = unique(ants.id);
% stC = cell(size(ids));
% for id = 1:numel(ids)
%     ant = ants(ants.id==ids(id),:);
% %     ant(isnan(ant.alpha),:) = [];
%     d = sqrt(sum([ant.x(1:end-w+1) - ant.x(w:end), ant.y(1:end-w+1) - ant.y(w:end)].^2, 2));
%     cs = cumsum(ant.s,'omitnan');
%     l = cs(w:end) - cs(1:end-w+1);
%     stC{id,1} = [pad; d./l; pad];
% end
% st = single(vertcat(stC{:}));

% Make table of x,y of pre point, focal point, post point

d = sqrt(sum([ants.x(1:end-w+1) - ants.x(w:end), ants.y(1:end-w+1) - ants.y(w:end)].^2, 2));
cs = cumsum(ants.s,'omitnan');
l = cs(w:end) - cs(1:end-w+1);
st = single([pad; d./l; pad]);
gapInd = [find(isnan(ants.s)); find([false; diff(ants.id)>0])];
gapIndC = cell(numel(gapInd)-2,1);
for g = 2:numel(gapInd)-1
    gapIndC{g,1} = (gapInd(g)-floor(w/2):gapInd(g)+floor(w/2))';
end
gapInds = vertcat(gapIndC{:});
gapInds(gapInds<1 | gapInds > numel(st)) = [];
st(gapInds) = nan;
