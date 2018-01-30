function [ggMatch,match] = matchPoints(x0,x1,knn)

% Create knn descriptor
knnfD = 16;
knnfM = 16;
[x0knnD,x0idx] = knnDescriptor(x0,knnfD,knnfM);
[x1knnD,x1idx] = knnDescriptor(x1,knnfD,knnfM);

% Match points based on descriptor
MdlKDT = KDTreeSearcher(x1knnD);
[match,D] = knnsearch(MdlKDT,x0knnD,'K',1);
match(D>=2) = 0;

% Find neighbors of unique match in original image idx
[u,ia,~] = unique(match);
n = histc(match,u);
n = n==1;
ia = nonzeros(ia.*n);
Tmatch = match(ia);
Fmatch = setdiff([1:length(x1)]',Tmatch);
ia = [ia;zeros(size(Fmatch))];
Tmatch = [Tmatch;Fmatch];
uMatch = [Tmatch,ia];
uMatch = sortrows(uMatch);
uMatch = uMatch(:,2);

sizex1idx = size(x1idx);
x10idx = uMatch(x1idx(:));
x10idx = reshape(x10idx,sizex1idx);

matchP = knnfM/2;;
% Compare condition of match
ggMatch = zeros(size(match));
for i = 1:length(x0)
    x0n = x0idx(i,:);
    x1n = x10idx(i,:);
    if length(intersect(x0n,x1n)) > matchP
        ggMatch(i) = match(i);
    end
end


end

