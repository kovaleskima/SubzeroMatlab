function [p] = dsearchn_test(x,tri)
    xi = tri;
    t = zeros(size(xi,1),1);
    d = zeros(size(xi,1),1);
    for i = 1:size(xi,1) 
        yi = repmat(xi(i,:),size(x,1),1);
        [d(i),t(i)] = min(sum((x-yi).^2,2));
    end
    d = sqrt(d); 
    p = x(t(d<1),:);
end