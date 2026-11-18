function rate = fig3_kill(CB, kind, p)
switch kind
    case 'linear'
        rate = p(1)*CB;
    case 'hill1'
        rate = p(1) ./ (1+p(2)./CB);
    case 'hill'
        rate = p(1) ./ (1+(p(2)./CB).^p(3));
end
end
