function x = znormalize(x)
    x = (x - mean(x,2))./std(x,0,2);
end

