function normalized = quantilenorm(X)
    % Simple quantile normalization for Octave
    % X: matrix where each column is a sample
    [n, m] = size(X);
    [X_sorted, idx] = sort(X);
    mean_sorted = mean(X_sorted, 2);
    normalized = zeros(n, m);
    for i = 1:m
        normalized(idx(:,i), i) = mean_sorted;
    end
end
