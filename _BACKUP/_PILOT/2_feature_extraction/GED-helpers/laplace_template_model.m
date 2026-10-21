function yhat = laplace_template_model(x, p)
% Normalized Laplace peak template used for PF scoring.
A = p(1);
mu = p(2);
b = p(3);
yhat = A .* exp(-abs(x - mu) ./ max(b, eps));
yhat = max(yhat, 0);
if max(yhat) > 0
    yhat = yhat ./ max(yhat);
end
end
