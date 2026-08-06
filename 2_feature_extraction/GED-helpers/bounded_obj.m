function v = bounded_obj(p, lb, ub, obj_fun)
% Evaluate obj_fun after clipping parameters to [lb, ub].
p_clip = min(max(p, lb), ub);
v = obj_fun(p_clip);
if ~isfinite(v)
    v = 1e6;
end
end
