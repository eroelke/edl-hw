function [c0, c1] = deadbConst(constr_i, constr_f)
%%%DEADBCONST Find deadband constants given initial and final constraints.
%
% Given the initial constraints constr_i and final constraints constr_f,
% the function computes the value of the constants c0, c1 specifying the
% velocity-dependent deadband.
% All units must be consistent, e.g. velocity in km/s, deadband values in
% deg, c0 in deg, c1 in deg/(km/s).
% 
% Inputs:
%    * constr_i, 2x1 double, initial velocity and deadband value.
%    * constr_f, 2x1 double, final velocity and deadband value.
% 
% Output:
%    * c0, double, deadband constant.
%    * c1, double, deadband constant.
% 
% Author:
%    Davide Amato, CU Boulder, davide.amato@colorado.edu.
% 
% History of changes:
%    20190702: Function created.

%% Unpack
v_i = constr_i(1); delta_i = constr_i(2);
v_f = constr_f(1); delta_f = constr_f(2);

%% c0, c1
c0 = (delta_i * v_f - delta_f * v_i) / ( v_f - v_i);
c1 = (delta_f - delta_i) / ( v_f - v_i);

end