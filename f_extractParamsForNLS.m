%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (C) 2016 Aaron Staranowicz and Gian Luca Mariottini
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%
%%

if f_Opt.method == 1    % 6pnt
    disp('6pnt')
    switch f_Opt.case
        case 1% K_D, K_R only
            disp('K_D, K_R only')
            costFunc = @(X)f_min6pntKinect(X,Sol_LS, f_Opt);
            X0 = [D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0];
        case 2 % R and t only
            disp(' R and t only')
            costFunc = @(X)f_min6pntKinect(X,Sol_LS, f_Opt);
            X0 = [r0,p0,y0, tx,ty,tz];
        case 3 % lens distortion only
            disp('lens distortion only')
            costFunc = @(X)f_min6pntKinect(X,Sol_LS, f_Opt);
            X0 = [k_c0, k_d0];
        case 4 % all but lens
            disp('all but lens')
            costFunc = @(X)f_min6pntKinect(X,Sol_LS, f_Opt);
            X0 = [r0,p0,y0, tx,ty,tz, D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0];
        case 5 % all include lens
            disp('all include lens')
            costFunc = @(X)f_min6pntKinect(X,Sol_LS, f_Opt);
            X0 = [r0,p0,y0, tx,ty,tz, D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0, k_c0, k_d0];
    end
else    % Conic
    disp('Conic')
    switch f_Opt.case
        case 1% K_D, K_R only
            disp('K_D, K_R only')
            costFunc = @(X)f_minConicKinect(X,Sol_LS, f_Opt);
            X0 = [D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0];
        case 2 % R and t only
            disp(' R and t only')
            costFunc = @(X)f_minConicKinect(X,Sol_LS, f_Opt);
            X0 = [r0,p0,y0, tx,ty,tz];
        case 3 % lens distortion only
            disp('lens distortion only')
            costFunc = @(X)f_minConicKinect(X,Sol_LS, f_Opt);
            X0 = [k_c0, k_d0];
        case 4 % all but lens
            disp('all but lens')
            costFunc = @(X)f_minConicKinect(X,Sol_LS, f_Opt);
            X0 = [r0,p0,y0, tx,ty,tz, D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0];
        case 5 % all include lens
            disp('all include lens')
            costFunc = @(X)f_minConicKinect(X,Sol_LS, f_Opt);
            X0 = [r0,p0,y0, tx,ty,tz, D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0, k_c0, k_d0];
    end
end