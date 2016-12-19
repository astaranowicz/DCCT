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

switch f_Opt.case
    case 1% K_D, K_R only
        %X0 = [D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0];
        K_R = [X(6) 0*X(8) X(9);
                0    X(7) X(10);
                0     0     1];
        K_D = [X(1) 0*X(3) X(4);
                0    X(2) X(5);
                0    0     1];
        R_R_D = Sol_LS.R_R_D;
        R_t_D = Sol_LS.R_t_D;
        k_c = Sol_LS.k_c;
        k_d = Sol_LS.k_d;
    case 2 % R and t only
        %X0 = [r0,p0,y0, tx,ty,tz];
        K_R = Sol_LS.K_R;
        K_D = Sol_LS.K_D;
        R_R_D = f_rpy2R([X(1), X(2), X(3)]);
        R_t_D = [X(4); X(5); X(6)];
        k_c = Sol_LS.k_c;
        k_d = Sol_LS.k_d;
    case 3 % lens distortion only
        %X0 = [k_c0, k_d0];
        K_R = Sol_LS.K_R;
        K_D = Sol_LS.K_D;
        R_R_D = Sol_LS.R_R_D;
        R_t_D = Sol_LS.R_t_D;
        k_c = [X(1); X(2); X(3); X(4);X(5)];
        k_d = [X(6); X(7); X(8); X(8); X(10)];
    case 4 % all but lens
        %X0 = [r0,p0,y0, tx,ty,tz, D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0];
        K_R = [X(12) 0*X(14) X(15);
                0     X(13) X(16);
                0     0     1];
        K_D = [X(7) 0*X(9) X(10);
                0    X(8) X(11);
                0    0     1];
        R_R_D = f_rpy2R([X(1), X(2), X(3)]);
        R_t_D = [X(4); X(5); X(6)];
        k_c = Sol_LS.k_c;
        k_d = Sol_LS.k_d;
    case 5 % all include lens
        %X0 = [r0,p0,y0, tx,ty,tz, D_fu,D_fv,D_s,D_u0,D_v0, R_fu,R_fv,R_s,R_u0,R_v0, k_c0, k_d0];
        K_R = [X(12) 0*X(14) X(15);
                0     X(13) X(16);
                0     0     1];
        K_D = [X(7) 0*X(9) X(10);
                0    X(8) X(11);
                0    0     1];
        R_R_D = f_rpy2R([X(1), X(2), X(3)]);
        R_t_D = [X(4); X(5); X(6)];
        k_c = [X(17); X(18); X(19); X(20);X(21)];
        k_d = [X(22); X(23); X(24); X(25); X(26)];
end
