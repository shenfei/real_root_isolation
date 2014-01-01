function [rootCell,nroot] = HOM4PS2(n,P,mthd);

%
%      Homotopy method for solving polynomial systems.
%
%  Syntax:   >> [rootCell,nroot] = HOM4PS2(n,P,mthd)
%  Input :
%           n -- the number of equations.
%           P -- polynomial system (string form).
%        mthd -- (optional)
%                = 1 : use the polyhedral homotopy method. (default)
%                = 2 : use the classical homotopy method.
%
%  output :
%    rootCell -- Matlab Cell containing information of solutions.
%             rootCell{:,1} : solutions.
%             rootCell{:,2} : the residue of solutions.
%             rootCell{:,3} : the condition number of solutions.
%
%        e.g. rootCell{1,1} is the first solution.
%             rootCell{1,2} is the residue of the first solution.
%             rootCell{1,3} is the condition number of the first solution.
%
%       nroot -- the number of solutions.
%
%  Example:
% >> P = '{3+x^3*y^2+x; 4*y*z^5+8*x^2*y^4*z^4-1; x+y+z-1;}'
% >> [rootCell,nroot] = HOM4PS2(3,P);
%

if (nargin < 2)
   error('Two inputs required.')
end

if (nargin < 3)
    mthd = 1;
end

if ~(mthd == 1 | mthd == 2)
    error('(The third input (mthd) should be 1 or 2.')
end

%----------------------------------------------

   fid = fopen('PS.sym','w');
   fprintf(fid,'%s',P);
   fclose(fid);

   !sym2num < PS.sym > PS.num

   [nroot,ncrv,sol_info,roots] = HOM4PS2_ker(mthd,rand(1,n));

   for ind = 1: nroot
       rootCell{ind,1} = roots(:,ind);
       rootCell{ind,2} = sol_info(1,ind);
       rootCell{ind,3} = sol_info(2,ind);
   end

   if (nroot >= 500000)
       disp('Only the first 500,000 roots will be output.')
   end

    fprintf(' The # of roots:            %g\n',nroot);
    fprintf(' The # of paths followed:   %g\n',ncrv);

end
