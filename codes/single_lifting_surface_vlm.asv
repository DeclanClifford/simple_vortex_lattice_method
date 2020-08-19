function[CL, CDi, delta]=single_lifting_surface_vlm(alpha, AR, N, sigma,...
                                                    sweep_angle,...
                                                    dihedral_angle)
% A simple vortex lattice method. Models a lifting surface as a flat plate
% finite wing, with control points placed to ensure agreement with slender
% wing theory at low aspect ratios and thin aerofoil theory at very large
% aspect ratios. The model can account for aspect ratio, taper ratio, 
% sweep angle and dihedral angle.

%%%NOTE%%%
% The model only works for even numbers of horseshoe vortex elements and 
% one horesehoe vortex elements. If you a enter a number of horseshoe
% vortex elements that is an odd number that is not one, it will add one
% element i.e. 51 elements becomes 52 elements.

% outputs
% CL = lift coefficient
% CDi = induced drag coefficient
% delta = induced drag factor

% inputs
% alpha = angle of attack
% AR = wing aspect ratio
% N = Number of horseshoe vortex elements
% sigma = wing taper ratio
% sweep_angle = wing sweep angle
% dihedral_angle = wing dihedral angle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vinf = 1;                           % free stream velocity
rho = 1.225;                        % fluid density
Q = [Vinf * cos(alpha);             % velocity (relative to zero lift line)
     0.0;
     Vinf * sin(alpha)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% wing geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this model sets the
% standard mean chord to 
% always == 1
croot = 2 / (1 + sigma);            % wing root chord
ctip = (2 * sigma) / (1 + sigma);   % wing tip chord
smc = 0.5 * (croot + ctip);         % wing standard mean chord
bw = AR * smc;                      % wing span


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% horseshoe vortex scaling factors %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enforcing even number of horseshoe vortex elements
if N > 1 && rem(N, 2) ~= 0
    N = N + 1;
end
    
large = 1.0e6;                      % length of the trailing vortices
fact = sqrt(N / (N + 1));           % horseshoe vortex size factor 
bp = bw * fact;                     % trailing vortex distance factor
dy = bp / N;                        % distance between trailing vortices

% NOTE - the control point distance factor, cp is defined later, after the
% taper ratio has been applied to each vortex element

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% vortex lattice method %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for zero aspect ratio wings i.e. no wing
if AR == 0
    
    % sets all outputs to zero
    CL = 0;
    CDi = 0;
    delta = 0;

% for non-zero aspect ratio wings    
else
    % locations of vortex segments, control and normal vectors   
    % one element used for wing
    if N == 1
        
        % control point distance factor for one horseshoe vortex
        cp = 0.5 * smc / fact;
        
        % coordinates of vortex segment on front wing
        xa(:, 1) = [0.0; - 0.5 * bp; 0.0];
        xb(:, 1) = [0.0; + 0.5 * bp; 0.0];
        xc(:, 1) = 0.5 * (xa(:, 1) + xb(:, 1)) + [cp; 0; 0];
        
        % normal vector for one vortex segment
        n(:, 1) = [0; 0; 1];
     
    % more than one element used for wing
    % NOTE - method only accepts even numbers of horseshoe vortex elements
    elseif N > 1 && rem(N, 2) == 0 
        
        % generates horseshoe vortex coordinates starting at the port wing
        % tip and ending at the starboard wing tip.
        
        % port wing
        xa(:, 1)=[(((N / 2) * - dy) + dy / 2) * tan(- sweep_angle);
                 (N / 2) * - dy;
                 ((N / 2) * dy) * sin(dihedral_angle)];
             
        xb(:, 1)=[(((N / 2) * - dy) + dy / 2) * tan(- sweep_angle);
                 ((N / 2) - 1) * - dy;
                 ((N / 2) - 1) * dy * sin(dihedral_angle)];

        for i = 2 : N / 2
            xa(:, i) = xb(:, i - 1) + [dy * tan(- sweep_angle); 0; 0];
            xb(:, i) = xa(:, i)+[0.0; dy; - dy * sin(dihedral_angle)];
        end
        
        % port wing normal vectors
        for i = 1 : N / 2
            n(:, i)=[0; sin(dihedral_angle); cos(dihedral_angle)];
        end
        
        % starboard wing
        xa(:, (N / 2) + 1) = [(dy / 2) * tan(sweep_angle); 0.0; 0.0];
        xb(:, (N / 2) + 1) = xa(:, (N / 2) + 1) +...
                             [0.0; dy ; dy * sin(dihedral_angle)];
        
        for i = (N / 2) + 2 : N
            xa(:, i) = xb(:, i - 1) + [dy * tan(sweep_angle); 0; 0];
            xb(:, i) = xa(:, i) + [0.0; dy; dy * sin(dihedral_angle)];
        end
        
        % starboard wing normal vectors
        for i = (N / 2) + 1 : N
            n(:, i) = [0; - sin(dihedral_angle); cos(dihedral_angle)];
        end       
               
    end 
    
    % establishing wing 
    for i = 1 : N
        % coords of centre of each horseshoe vortex
        centre_coords(i) = 0.5 * (xa(2, i) + xb(2, i));
        
        % assigning wing coords to horseshoe vortices using an equation
        % which gives the local chord for a straight tapered wing
        cw(i)=croot + (2 * abs(centre_coords(i)) / bw) * (ctip - croot);
        
        % control point distance factors
        cp(i) = 0.5 * cw(i) / fact;
        
        % control point locations
        xc(:, i) = 0.5 * (xa(:, i) + xb(:, i)) + [cp(i); 0; 0];
    end

    % set matrix and right hand side elements and solve for circulations
    for i = 1 : N
        for j = 1 : N
            I = vfil(xa(:, j),xb(:, j),xc(:, i));
            I = I + vfil(xa(:, j) + [large; 0; 0], xa(:, j), xc(:, i));
            I = I + vfil(xb(:, j), xb(:, j) + [large; 0; 0], xc(:, i));
            A(i, j) = I(1) * n(1, i) + I(2) * n(2, i) + I(3) * n(3, i);
        end
    rhs(i) = - (Q(1) * n(1, i) + Q(2) * n(2, i) + Q(3) * n(3, i));
    end
    gamma=A \ rhs';

    % forces at centres of bound vortices
    bc(:, :) = 0.5 * (xa(:, :) + xb(:, :));
    for i = 1 : N
    % first find the local velocity vector u on each load element i
    u = Q;

        for j = 1 : N
            u = u + vfil(xa(:, j), xb(:, j),bc(:, i)) * gamma(j);
            u = u + vfil(xa(:, j) +...
                    [large; 0; 0], xa(:, j), bc(:, i)) * gamma(j);
            u = u + vfil(xb(:, j), xb(:, j) +...
                    [large; 0; 0], bc(:, i)) * gamma(j);
        end
        
    % vector cross product of the local velocity and bound vortex direction
    % gives the direction of the force - multiplied by rho to give
    % the force in its correct units
    s = xb(:, i) - xa(:, i);
    Fx(i) = rho * (u(2) * s(3) - u(3) * s(2)) * gamma(i);
    Fz(i) = rho * (u(1) * s(2) - u(2) * s(1)) * gamma(i);
    end
    
    % lift coefficient
    CL = (sum(Fz(1 : N)) * cos(alpha) -...
          sum(Fx(1 : N)) * sin(alpha)) / (0.5 * rho * bw * smc);
    
    % induced drag coefficient
    CDi = (sum(Fx(1 : N)) * cos(alpha) +...
           sum(Fz(1 : N)) * sin(alpha)) / (0.5 * rho * bw * smc);
    
    % induced drag factor   
    delta = (CDi * pi * AR / CL^2) - 1;
end
