
%  Scalar model
%  Computing the focused field (eq 6): 
%    - Depends on:
%         - position,
%         - wavelength
%         - NA_c
%         - order-l spherical Bessel function of the first kind
%             - time-consuming
%             - depends on distance from p_f
%             - pre-compute and store in a N_lxR table => one sample per pixel.
%       
%         - order-l Legendre polynomial
%             - recursive definition
%             - depend on condenser angle alpha
%             - are constant for all points in E_f => precompute
%             - P_l(cos theta) depends on position 
%                 -compute independently for each point using a recursive
%                     definition for each order l=[0 N_l] 4 fpos
%     - Take advantage of the symmetry provided using a spherical condenser.
%     - Reconstruct the incident field between the condenser and the object
%       from a single cross-section of a cylindrical volume.
%         - Store the cross section as a RxR float array.
%         - Reference it in terms of cylindrical coordinates p \in [u,v] with
%             the origin located at p_f.
%         - Store only 1/4 of this cross-section by taking advantage of the
%             symmetry along the cylinder axis v.
%         - Correct for the phase shift by swapping the sign of the imaginary
%             component for v<0.
%         - Controll the accuaracy of the simulation by defining the desired
%             array resolution R. 
%     - Precompute a high-order representation of the entire field and store
%         efficiently  when the user makes a change to the condenser NA or incident wavelength.
%     - independent of the position and material of particles.
%     
% -Compute the internal and external scattered field
%     - cylindrically symmetric when the position of the sphere is at the focal point.
%         or when the incident light is represented by a single plane wave.
%     - utilize the symmetry in the plane-wave solution
%     - compute image and resample using Monte-Carlo integration 
%     - Components:
%         - the scattering coefficients given by B() and A() - eqs 8 and 11
%             -independent of position and therefore constant for any wavelength
%                 lambda and material m
%             - precompute for each scatterer
%         - propogation function given by h_l in E_s (dependent on sperhical Bessel
%             functions of the first and second kinds) and j_l in E_i.
%             -dependent on distance from sphere center
%             -1D and bounded by R
%             -precompute for a range of distance values
%             -the Hankel function in E_s is independent of any material properties
%                 and it is stored with j_l(ka) for the focused field E_f
%             -Hankel function = linear combination of Bessel functions of first
%                 and second kind
%                 -store j_l(ka) and y_l(ka) for both E_f and E_s
%             -Bessel function for E_i dependent on the index of refraction n(k)
%                     of the particle
%                 -store this table for each sphere
%                 -parameter range requires only d<r
%         -compute the Legendre polynomials recursively
%         - the Legendre polynomial P_l(cos theta)
%         - the phase shift exponentinal dependent on the light direction and 
%             sphere position
%         - how to compute the derivative of Bessel and Hankel functions depending
%             on boundary conditions - [13]
% -Monte-Carlo Integration
%     -determine the solutions for E_i and E_s using Monte-Carlo integration based
%         on a uniform distribution of sample points within the solid angle alpha
%         defined by NA_c
%     -determine the source points for plane waves using a stratified uniform
%         distribution projected onto the unit sphere based on Archimede's principle
%     -compute a uniform distribution of points in cylindrical coordinates in the
%         range phi=[0 2pi], z=[cos(alpha) 1]
% -Surface fields
%     -visualize the field using a geometric surface
% -Objective and detector
%     -the user specifies the field of view S and field resolution R
%     -perform fft of the field slice
%     -eliminate the field components above the cutoff frequency f_c
%     -user specified detector parameters
%     -specify material properties of the spheres at run-time or as a wavelength
%         dependent set of refractive indices.
%     -if the extinction spectrum is known use Kramers-Kronig relation to determine the phase speed
% -Accuracy and timing details
%     -time consuming operations:
%         -changes to the incident wavelength lambda and the condenser NA
%             -requires reevaluation of the focused field and all scattered fields
%     -allow the user to select the resolution of the simulation
%     -use stratified sampling to limit variance across the condenser aperture
%     -reconstruct randomly sample points when NA_c is changed
%     -use the same random seed
%    
%             

function runBimSim(varargin)




