% Fitting ellipses, ellipsoids and other quadratic curves and surfaces
% Copyright 2011-2013 Levente Hunyadi
%
% Quadratic fitting functions relate to Thesis 1 in the PhD dissertation "Estimation methods in the
% errors-in-variables context" by Levente Hunyadi, submitted in 2013 to the Budapest University of
% Technology and Economics.
%
% Distance of a point and a quadratic curve
%   quad2dproj                - Projects a set of points onto a quadratic curve.
%
% Generating points along an ellipse
%   ellipse                   - Generates points along an ellipse.
%   ellipse_orbital           - Generates points along an ellipse using Kepler's parameters.
%
% Different representation of an ellipse
%   is_ellipse                - Test if implicit parameters represent an ellipse.
%   ellipse_ex2im             - Cast ellipse defined with explicit parameters to implicit form.
%   ellipse_im2ex             - Cast ellipse defined with implicit parameter vector to explicit form.
%   ellipse_im2kepler         - Cast ellipse defined with standard parameter vector to Kepler form.
%   ellipse_kepler2im         - Cast ellipse defined in Kepler form to standard parameter vector form.
%
% Distance of a point and an ellipse
%   ellipse_distance          - Distance of points projected onto an ellipse.
%
% Fitting ellipses to data
%   ellipsefit                - Fit an ellipse to data by minimizing point-to-curve distance.
%   ellipsefit_direct         - Direct least squares fitting of ellipses.
%   ellipsefit_foot           - Maximum likelihood ellipse fit with foot points.
%   ellipsefit_koopmans       - Fit an ellipse to data using the nonlinear Koopmans method.
%   ellipsefit_robust         - Constrained ellipse fit by solving a modified eigenvalue problem.
%
% Fitting circles to data
%   circlefit                 - Fit a circle to data using the least squares approach.
%
% Fitting parabolas to data
%   parabolafit_cals          - Fit a parabola using consistent algebraic least squares.
%   parabolafit_direct        - Direct least squares fitting of parabolas.
%   parabolafit_directm       - Direct least squares fitting of parabolas using a pre-processed scatter matrix.
%
% Different representations of an ellipsoid
%   is_ellipsoid              - Test if implicit parameters represent an ellipsoid.
%   ellipsoid_ex2im           - Cast ellipsoid defined with explicit parameters to implicit vector form.
%   ellipsoid_im2ex           - Cast ellipsoid defined with implicit parameter vector to explicit form.
%
% Distance of a point and an ellipsoid
%   ellipsoid_distance        - Distance of points projected onto an ellipsoid.
%
% Fitting ellipsoids to data
%   ellipsoidfit              - Fit an ellipsoid to data by minimizing point-to-surface distance.
%   ellipsoidfit_leastsquares - Fit an ellipsoid to a set of 3D data points.
%   ellipsoidfit_aml          - Approximated maximum likelihood fit of ellipsoids.
%   ellipsoidfit_direct       - Direct least squares fitting of ellipsoids under the constraint 4J - I^2 > 0.
%   ellipsoidfit_simple       - Simple iterative least squares fitting of ellipsoids under the constraint k*J - I^2 > 0.
%   ellipsoidfit_koopmans     - Fit an ellipsoid to data using the nonlinear Koopmans method.
%   ellipsoidfit_iterative    - Iterative least squares fitting of ellipsoids under the constraint k*J - I^2 > 0.
%   ellipsoidfit_residuals    - Finds the distance of the point (x,y,z) to the nearest point on the ellipsoid.
%
% Plotting ellipsoids
%   plot_ellipsoid            - Plot ellipsoid specified with center, radii and rotation matrix.
%   plot_ellipsoid_im         - Plot ellipsoid specified with implicit parameters.
%   plot_ellipsoid_part       - Plot ellipsoid specified with center, radii and rotation matrix.
%   ellipsoid_projections     - Plots ellipsoid projections on the xy, xz and yz planes.
%
% Generating points along a sphere
%   sphere_gd                 - Construct a geodesic sphere.
%
% Fitting a sphere to data
%   spherefit                 - Fit a sphere to data using the least squares approach.
%
% Plotting spheres
%   plot_sphere_part          - Partial sphere surface.
%
% Fitting quadratic curves to data
%   quad2dfit_cals            - Fit a quadratic curve to data using consistent algebraic least squares.
%   quad2dfit_koopmans        - Fit a quadratic curve to data using the nonlinear Koopmans method.
%   quad2dfit_leastsquares    - Fit a quadratic curve to a set of 2D data points using least-squares.
%   quad2dfit_lsnormal        - Fit a quadratic curve to a set of 2D data points with normals.
%   quad2dfit_taubin          - General quadratic curve fit with Taubin's method.
%   quad2dfit_hyperaccurate   - General quadratic curve fit with Kanatani's hyperaccurate fit method.
%
% Fitting quadratic curves to data with constraints
%   quad2dconfit_koopmans     - Fit a constrained quadratic curve to data using the nonlinear Koopmans method.
%
% Fitting quadratic surfaces to data
%   quad3dfit_koopmans        - Fit a quadratic curve to data using the nonlinear Koopmans method.
%   quad3dfit_taubin          - General quadric surface fit with Taubin's method.
%
% Projection of a point to a quadratic curve or surface
%   ellipseproj               - Projects a set of points onto an ellipse.
%   quadprojfun               - The projection function Q(t) for various conics and quadrics.
%
% Covariance matrices for quadratic curve and surface fitting
%   quadmake                  - Builds covariance matrices for quadratic forms.
%   quad2dcov                 - Theoretical noise covariance structure corresponding to a 2D quadratic function.
%   quad2dcovred              - Reduced-size noise covariance structure corresponding to a 2D quadratic function.
%   quad3dcovmat              - Theoretical noise covariance matrix corresponding to a 3D quadratic function.
%   quad3dcovpoly             - Theoretical noise covariance structure corresponding to a 3D quadratic function.
%
% Fitting planes to data
%   planefit                  - Plane fit to noisy 3D data points.
%
% Plotting planes
%   plotconvhull              - Plot convex hull of a set of data points.
%
% Miscellaneous plotting functions
%   lineclip                  - Draw a line with clipping using the Cohen-Sutherland algorithm.
%   nntrifaces                - Nearest-neighbor triangle face vertex indices.
%
% Similarity transformations
%   quad2d_similarity         - Translate data to place centroid at origin and apply isotropic scaling.
%
% Rotations
%   rot3d                     - Rotate points in three dimensions.
%   ang2rot                   - Convert Euler angles to rotation matrix.
%   quat2rot                  - Convert quaternion to rotation matrix.
%   rot2ang                   - Convert rotation matrix to Euler angles.
%   rot2quat                  - Convert rotation matrix to quaternion.
%
% Translations
%   quad3d_center             - Center of a central quadric (quadratic surface).
%   quad3d_translate          - Translates a quadratic surface in implicit form by the given coordinates.
%   quad2d_translate          - Translates a quadratic curve in implicit form by the given coordinates.
%
% Plot functions
%   plot_results              - Plot several data series against a single data series.
%
% Utility functions
%   quad2dfit_paramchk        - Quadratic 2-D fitting parameter check.
%   quad3dfit_paramchk        - Quadratic 2-D fitting parameter check.
%
% 2D examples
%   example_ellipsefit1       - Demonstration of various ellipse (and general quadratic curve) fits.
%   example_ellipsefit2       - Various ellipse fits to some data points.
%   example_ellipsefitnoise   - Comparing ellipse fits under various noise conditions.
%   example_ellipsefitsector  - Ellipse fitting to data from a limited sector.
%   example_ellipsenormfit    - Comparing the accuracy of various ellipse fits with and without normals.
%   example_ellipseproj       - Demonstration of projecting points to an ellipse.
%   example_parabolafit       - Demonstration of parabola fitting.
%   example_quad2dcon         - Demonstration of constrained quadratic curve fitting.
%   example_quad2dcomp        - Comparison of various ellipse (and general quadratic curve) fits.
%   example_kcrlb             - Demonstration of fitting algoritm performance compared to KCR-LB.
%
% 3D examples
%   example_ellipsoid         - A sample ellipsoid.
%   example_ellipsoidfit1     - Demonstration of various ellipsoid fits.
%   example_ellipsoidfit2     - Demonstration of fits to various special ellipsoids.
%   example_ellipsoidcomp     - Comparative demonstration of fits to random ellipsoids.
%   example_spherefit         - Demonstration of least-squares sphere fit.
%   example_planefit          - Demonstration of plane fitting.
%
% Eigenvalues and eigenvectors
%   eigsm                     - Smallest eigenvalue with real eigenvector.
%   gsvd_min                  - Singular vector belonging to smallest singular value.
%   mpolyeig                  - Eigenvalues and eigenvectors for matrix polynomial.
%   mpolyval                  - Evaluates a matrix polynomial T(s) at the specified value of s.
%   mpolycomp                 - Matrix polynomial companion form.
%   mpolycomps                - Transforms a matrix polynomial into a (symmetric) linearized form.
%
% Symbolic manipulation
%   ellipsoidsym              - Symbolic expression of ellipsoid in implicit form with center, axes and quaternions.
%   symcrosscov               - Compensation for observation cross-covariance matrix.
%
%
% Quadratic curves and quadric surfaces in implicit form
% Copyright 2010 Levente Hunyadi
%
% Examples
%   example_conic      - Demonstrates plotting conic sections given as an implicit equation.
%   example_pqdist     - Sample code for distance of quadratic curve and point.
%
% Conic sections as an implicit equation
%   imconic            - Compute parameters of or plot conic section given in implicit form.
%   imconicdiscr       - Discriminant for conic section.
%   imconicisect       - Intersections of conic with window bounds.
%   imconicrotate      - Rotates a conic section given as an implicit equation.
%   imconicrotation    - Rotation matrix of a conic section given as an implicit equation.
%   imconictranslate   - Translates a conic section given as an implicit equation.
%   imconictranslation - Translation vector of a conic section given as an implicit equation.
%
% Distance of point and quadratic curve or quadric surface
%   pqdist             - Distance of point(s) and quadratic curve or quadric surface.
%   pqdistpoly         - Polynomial in parameter t for distance between point y and foot point x.
%
% Quadratic surface classification
%   imquad             - Identify type of quadratic surface.
%
% Geometric transformations
%   invtransrot        - Apply translation and inverse rotation to a set of points.
%   transrot           - Apply translation and rotation to a set of points.
%
% Matrices in symbolic variables
%   example_symm       - Sample code for matrices in symbolic variables
%   symm               - Create a matrix of symbolic variables.
%   symmd              - Create a diagonal matrix of symbolic variables.
%   symms              - Create a symmetric matrix of symbolic variables.
%   sym2matlab         - Convert symbolic expression containing matrix entries into MatLab.
%   symmvars           - Parse parameters for symbolic matrix creation.
%
% Utility functions
%   intfilter          - Filters those elements of x that are outside a set of intervals.
%
%
% Polynomials in symbolic variables
% Copyright 2011 Levente Hunyadi
%
% Symbolic polynomials
%   sympoly      - Symbolic polynomial array of multiple variables.
%   sympolys     - Create multiple symbolic polynomials in the caller workspace.
%
% Symbolic variables
%   symvariable  - Symbolic variable in a symbolic expression.
%   symvard      - Symbolic variable in a symbolic expression with time dynamics.
%
% Unit tests
%   test_sympoly - Test suite for symbolic polynomial manipulation.
%
% Utility functions
%   latexgreek   - Replace all recognized Greek letters with their LaTeX equivalents.
%   num2string   - Convert a numeric value to a string in the specified format.
%   strjoin      - Concatenate a cell array of strings.
