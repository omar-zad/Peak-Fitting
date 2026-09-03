# Supervisor-facing methods summary

One-dimensional GIWAXS FR, IP and OOP line profiles were fitted without
smoothing or vertical intensity offsets. Candidate peaks were represented by
area-normalised Gaussian, Lorentzian or pseudo-Voigt components plus a bounded
local polynomial background. Peak centres and widths were constrained to
predefined physically and visually justified intervals. Fits used bounded
nonlinear least squares with a small multistart set; no MCMC or unconstrained
exhaustive search was used.

The profile family was selected separately for each q-window by comparing BIC
across all active frames in that window. A simpler Gaussian or Lorentzian model
was retained when pseudo-Voigt improved BIC by less than two units. Component
evidence was evaluated by removing each peak and refitting, combining the
resulting delta BIC with integrated-area signal-to-noise and parameter-bound
checks.

For the Drop40 low-q window, component number was examined independently of
profile family by fitting matched one-, two- and three-component candidates to
the same data points and background form. Global BIC supported three components
for FR, two for IP and one for OOP. This model-level comparison is distinct
from the leave-one-component-out test used to classify individual peaks within
each frame.

Every accepted result required numerical convergence and visual agreement of
the raw data, total model, individual components and residuals. Reported q
positions and d-spacings were separated from internal numerical fit values so
that tentative, bound-limited and reference-stabilised components could remain
in the audit trail without being presented as independent measurements.

Apparent FWHM values are conditional on the selected line-shape/background
model and are not instrument-deconvolved. They were therefore treated as fit
diagnostics for supporting samples rather than converted directly into CCL,
crystallite size or strain.
