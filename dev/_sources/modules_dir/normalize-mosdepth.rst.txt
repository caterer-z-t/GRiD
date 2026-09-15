Coverage Normalization
======================

Normalizes mosdepth coverage in two steps: within-individual (divide by sample mean)
and across-individual (z-score transformation). High-variance regions are then selected
to reduce noise before neighbor computation.

.. automodule:: grid.utils.normalize_mosdepth
   :members:
   :undoc-members:
   :show-inheritance:
