CRAM Index Management
=====================

Verifies that ``.crai`` index files exist for all input CRAMs, creating them
with ``samtools index`` if missing. Required before any step that accesses
specific genomic regions.

.. automodule:: grid.utils.ensure_crai
   :members:
   :undoc-members:
   :show-inheritance:
