Computing IBS/IBD Neighbors
===========================

Step 7 of GRiD (haplotype inference) requires a pre-computed haplotype neighbors
file as input. Two methods are supported, each producing a different file format
that GRiD's ``hi_inference`` module reads automatically:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Method
     - Tool
     - Best for
   * - **IBS**
     - `computeIBSpbwt <https://github.com/mhujoel/STRs/tree/main/cpp_files>`_ (Hujoel et al. 2026)
     - Cohorts without long-range IBD data; fast PBWT-based matching
   * - **IBD**
     - `iLASH <https://github.com/roohy/iLASH>`_
     - Larger cohorts where true IBD segments improve phasing specificity

Select the page for the method you want to use:

.. toctree::
   :maxdepth: 1

   ibs
   ibd