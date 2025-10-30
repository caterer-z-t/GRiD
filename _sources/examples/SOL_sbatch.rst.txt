SLURM Job Example
=================

This example shows how to run the SOL pipeline using an HPC SLURM submission script.

The job supports flexible step selection and can run individual pipeline steps or ranges.

.. literalinclude:: ../../../examples/SOL.sh
   :language: bash
   :linenos:
   :caption: "SOL.sh - Example SLURM submission script"

Usage
-----

.. code-block:: bash

   sbatch SOL.sh        # run all steps (1–9)
   sbatch SOL.sh 5      # run only step 5
   sbatch SOL.sh 3-7    # run steps 3 through 7
