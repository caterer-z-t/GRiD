Python Pipeline Example
=======================

This example shows how to run the SOL pipeline directly in Python,
using functions imported from the GRiD CLI.  

This is useful for local debugging, Jupyter notebooks, or workflow integration.

.. literalinclude:: ../../../examples/SOL.py
   :language: python
   :linenos:
   :caption: "SOL.py - Python-based pipeline"

Usage
-----

.. code-block:: bash

   python SOL.py        # run all steps (1–9)
   python SOL.py 5      # run only step 5
   python SOL.py 3-7    # run steps 3 through 7
