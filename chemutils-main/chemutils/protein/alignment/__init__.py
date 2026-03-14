"""
Code for sequence alignment and mutation nomenclature.

This standardised nomenclature includes:
* Insertions, such as:
  * `A15_B16insGPP`, the insertion of GPP between an A at the 15th position and a B at the 16th position
  * `_A1insGPP`, the insertion of GPP before an A at the 1st position
* Deletions, such as:
    * `A15_B27del`, the deletion of the region between an A at the 15th position and a B at the 27th position
* Mutations, such as:
    * `G12D`, the replacement of the G at the 12th position by a D
"""

from .alignment import SequenceAlignment

__all__ = ["SequenceAlignment"]
