from .node import DigraphNode


class PhantomNode(DigraphNode):
    """
    Node that is added to pad other node's such that they all have the same number of children.

    Section [P-92.1.4.1](https://iupac.qmul.ac.uk/BlueBook/P9.html#92010401) of the IUPAC Blue Book
    describes the use of phantom atoms.
    """

    @property
    def atomic_number(self) -> int:
        return 0

    @property
    def atomic_mass(self) -> float:
        return 0

    @property
    def name(self) -> str:
        return "0"
