from ..equal import AtomProperty, BondProperty

SubstructureAtomProperties = AtomProperty.AtomicNumber | AtomProperty.Aromatic

SubstructureBondProperties = BondProperty.Order | BondProperty.Aromatic

AutomorphAtomProperties = (
    AtomProperty.AtomicNumber
    | AtomProperty.Aromatic
    | AtomProperty.IsInRing
    | AtomProperty.HeavyDegree
)

AutomorphBondProperties = BondProperty.Aromatic
