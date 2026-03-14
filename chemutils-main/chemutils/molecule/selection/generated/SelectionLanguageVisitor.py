# Generated from SelectionLanguage.g4 by ANTLR 4.13.2
from antlr4 import *
if "." in __name__:
    from .SelectionLanguageParser import SelectionLanguageParser
else:
    from SelectionLanguageParser import SelectionLanguageParser

# This class defines a complete generic visitor for a parse tree produced by SelectionLanguageParser.

class SelectionLanguageVisitor(ParseTreeVisitor):

    # Visit a parse tree produced by SelectionLanguageParser#full_expr.
    def visitFull_expr(self, ctx:SelectionLanguageParser.Full_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#parenthesized_expr.
    def visitParenthesized_expr(self, ctx:SelectionLanguageParser.Parenthesized_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#not.
    def visitNot(self, ctx:SelectionLanguageParser.NotContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#same.
    def visitSame(self, ctx:SelectionLanguageParser.SameContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#or.
    def visitOr(self, ctx:SelectionLanguageParser.OrContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#within.
    def visitWithin(self, ctx:SelectionLanguageParser.WithinContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#and.
    def visitAnd(self, ctx:SelectionLanguageParser.AndContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#parenthesized.
    def visitParenthesized(self, ctx:SelectionLanguageParser.ParenthesizedContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#query.
    def visitQuery(self, ctx:SelectionLanguageParser.QueryContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#xor.
    def visitXor(self, ctx:SelectionLanguageParser.XorContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#bonded_to.
    def visitBonded_to(self, ctx:SelectionLanguageParser.Bonded_toContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#around.
    def visitAround(self, ctx:SelectionLanguageParser.AroundContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#beyond.
    def visitBeyond(self, ctx:SelectionLanguageParser.BeyondContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#simple_expr.
    def visitSimple_expr(self, ctx:SelectionLanguageParser.Simple_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#all.
    def visitAll(self, ctx:SelectionLanguageParser.AllContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#none.
    def visitNone(self, ctx:SelectionLanguageParser.NoneContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#pdb_expr.
    def visitPdb_expr(self, ctx:SelectionLanguageParser.Pdb_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#chain.
    def visitChain(self, ctx:SelectionLanguageParser.ChainContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#residue_name.
    def visitResidue_name(self, ctx:SelectionLanguageParser.Residue_nameContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#residue_number.
    def visitResidue_number(self, ctx:SelectionLanguageParser.Residue_numberContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#alternate_location.
    def visitAlternate_location(self, ctx:SelectionLanguageParser.Alternate_locationContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#atom_serial.
    def visitAtom_serial(self, ctx:SelectionLanguageParser.Atom_serialContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#b_factor.
    def visitB_factor(self, ctx:SelectionLanguageParser.B_factorContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#occupancy.
    def visitOccupancy(self, ctx:SelectionLanguageParser.OccupancyContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#hetatm.
    def visitHetatm(self, ctx:SelectionLanguageParser.HetatmContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#atom_property.
    def visitAtom_property(self, ctx:SelectionLanguageParser.Atom_propertyContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#property_expr.
    def visitProperty_expr(self, ctx:SelectionLanguageParser.Property_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#index.
    def visitIndex(self, ctx:SelectionLanguageParser.IndexContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#element.
    def visitElement(self, ctx:SelectionLanguageParser.ElementContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#atom_name.
    def visitAtom_name(self, ctx:SelectionLanguageParser.Atom_nameContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#formal_charge.
    def visitFormal_charge(self, ctx:SelectionLanguageParser.Formal_chargeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#degree.
    def visitDegree(self, ctx:SelectionLanguageParser.DegreeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#heavy_degree.
    def visitHeavy_degree(self, ctx:SelectionLanguageParser.Heavy_degreeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#valence.
    def visitValence(self, ctx:SelectionLanguageParser.ValenceContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#hydrogen_count.
    def visitHydrogen_count(self, ctx:SelectionLanguageParser.Hydrogen_countContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#ring_size.
    def visitRing_size(self, ctx:SelectionLanguageParser.Ring_sizeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#molecule.
    def visitMolecule(self, ctx:SelectionLanguageParser.MoleculeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#membership_expr.
    def visitMembership_expr(self, ctx:SelectionLanguageParser.Membership_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#aromatic.
    def visitAromatic(self, ctx:SelectionLanguageParser.AromaticContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#aliphatic.
    def visitAliphatic(self, ctx:SelectionLanguageParser.AliphaticContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#halogen.
    def visitHalogen(self, ctx:SelectionLanguageParser.HalogenContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#metal.
    def visitMetal(self, ctx:SelectionLanguageParser.MetalContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#in_ring.
    def visitIn_ring(self, ctx:SelectionLanguageParser.In_ringContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#water.
    def visitWater(self, ctx:SelectionLanguageParser.WaterContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#backbone.
    def visitBackbone(self, ctx:SelectionLanguageParser.BackboneContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#alpha_carbon.
    def visitAlpha_carbon(self, ctx:SelectionLanguageParser.Alpha_carbonContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#chiral.
    def visitChiral(self, ctx:SelectionLanguageParser.ChiralContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#smarts_expr.
    def visitSmarts_expr(self, ctx:SelectionLanguageParser.Smarts_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#smarts.
    def visitSmarts(self, ctx:SelectionLanguageParser.SmartsContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#openeye_atombondset_expr.
    def visitOpeneye_atombondset_expr(self, ctx:SelectionLanguageParser.Openeye_atombondset_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#openeye_atom_expr.
    def visitOpeneye_atom_expr(self, ctx:SelectionLanguageParser.Openeye_atom_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#openeye_bond_expr.
    def visitOpeneye_bond_expr(self, ctx:SelectionLanguageParser.Openeye_bond_exprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#string_selector.
    def visitString_selector(self, ctx:SelectionLanguageParser.String_selectorContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#integer_selector.
    def visitInteger_selector(self, ctx:SelectionLanguageParser.Integer_selectorContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#float_selector.
    def visitFloat_selector(self, ctx:SelectionLanguageParser.Float_selectorContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#integer_slices.
    def visitInteger_slices(self, ctx:SelectionLanguageParser.Integer_slicesContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#integer_comparison.
    def visitInteger_comparison(self, ctx:SelectionLanguageParser.Integer_comparisonContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#integer.
    def visitInteger(self, ctx:SelectionLanguageParser.IntegerContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#integer_range.
    def visitInteger_range(self, ctx:SelectionLanguageParser.Integer_rangeContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#float_slices.
    def visitFloat_slices(self, ctx:SelectionLanguageParser.Float_slicesContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#float_comparison.
    def visitFloat_comparison(self, ctx:SelectionLanguageParser.Float_comparisonContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#float.
    def visitFloat(self, ctx:SelectionLanguageParser.FloatContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionLanguageParser#float_range.
    def visitFloat_range(self, ctx:SelectionLanguageParser.Float_rangeContext):
        return self.visitChildren(ctx)



del SelectionLanguageParser