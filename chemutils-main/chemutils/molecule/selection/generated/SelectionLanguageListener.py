# Generated from SelectionLanguage.g4 by ANTLR 4.13.2
from antlr4 import *
if "." in __name__:
    from .SelectionLanguageParser import SelectionLanguageParser
else:
    from SelectionLanguageParser import SelectionLanguageParser

# This class defines a complete listener for a parse tree produced by SelectionLanguageParser.
class SelectionLanguageListener(ParseTreeListener):

    # Enter a parse tree produced by SelectionLanguageParser#full_expr.
    def enterFull_expr(self, ctx:SelectionLanguageParser.Full_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#full_expr.
    def exitFull_expr(self, ctx:SelectionLanguageParser.Full_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#parenthesized_expr.
    def enterParenthesized_expr(self, ctx:SelectionLanguageParser.Parenthesized_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#parenthesized_expr.
    def exitParenthesized_expr(self, ctx:SelectionLanguageParser.Parenthesized_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#not.
    def enterNot(self, ctx:SelectionLanguageParser.NotContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#not.
    def exitNot(self, ctx:SelectionLanguageParser.NotContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#same.
    def enterSame(self, ctx:SelectionLanguageParser.SameContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#same.
    def exitSame(self, ctx:SelectionLanguageParser.SameContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#or.
    def enterOr(self, ctx:SelectionLanguageParser.OrContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#or.
    def exitOr(self, ctx:SelectionLanguageParser.OrContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#within.
    def enterWithin(self, ctx:SelectionLanguageParser.WithinContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#within.
    def exitWithin(self, ctx:SelectionLanguageParser.WithinContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#and.
    def enterAnd(self, ctx:SelectionLanguageParser.AndContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#and.
    def exitAnd(self, ctx:SelectionLanguageParser.AndContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#parenthesized.
    def enterParenthesized(self, ctx:SelectionLanguageParser.ParenthesizedContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#parenthesized.
    def exitParenthesized(self, ctx:SelectionLanguageParser.ParenthesizedContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#query.
    def enterQuery(self, ctx:SelectionLanguageParser.QueryContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#query.
    def exitQuery(self, ctx:SelectionLanguageParser.QueryContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#xor.
    def enterXor(self, ctx:SelectionLanguageParser.XorContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#xor.
    def exitXor(self, ctx:SelectionLanguageParser.XorContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#bonded_to.
    def enterBonded_to(self, ctx:SelectionLanguageParser.Bonded_toContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#bonded_to.
    def exitBonded_to(self, ctx:SelectionLanguageParser.Bonded_toContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#around.
    def enterAround(self, ctx:SelectionLanguageParser.AroundContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#around.
    def exitAround(self, ctx:SelectionLanguageParser.AroundContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#beyond.
    def enterBeyond(self, ctx:SelectionLanguageParser.BeyondContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#beyond.
    def exitBeyond(self, ctx:SelectionLanguageParser.BeyondContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#simple_expr.
    def enterSimple_expr(self, ctx:SelectionLanguageParser.Simple_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#simple_expr.
    def exitSimple_expr(self, ctx:SelectionLanguageParser.Simple_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#all.
    def enterAll(self, ctx:SelectionLanguageParser.AllContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#all.
    def exitAll(self, ctx:SelectionLanguageParser.AllContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#none.
    def enterNone(self, ctx:SelectionLanguageParser.NoneContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#none.
    def exitNone(self, ctx:SelectionLanguageParser.NoneContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#pdb_expr.
    def enterPdb_expr(self, ctx:SelectionLanguageParser.Pdb_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#pdb_expr.
    def exitPdb_expr(self, ctx:SelectionLanguageParser.Pdb_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#chain.
    def enterChain(self, ctx:SelectionLanguageParser.ChainContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#chain.
    def exitChain(self, ctx:SelectionLanguageParser.ChainContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#residue_name.
    def enterResidue_name(self, ctx:SelectionLanguageParser.Residue_nameContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#residue_name.
    def exitResidue_name(self, ctx:SelectionLanguageParser.Residue_nameContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#residue_number.
    def enterResidue_number(self, ctx:SelectionLanguageParser.Residue_numberContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#residue_number.
    def exitResidue_number(self, ctx:SelectionLanguageParser.Residue_numberContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#alternate_location.
    def enterAlternate_location(self, ctx:SelectionLanguageParser.Alternate_locationContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#alternate_location.
    def exitAlternate_location(self, ctx:SelectionLanguageParser.Alternate_locationContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#atom_serial.
    def enterAtom_serial(self, ctx:SelectionLanguageParser.Atom_serialContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#atom_serial.
    def exitAtom_serial(self, ctx:SelectionLanguageParser.Atom_serialContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#b_factor.
    def enterB_factor(self, ctx:SelectionLanguageParser.B_factorContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#b_factor.
    def exitB_factor(self, ctx:SelectionLanguageParser.B_factorContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#occupancy.
    def enterOccupancy(self, ctx:SelectionLanguageParser.OccupancyContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#occupancy.
    def exitOccupancy(self, ctx:SelectionLanguageParser.OccupancyContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#hetatm.
    def enterHetatm(self, ctx:SelectionLanguageParser.HetatmContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#hetatm.
    def exitHetatm(self, ctx:SelectionLanguageParser.HetatmContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#atom_property.
    def enterAtom_property(self, ctx:SelectionLanguageParser.Atom_propertyContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#atom_property.
    def exitAtom_property(self, ctx:SelectionLanguageParser.Atom_propertyContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#property_expr.
    def enterProperty_expr(self, ctx:SelectionLanguageParser.Property_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#property_expr.
    def exitProperty_expr(self, ctx:SelectionLanguageParser.Property_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#index.
    def enterIndex(self, ctx:SelectionLanguageParser.IndexContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#index.
    def exitIndex(self, ctx:SelectionLanguageParser.IndexContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#element.
    def enterElement(self, ctx:SelectionLanguageParser.ElementContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#element.
    def exitElement(self, ctx:SelectionLanguageParser.ElementContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#atom_name.
    def enterAtom_name(self, ctx:SelectionLanguageParser.Atom_nameContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#atom_name.
    def exitAtom_name(self, ctx:SelectionLanguageParser.Atom_nameContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#formal_charge.
    def enterFormal_charge(self, ctx:SelectionLanguageParser.Formal_chargeContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#formal_charge.
    def exitFormal_charge(self, ctx:SelectionLanguageParser.Formal_chargeContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#degree.
    def enterDegree(self, ctx:SelectionLanguageParser.DegreeContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#degree.
    def exitDegree(self, ctx:SelectionLanguageParser.DegreeContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#heavy_degree.
    def enterHeavy_degree(self, ctx:SelectionLanguageParser.Heavy_degreeContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#heavy_degree.
    def exitHeavy_degree(self, ctx:SelectionLanguageParser.Heavy_degreeContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#valence.
    def enterValence(self, ctx:SelectionLanguageParser.ValenceContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#valence.
    def exitValence(self, ctx:SelectionLanguageParser.ValenceContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#hydrogen_count.
    def enterHydrogen_count(self, ctx:SelectionLanguageParser.Hydrogen_countContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#hydrogen_count.
    def exitHydrogen_count(self, ctx:SelectionLanguageParser.Hydrogen_countContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#ring_size.
    def enterRing_size(self, ctx:SelectionLanguageParser.Ring_sizeContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#ring_size.
    def exitRing_size(self, ctx:SelectionLanguageParser.Ring_sizeContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#molecule.
    def enterMolecule(self, ctx:SelectionLanguageParser.MoleculeContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#molecule.
    def exitMolecule(self, ctx:SelectionLanguageParser.MoleculeContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#membership_expr.
    def enterMembership_expr(self, ctx:SelectionLanguageParser.Membership_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#membership_expr.
    def exitMembership_expr(self, ctx:SelectionLanguageParser.Membership_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#aromatic.
    def enterAromatic(self, ctx:SelectionLanguageParser.AromaticContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#aromatic.
    def exitAromatic(self, ctx:SelectionLanguageParser.AromaticContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#aliphatic.
    def enterAliphatic(self, ctx:SelectionLanguageParser.AliphaticContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#aliphatic.
    def exitAliphatic(self, ctx:SelectionLanguageParser.AliphaticContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#halogen.
    def enterHalogen(self, ctx:SelectionLanguageParser.HalogenContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#halogen.
    def exitHalogen(self, ctx:SelectionLanguageParser.HalogenContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#metal.
    def enterMetal(self, ctx:SelectionLanguageParser.MetalContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#metal.
    def exitMetal(self, ctx:SelectionLanguageParser.MetalContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#in_ring.
    def enterIn_ring(self, ctx:SelectionLanguageParser.In_ringContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#in_ring.
    def exitIn_ring(self, ctx:SelectionLanguageParser.In_ringContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#water.
    def enterWater(self, ctx:SelectionLanguageParser.WaterContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#water.
    def exitWater(self, ctx:SelectionLanguageParser.WaterContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#backbone.
    def enterBackbone(self, ctx:SelectionLanguageParser.BackboneContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#backbone.
    def exitBackbone(self, ctx:SelectionLanguageParser.BackboneContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#alpha_carbon.
    def enterAlpha_carbon(self, ctx:SelectionLanguageParser.Alpha_carbonContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#alpha_carbon.
    def exitAlpha_carbon(self, ctx:SelectionLanguageParser.Alpha_carbonContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#chiral.
    def enterChiral(self, ctx:SelectionLanguageParser.ChiralContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#chiral.
    def exitChiral(self, ctx:SelectionLanguageParser.ChiralContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#smarts_expr.
    def enterSmarts_expr(self, ctx:SelectionLanguageParser.Smarts_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#smarts_expr.
    def exitSmarts_expr(self, ctx:SelectionLanguageParser.Smarts_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#smarts.
    def enterSmarts(self, ctx:SelectionLanguageParser.SmartsContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#smarts.
    def exitSmarts(self, ctx:SelectionLanguageParser.SmartsContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#openeye_atombondset_expr.
    def enterOpeneye_atombondset_expr(self, ctx:SelectionLanguageParser.Openeye_atombondset_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#openeye_atombondset_expr.
    def exitOpeneye_atombondset_expr(self, ctx:SelectionLanguageParser.Openeye_atombondset_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#openeye_atom_expr.
    def enterOpeneye_atom_expr(self, ctx:SelectionLanguageParser.Openeye_atom_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#openeye_atom_expr.
    def exitOpeneye_atom_expr(self, ctx:SelectionLanguageParser.Openeye_atom_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#openeye_bond_expr.
    def enterOpeneye_bond_expr(self, ctx:SelectionLanguageParser.Openeye_bond_exprContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#openeye_bond_expr.
    def exitOpeneye_bond_expr(self, ctx:SelectionLanguageParser.Openeye_bond_exprContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#string_selector.
    def enterString_selector(self, ctx:SelectionLanguageParser.String_selectorContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#string_selector.
    def exitString_selector(self, ctx:SelectionLanguageParser.String_selectorContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#integer_selector.
    def enterInteger_selector(self, ctx:SelectionLanguageParser.Integer_selectorContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#integer_selector.
    def exitInteger_selector(self, ctx:SelectionLanguageParser.Integer_selectorContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#float_selector.
    def enterFloat_selector(self, ctx:SelectionLanguageParser.Float_selectorContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#float_selector.
    def exitFloat_selector(self, ctx:SelectionLanguageParser.Float_selectorContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#integer_slices.
    def enterInteger_slices(self, ctx:SelectionLanguageParser.Integer_slicesContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#integer_slices.
    def exitInteger_slices(self, ctx:SelectionLanguageParser.Integer_slicesContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#integer_comparison.
    def enterInteger_comparison(self, ctx:SelectionLanguageParser.Integer_comparisonContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#integer_comparison.
    def exitInteger_comparison(self, ctx:SelectionLanguageParser.Integer_comparisonContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#integer.
    def enterInteger(self, ctx:SelectionLanguageParser.IntegerContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#integer.
    def exitInteger(self, ctx:SelectionLanguageParser.IntegerContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#integer_range.
    def enterInteger_range(self, ctx:SelectionLanguageParser.Integer_rangeContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#integer_range.
    def exitInteger_range(self, ctx:SelectionLanguageParser.Integer_rangeContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#float_slices.
    def enterFloat_slices(self, ctx:SelectionLanguageParser.Float_slicesContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#float_slices.
    def exitFloat_slices(self, ctx:SelectionLanguageParser.Float_slicesContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#float_comparison.
    def enterFloat_comparison(self, ctx:SelectionLanguageParser.Float_comparisonContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#float_comparison.
    def exitFloat_comparison(self, ctx:SelectionLanguageParser.Float_comparisonContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#float.
    def enterFloat(self, ctx:SelectionLanguageParser.FloatContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#float.
    def exitFloat(self, ctx:SelectionLanguageParser.FloatContext):
        pass


    # Enter a parse tree produced by SelectionLanguageParser#float_range.
    def enterFloat_range(self, ctx:SelectionLanguageParser.Float_rangeContext):
        pass

    # Exit a parse tree produced by SelectionLanguageParser#float_range.
    def exitFloat_range(self, ctx:SelectionLanguageParser.Float_rangeContext):
        pass



del SelectionLanguageParser