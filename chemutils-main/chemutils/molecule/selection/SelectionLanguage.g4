// ANTLR Grammar file defining the selection language used by chemutils.

grammar SelectionLanguage;

// Full expression, consisting of a single selection expresion followed by the end of the input.
full_expr: WS? expr WS? EOF;

// An expression wrapped in parentheses.
parenthesized_expr: OPEN_PAREN WS? expr WS? CLOSE_PAREN;

// A general expression
// The order here defines the precedence of the operators (for example, `and` before `or`).
expr: 
    parenthesized_expr  #parenthesized
    // `not {expr}` to invert a selection
    | NOT WS expr #not
    // Simple queries like `chain A`
    | simple_expr #query
    // `same chain as {expr}` to expand a selection to a chain
    | SAME WS atom_property WS AS WS expr #same
    // `bonded to {expr}` to select atoms bonded to the selection
    | BONDED WS TO WS expr #bonded_to
    // `within {dist} of {expr2}` for all atoms within `dist` angstroms of `expr2`.
    | WITHIN WS float WS OF WS expr  #within
    // `around {dist} of {expr2}` for all atoms within `dist` angstroms of `expr2`, but not in `{expr2}`.
    | AROUND WS float WS OF WS expr  #around
    // `beyond {dist} of {expr2}` for all atoms beyond `dist` angstroms of `expr2`.
    | BEYOND WS float WS OF WS expr  #beyond
    // `{expr1} and {expr2} and ...` for conditional and.
    | expr WS (AND WS expr)+  #and
    // `{expr1} xor {expr2} xor ...` for conditional xor.
    | expr WS (XOR WS expr)+  #xor
    // `{expr1} or {expr2} or ...` for conditional or.
    | expr WS (OR WS expr)+  #or
    ;

// Simple expression (non-condtionals)
simple_expr: all | none | property_expr | pdb_expr | membership_expr | smarts_expr | openeye_atom_expr | openeye_atombondset_expr;

// `all` expression that selects everything
all: ALL;
// `none` expression that selects nothing
none: NONE;

// Expressions based on PDB fields
pdb_expr: chain | residue_name | residue_number | alternate_location | atom_serial | b_factor | occupancy | hetatm;

// `chain {chain1} {chain2} ...` matches atoms with the specific chain IDs
// `chain` matches any atom with a chain ID
chain: CHAIN (WS string_selector)?;
// `resname {resname1} {resname2} ...` matches atoms with the specific residue names
// `resname` matches any atom with a residue name
residue_name: RESIDUE_NAME (WS string_selector)?;
// `resnum {resnum1} {resnum2} ...` matches atoms with the specific residue numbers
// `resnum` matches any atom with a residue number
residue_number: RESIDUE_NUMBER (WS integer_selector)?;
// `altloc {alt1} {alt2} ...` matches atoms with the specific alternate locations
// `altloc` matches any atom with an alternate location
alternate_location: ALTERNATE_LOCATION (WS string_selector)?;
// `serial {serial1} {serial2} ...` matches atoms with the specific serial numbers. This also supports ranges and comparisons.
atom_serial: ATOM_SERIAL (WS integer_selector)?;
// `b_factor {b_factor1} {b_factor2} ...` matches atoms with specific B-factors. This also supports ranges and comparisons.
b_factor: B_FACTOR (WS float_selector)?;
// `occupancy {occupancy1} {occupancy2} ...` matches atoms with specific occupancies. This also supports ranges and comparisons.
occupancy: OCCUPANCY (WS float_selector)?;
// `hetatm` matches any atom that originates from a HETATM record
hetatm: HETATM;

atom_property: CHAIN | RESIDUE | MOLECULE | RING_SYSTEM | ATOM_NAME | ELEMENT | FORMAL_CHARGE | DEGREE | HEAVY_DEGREE | VALENCE | HYDROGEN_COUNT | RESIDUE_NAME | RESIDUE_NUMBER | ALTERNATE_LOCATION | B_FACTOR | OCCUPANCY;

// Expressions based on atom properties
property_expr: index | element | atom_name | formal_charge | degree | heavy_degree | valence | hydrogen_count | ring_size | molecule;

// `index {index1} {index2} ...` matches atoms with the specific indices. This also supports ranges and comparisons.
index: INDEX WS integer_selector;
element: ELEMENT (WS IDENTIFIER)+;
atom_name: ATOM_NAME (WS string_selector)?;
formal_charge: FORMAL_CHARGE WS integer_selector;
degree: DEGREE WS integer_selector;
heavy_degree: HEAVY_DEGREE WS integer_selector;
valence: VALENCE WS integer_selector;
hydrogen_count: HYDROGEN_COUNT WS integer_selector;
// `ring_size {size1} {size2} ...` matches atoms in rings of the specific sizes. This also supports ranges.
ring_size: RING_SIZE WS integer_slices;
molecule: MOLECULE WS smarts;

membership_expr: aromatic | aliphatic | halogen | metal | in_ring | water | backbone | alpha_carbon | chiral;

aromatic: AROMATIC;
aliphatic: ALIPHATIC;
halogen: HALOGEN;
metal: METAL;
in_ring: IN_RING;
water: WATER;
backbone: BACKBONE;
alpha_carbon: ALPHA_CARBON;
chiral: CHIRAL;

smarts_expr: smarts;

smarts: IDENTIFIER | SMARTS;

openeye_atombondset_expr: NUMERIC_INDEX WS openeye_atom_expr (WS openeye_atom_expr)* (WS NUMERIC_INDEX WS openeye_bond_expr (WS openeye_bond_expr)*)?;
openeye_atom_expr: IDENTIFIER | OPENEYE_ATOM_IDENTIFIER;
openeye_bond_expr: OPENEYE_BOND_IDENTIFIER;

// General selector for a number
// Can be list of numbers & ranges (`1 3 7`, `1-3 5 8-11`) or a comparison (`= 3`, `> 5`, `<= 9`)
string_selector: (IDENTIFIER | POSITIVE_INTEGER | STRING) (WS (IDENTIFIER | POSITIVE_INTEGER | STRING))*;
integer_selector: integer_slices | integer_comparison;
float_selector: integer_selector | float_slices | float_comparison;

// Combination of 1 or more integers or integer ranges
integer_slices: (integer | integer_range) (WS (integer | integer_range))*;
integer_comparison: OPERATOR WS integer; 
integer: INTEGER | POSITIVE_INTEGER;
integer_range: INTEGER_RANGE;

// Combination of 1 or more floats or float ranges
float_slices: (float | float_range) (WS (float | float_range))*;
float_comparison: OPERATOR WS float | INTEGER;
float: FLOAT | POSITIVE_FLOAT | INTEGER | POSITIVE_INTEGER;
float_range: FLOAT_RANGE | INTEGER_RANGE;

// Keywords
AROUND: 'around';
WITHIN: 'within';
BEYOND: 'beyond';
SAME: 'same';
BONDED: 'bonded';
OF: 'of';
AS: 'as';
TO: 'to';
ALL: 'all';
NONE: 'none';
CHAIN: 'chain';
RESIDUE: 'residue';
MOLECULE: 'molecule';
RESIDUE_NAME: 'resn' | 'resname';
RESIDUE_NUMBER: 'resi' | 'resnum';
ALTERNATE_LOCATION: 'alt' | 'altloc';
ATOM_SERIAL: 'serial';
B_FACTOR: 'b_factor';
OCCUPANCY: 'occ';
HETATM: 'hetatm';
INDEX: 'index';
ELEMENT: 'elem' | 'element';
ATOM_NAME: 'name';
FORMAL_CHARGE: 'formal_charge';
CHARGE: 'charge';
HEAVY_DEGREE: 'hvy_degree';
DEGREE: 'degree';
VALENCE: 'valence';
HYDROGEN_COUNT: 'hcount';
RING_SIZE: 'ring_size';
RING_SYSTEM: 'ring_system';
AROMATIC: 'aromatic';
ALIPHATIC: 'aliphatic';
HALOGEN: 'halogen';
METAL: 'metal';
IN_RING: 'ring';
WATER: 'water';
BACKBONE: 'backbone';
CHIRAL: 'chiral';
ALPHA_CARBON: 'alpha_carbon';
NOT  : 'not';
OR   : 'or';
AND  : 'and';
XOR  : 'xor';
OPERATOR  : '!=' | '=' | '<' | '>' | '<=' | '>=';
WS              :   [ \t\r\n]+;
OPEN_PAREN: '(';
CLOSE_PAREN: ')';

FLOAT_RANGE: ((POSITIVE_FLOAT | POSITIVE_INTEGER) '-' POSITIVE_FLOAT) |  (POSITIVE_FLOAT '-' (POSITIVE_FLOAT | POSITIVE_INTEGER)) | ((FLOAT | INTEGER) ':' FLOAT) | (FLOAT ':' (FLOAT | INTEGER));
INTEGER_RANGE: (POSITIVE_INTEGER '-' POSITIVE_INTEGER) | (INTEGER ':' INTEGER);

FLOAT: '-'? [0-9]+ '.' [0-9]+;
POSITIVE_FLOAT: '-'? [0-9]+ '.' [0-9]+;
POSITIVE_INTEGER: [0-9]+;
INTEGER: '-'? [0-9]+;

NUMERIC_INDEX: '[' INTEGER ']';

OPENEYE_ATOM_IDENTIFIER: [0-9]+ WS [A-Z] [a-z]?;
OPENEYE_BOND_IDENTIFIER: [0-9]+ [A-Z] [a-z]? '-' [0-9]+ [A-Z] [a-z]?;

IDENTIFIER: [A-Za-z0-9]+;

// Matches anything in double or single quotes.
STRING : '"' ( '\\"' | ~["] )* '"' | ('\'' ('\\\'' | ~['] )* '\'');

// Matches any valid SMILES/SMARTS string.
SMARTS: (('$(' SMARTS ')') | ('(' SMARTS ')') | ([BCNOPSFIAbcnopsa\-=#:!&,;0-9*]+ | 'Cl' | 'Br' | '[' ([A-Za-z0-9-=#+!&,;@]+ | '$(' SMARTS ')')+ ']'))+;