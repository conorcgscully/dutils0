from typing import NoReturn

from antlr4 import BailErrorStrategy, CommonTokenStream, InputStream, Parser, RecognitionException
from openeye import oechem

from .generated.SelectionLanguageLexer import SelectionLanguageLexer
from .generated.SelectionLanguageParser import SelectionLanguageParser
from .visitor import SelectionVisitor  # type: ignore


class SelectionLanguagError(ValueError):
    pass


class ErrorStrategy(BailErrorStrategy):  # type: ignore
    def recover(self, recognizer: Parser, e: RecognitionException) -> NoReturn:
        raise e


def evaluate_selection(*, oemol: oechem.OEMolBase, query: str) -> oechem.OEAtomBondSet:
    input_stream = InputStream(query)
    lexer = SelectionLanguageLexer(input_stream)
    token_stream = CommonTokenStream(lexer)

    token_stream.fill()

    parser = SelectionLanguageParser(token_stream)
    parser._errHandler = ErrorStrategy()
    try:
        tree = parser.full_expr()
    except RecognitionException as e:
        raise SelectionLanguagError(f"Failed to parse selection language `{query}`.") from e

    return SelectionVisitor(oemol).visit(tree)
