#!/usr/bin/env python3

"""A Pack Rat Parsing Expression Grammer parser."""

import re
from os.path import dirname, join as join_path
from textwrap import indent

EBNF_GRAMMAR = join_path(dirname(__file__), 'ebnf.ebnf')

# pylint: disable=line-too-long
EBNF_DEFS = {
    'Syntax': ('CONJUNCT', ('ONEORMORE', ('CONJUNCT', ('ZEROORMORE', ('CONJUNCT', 'EmptyLine')), 'Definition', 'newline')), ('ZEROORMORE', ('CONJUNCT', 'EmptyLine'))),
    'Definition': ('CONJUNCT', 'Identifier', 'Whitespace', '"= "', 'Expression', '";"'),
    'Expression': ('DISJUNCT', 'Disjunct', 'Except', 'Conjunct'),
    'Conjunct': ('CONJUNCT', 'Item', ('ZEROORMORE', ('CONJUNCT', '" "', 'Item'))),
    'Disjunct': ('CONJUNCT', 'Atom', ('ONEORMORE', ('CONJUNCT', 'newline', 'Whitespace', '"| "', 'Atom'))),
    'Except': ('CONJUNCT', 'Atom', ('ONEORMORE', ('CONJUNCT', 'newline', 'Whitespace', '"- "', 'Atom'))),
    'Item': ('DISJUNCT', 'Repetition', 'Atom'),
    'Repetition': ('DISJUNCT', 'ZeroOrMore', 'ZeroOrOne', 'OneOrMore'),
    'ZeroOrMore': ('CONJUNCT', '"( "', 'Conjunct', '" )*"'),
    'ZeroOrOne': ('CONJUNCT', '"( "', 'Conjunct', '" )?"'),
    'OneOrMore': ('CONJUNCT', '"( "', 'Conjunct', '" )+"'),
    'Atom': ('DISJUNCT', 'Identifier', 'Reserved', 'Literal'),
    'Identifier': ('CONJUNCT', ('ONEORMORE', ('CONJUNCT', 'upper', ('ZEROORMORE', ('CONJUNCT', 'lower'))))),
    'Reserved': ('CONJUNCT', ('ONEORMORE', ('CONJUNCT', 'lower'))),
    'Literal': ('DISJUNCT', 'DString', 'SString'),
    'DString': ('CONJUNCT', '\'"\'', ('ZEROORMORE', ('CONJUNCT', 'NoDQuote')), '\'"\''),
    'SString': ('CONJUNCT', '"\'"', ('ZEROORMORE', ('CONJUNCT', 'NoSQuote')), '"\'"'),
    'NoDQuote': ('EXCEPT', 'print', '\'"\''),
    'NoSQuote': ('EXCEPT', 'print', '"\'"'),
    'Whitespace': ('CONJUNCT', ('ONEORMORE', ('CONJUNCT', 'blank'))),
    'EmptyLine': ('CONJUNCT', ('ZEROORONE', ('CONJUNCT', '"#"', ('ZEROORMORE', ('CONJUNCT', 'print')))), 'newline'),
}


def create_parser_from_file(filepath, debug=False):
    """Create parser from a EBNF grammar file.

    Arguments:
        filepath (str): Path to EBNF file.
        debug (bool): Print debugging information. Defaults to False.

    Returns:
        PEGParser: A parser for the grammar.
    """
    with open(filepath) as fd:
        ebnf = fd.read()
    return create_parser(ebnf, debug=debug)


def create_parser(ebnf, debug=False):
    """Create parser from a EBNF grammar.

    Arguments:
        ebnf (str): A EBNF grammar.
        debug (bool): Print debugging information. Defaults to False.

    Returns:
        PEGParser: A parser for the grammar.
    """
    return PEGParser(EBNFWalker().parse(ebnf), debug=debug)


def one_line_format(string):
    """Escape tabs and spaces in a string.

    Arguments:
        string (str): The string to escape.

    Returns:
        str: The escaped string.
    """
    string = re.sub(r'\t', r'\\t', string)
    if '\n' in string:
        string = string[:string.index('\n')]
    return string


def index_to_line_col(string, index):
    """Convert an index in a string to line and column number.

    Arguments:
        string (str): The string.
        index (int): The index.

    Returns:
        int: The line number of that index.
        int: The column number of that index.
    """
    line_num = string.count('\n', 0, index) + 1
    prev_newline = string.rfind('\n', 0, index)
    if prev_newline == -1:
        column = index
    else:
        column = index - prev_newline
    return line_num, column


class ASTNode:
    """Abstract Syntax Tree (AST) node."""

    def __init__(self, term, children, string, start_pos, end_pos):
        """Initialize the ASTNode.

        The string, start_pos, and end_pos arguments matches the arguments to
        range(), ie. the substring
            string[start_pos:end_pos]

        Arguments:
            term (str): The term this node matches.
            children (list[ASTNode]): Children nodes of this term in the grammar.
            string (str): The complete string being parsed.
            start_pos (int): Index of the first character matched by this node.
            end_pos (int): Index of the last character matched by this node.
        """
        self.term = term
        self.children = children
        self.string = string
        self.start_pos = start_pos
        self.end_pos = end_pos

    @property
    def match(self):
        """Return the substring matched by this node."""
        return self.string[self.start_pos:self.end_pos]

    def first_descendant(self, path=None):
        """Get the first ASTNode descendant that matches the path.

        See the docstring for descendants() for a description of the path
        argument.

        Arguments:
            path (str): The path to the desired descendant.

        Returns:
            ASTNode: The node representing the descendant.
        """
        path = ('*' if path is None else path).split('/')
        result = self
        for term in path:
            if term == '*':
                result = result.children[0]
            else:
                children = tuple(
                    child for child in result.children
                    if child.term == term
                )
                if children:
                    result = children[0]
                else:
                    return ()
        return result

    def descendants(self, path=None):
        """Get all ASTNode descendants that match the path.

        The path describes the term of each descendant separated by a '/'.
        Where the term is irrelevant, it can either be represented by '*' or
        left empty. For example, given the following grammar:

            Expression = Operand ( Operator Operand )*;
            Operand = ParenExpression
                    | Number;
            ParenExpression = "(" Expression ")";
            Operator = "+"
                     | "-"
                     | "*"
                     | "/";
            Number = ( digit )+;

        And an ASTNode representing the following Expression:

            (1+(2*3))-5

        The call `.descendants('Operand/ParenExpression/Expression/Operand')`
        would return the Operand ASTNodes with for '1' and '(2*3)'. Both of the
        paths 'Operand' and 'Operand/*' would give the ASTNodes for '(1+(2*3))'
        and '5', but while the first path would return two Operand ASTNodes,
        the second path would return a ParenExpression and a Number ASTNode.

        Arguments:
            path (str): The path to the desired descendant.

        Returns:
            ASTNode: All descendant ASTNodes that match the path.
        """
        path = ('*' if path is None else path).split('/')
        cur_gen = (self, )
        for term in path:
            next_gen = []
            for adult in cur_gen:
                next_gen.extend(adult.children)
            if term == '*':
                cur_gen = tuple(next_gen)
            else:
                cur_gen = tuple(child for child in next_gen if child.term == term)
        return cur_gen

    def pretty_print(self, indent_level=0):
        """Print the ASTNode using indentation to denote ancestry."""
        print('{}{}: {}'.format(indent_level * 4 * ' ', self.term, one_line_format(self.match)))
        for child in self.children:
            child.pretty_print(indent_level + 1)


class PEGParser:
    """Parser for Parsing Expression Grammars (PEGs).

    A fairly standard packrat parser. The core definitions are stored as
    constants, while the custom definitions are provided to the constructor.
    """

    CORE_DEFS = {
        'empty': r'',
        'blank': r'[ \t]',
        'digit': r'[0-9]',
        'upper': r'[A-Z]',
        'lower': r'[a-z]',
        'alpha': r'[A-Za-z]',
        'alnum': r'[0-9A-Za-z]',
        'punct': r"[-!\"#$%&'()*+,./:;<=>?@[\\\]^_`{|}~]",
        'print': r'[ -~]',
        'unicode': r'[^\x00-\x7F]',
        'newline': r'\n',
        'tab': r'\t',
    }

    def __init__(self, syntax, debug=False):
        """Initialize the Parser.

        Arguments:
            syntax (dict[str]): Dictionary of term definitions. This is usually
                produced by an EBNFWalker instance.
            debug (bool): Whether to print parsing information.
                Defaults to False.
        """
        self.custom_defs = syntax
        self.debug = debug
        self.cache = {}
        self.depth = 0
        self.trace = []
        self.max_position = 0

    def parse_file(self, filepath, term):
        """Parse the contents of a file as a given term.

        Arguments:
            filepath (str): The path to the file.
            term (str): The term to parse the string as.

        Returns:
            ASTNode: The root node of the AST.
        """
        with open(filepath) as fd:
            return self.parse(fd.read(), term)

    def parse(self, string, term):
        """Parse a string as a given term.

        Arguments:
            string (str): The string to parse.
            term (str): The term to parse the string as.

        Returns:
            ASTNode: The root node of the AST.
        """
        ast, parsed = self.parse_partial(string, term)
        if parsed == len(string):
            return ast
        else:
            return self._fail_parse(string, parsed)

    def parse_partial(self, string, term):
        """Parse a string as a given term.

        Arguments:
            string (str): The string to parse.
            term (str): The term to parse the string as.

        Returns:
            ASTNode: The root node of the AST.
            int: The number of characters parsed
        """
        self.cache = {}
        self.depth = 0
        self.trace = []
        self.max_position = 0
        ast, parsed = self._dispatch(string, term, 0)
        return ast, parsed

    def _fail_parse(self, string, parsed):
        """Fail a parse by raising SyntaxError with a trace.

        Arguments:
            string (str): The string to parse.
            parsed (int): The number of characters successfully parsed.

        Raises:
            SyntaxError: If the term is not defined, or if no parse was found
                before the end of the string.
        """
        trace = []
        for position, term in reversed(self.trace):
            line, col = index_to_line_col(string, position)
            trace.append('Failed to match {} at line {} column {} (position {})'.format(term, line, col, position))
            trace.append('  ' + string.splitlines()[line - 1].replace('\t', ' '))
            trace.append('  ' + (col - 1) * '-' + '^')
        message = 'only parsed {} of {} characters:\n'.format(parsed, len(string)) + indent('\n'.join(trace), '  ')
        raise SyntaxError(message)

    def _dispatch(self, string, term, position=0):
        """Dispatch the parsing to specialized functions.

        Arguments:
            string (str): The string to parse.
            term (str): The term to parse the string as.
            position (int): The position which with to start the parse.
                Defaults to 0.

        Returns:
            ASTNode: The root node of the AST.
            int: The number of characters successfully parsed.

        Raises:
            NameError: If the terminal is unknown.
        """
        if isinstance(term, tuple) and hasattr(self, '_match_{}'.format(term[0].lower())):
            return getattr(self, '_match_{}'.format(term[0].lower()))(string, term, position)
        elif isinstance(term, str):
            ast, pos = self._get_cached(term, position)
            if ast:
                return ast, pos
            elif term in self.custom_defs:
                return self._match_custom(string, term, position)
            elif term in PEGParser.CORE_DEFS:
                return self._match_core(string, term, position)
            elif re.match(r"^'[^']*'$", term) or re.match(r'^"[^"]*"$', term):
                return self._match_literal(string, term, position)
            else:
                raise NameError('Unknown terminal {}'.format(term))
        self._debug_print('unknown non-terminal: {}'.format(term))
        return self._fail(term, position)

    def _match_zeroormore(self, string, terms, position):
        """Parse zero-or-more of a term (the * operator).

        Arguments:
            string (str): The string to parse.
            terms (list[str]): The terms to repeat, as a syntax definition.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        terms = terms[1]
        children = []
        last_pos = position
        ast, pos = self._dispatch(string, terms, position)
        while ast:
            children.extend(ast.children)
            last_pos = pos
            ast, pos = self._dispatch(string, terms, pos)
        return ASTNode('ZEROORMORE', children, string, position, last_pos), last_pos

    def _match_zeroorone(self, string, terms, position):
        """Parse zero-or-one of a term (the ? operator).

        Arguments:
            string (str): The string to parse.
            terms (list[str]): The terms to repeat, as a syntax definition.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        terms = terms[1]
        ast, pos = self._dispatch(string, terms, position)
        if ast:
            return ast, pos
        return self._dispatch(string, 'empty', position)

    def _match_oneormore(self, string, terms, position):
        """Parse one-or-more of a term (the + operator).

        Arguments:
            string (str): The string to parse.
            terms (list[str]): The terms to repeat, as a syntax definition.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        terms = terms[1]
        ast, pos = self._dispatch(string, terms, position)
        if not ast:
            return ast, pos
        children = ast.children
        last_pos = pos
        ast, pos = self._dispatch(string, terms, pos)
        while ast:
            children.extend(ast.children)
            last_pos = pos
            ast, pos = self._dispatch(string, terms, pos)
        return ASTNode('ONEORMORE', children, string, position, last_pos), last_pos

    def _match_conjunct(self, string, terms, position):
        """Parse the concatenation of multiple terms.

        Arguments:
            string (str): The string to parse.
            terms (list[str]): The terms that are concatenated.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        children = []
        pos = position
        for term in terms[1:]:
            child_ast, child_pos = self._dispatch(string, term, pos)
            if child_ast:
                if isinstance(term, tuple) and (term[0] in ('ZEROORONE', 'ZEROORMORE', 'ONEORMORE')):
                    children.extend(child_ast.children)
                else:
                    children.append(child_ast)
                pos = child_pos
            else:
                return None, child_pos
        return ASTNode('AND', children, string, position, pos), pos

    def _match_disjunct(self, string, terms, position):
        """Parse the disjunction of/any of multiple terms.

        Arguments:
            string (str): The string to parse.
            terms (list[str]): The terms to attempt to parse as.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        for term in terms[1:]:
            ast, pos = self._dispatch(string, term, position)
            if ast:
                return ast, pos
        return None, position

    def _match_except(self, string, terms, position):
        """Parse the negation of multiple terms.

        Arguments:
            string (str): The string to parse.
            terms (list[str]): The first item (index 0) is the term to match;
                all subsequent items are terms to *not* match.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        ast, pos = self._dispatch(string, terms[1], position)
        if not ast:
            return self._fail(terms[1], position)
        for term in terms[2:]:
            nast = self._dispatch(string, term, position)[0]
            if nast and ast.match == nast.match:
                return self._fail(term, position)
        return ast, pos

    def _match_custom(self, string, term, position):
        """Dispatch a parse to the custom syntax definition.

        Arguments:
            string (str): The string to parse.
            term (str): The term to parse the string as.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        expression = self.custom_defs[term]
        self._debug_print('parse called at position {} with {} >>>{}'.format(
            position, term, one_line_format(string[position:position+32])
        ))
        self.depth += 1
        ast = self._dispatch(string, expression, position)[0]
        self.depth -= 1
        if not self.trace or position > self.trace[-1][0]:
            self.trace = [(position, term)]
        else:
            self.trace.append((position, term))
        if not ast:
            return self._fail(term, position)
        if isinstance(expression, tuple) and expression[0] == 'DISJUNCT':
            ast = ASTNode(term, [ast], string, position, position + len(ast.match))
        else:
            ast.term = term
        return self._cache_and_return(term, position, ast)

    def _match_core(self, string, term, position):
        """Parse a core syntax definition.

        Arguments:
            string (str): The string to parse.
            term (str): The term to parse the string as.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        match = re.match(PEGParser.CORE_DEFS[term], string[position:])
        if match:
            ast = ASTNode(term, [], string, position, position + len(match.group(0)))
            return self._cache_and_return(term, position, ast)
        return self._fail(term, position)

    def _match_literal(self, string, term, position):
        """Parse a literal.

        Arguments:
            string (str): The string to parse.
            term (str): The literal to parse the string as.
            position (int): The position which with to start the parse.

        Returns:
            ASTNode: The root node of this abstract syntax sub-tree.
            int: The index of the last character parsed.
        """
        if string[position:].find(term[1:-1]) == 0:
            ast = ASTNode(term, [], string, position, position + len(term[1:-1]))
            return self._cache_and_return(term, position, ast)
        return self._fail(term, position)

    def _fail(self, term, position):
        """Fails a particular parse to allow backtracking.

        Arguments:
            term (str): The literal to parse the string as.
            position (int): The position which with to start the parse.

        Returns:
            None: The absence of an ASTNode.
            int: The index of the last character parsed.
        """
        if term in self.custom_defs:
            self._debug_print('failed to match {} at position {}'.format(term, position))
        return None, position

    def _cache_and_return(self, term, position, ast):
        """Cache a successful parse and return the result.

        Arguments:
            term (str): The literal to parse the string as.
            position (int): The position which with to start the parse.
            ast (ASTNode): The root of the parsed abstract syntax sub-tree.

        Returns:
            None: The absence of an ASTNode.
            int: The index of the last character parsed.
        """
        self.cache[(term, position)] = ast
        return self._get_cached(term, position)

    def _get_cached(self, term, position):
        """Retrieve a parse from cache, if it exists.

        Arguments:
            term (str): The literal to parse the string as.
            position (int): The position which with to start the parse.

        Returns:
            None: The absence of an ASTNode.
            int: The index of the last character parsed.
        """
        if (term, position) in self.cache:
            if term in self.custom_defs:
                self._debug_print('matched {} at position {}'.format(term, position))
            ast = self.cache[(term, position)]
            return ast, position + len(ast.match)
        return None, position

    def _debug_print(self, obj):
        """Print debugging information with indentation."""
        if self.debug:
            print('    ' * self.depth + str(obj))


class ASTWalker:
    """A traversal of an AST.

    This is a base class for any processing that requires the bottom-up
    building of structures from an AST. Subclass functions with the name
    `parse_Term` - where Term is the name of the node in the grammar - are
    called when a that term is encountered on the way back up the tree, ie.
    in a post-order traversal. Each such function takes two arguments:

    * ast (ASTNode): The AST rooted at that node.
    * results (list[any]): The results from the descendants of the node.
    """

    PARSE_FUNCTION_PREFIX = '_parse_'

    class EmptySentinel:
        """Sentinel to indicate that no processing was done."""

        pass

    def __init__(self, parser, root_term):
        """Initialize the traversal.

        Arguments:
            parser (ArgumentParser): The parser to use.
            root_term (str): The term to start parsing on.
        """
        self.parser = parser
        self.root_term = root_term
        self._terms_to_expand = set(
            term[len('_parse_'):] for term in dir(self)
            if term.startswith(self.PARSE_FUNCTION_PREFIX)
        )
        noskips = list(self._terms_to_expand)
        while noskips:
            noskip = noskips.pop()
            for term, definition in self.parser.custom_defs.items():
                if term in self._terms_to_expand:
                    continue
                if ASTWalker.term_in_definition(noskip, definition):
                    noskips.append(term)
                    self._terms_to_expand.add(term)

    def _postorder_traversal(self, ast, depth=0):
        """Traverses the AST in post-order.

        Arguments:
            ast (ASTNode): The AST to traverse.
            depth (int): The current depth, for printing purposes.
                Defaults to 0.

        Returns:
            any: Whatever the parse_* functions return, or an EmptySentinel.
        """
        results = []
        for child in ast.descendants('*'):
            if child.term not in self._terms_to_expand:
                continue
            result, parsed = self._postorder_traversal(child, depth=depth + 1)
            if not isinstance(result, ASTWalker.EmptySentinel):
                if parsed:
                    results.append(result)
                else:
                    results.extend(result)
        function = self.PARSE_FUNCTION_PREFIX + ast.term
        if hasattr(self, function):
            return getattr(self, function)(ast, tuple(results)), True
        elif results:
            return tuple(results), False
        else:
            return ASTWalker.EmptySentinel(), False

    def parse_file(self, filepath, term=None):
        """Parse a file with the traversal.

        Arguments:
            filepath (str): The path to the file.
            term (str): The term to start parsing on. Defaults to the term from
                the constructor.

        Returns:
            any: Whatever the parse_* functions return.
        """
        with open(filepath) as fd:
            return self.parse(fd.read(), term)

    def parse(self, text, term=None):
        """Parse a string with the traversal.

        Arguments:
            text (str): The text to parse.
            term (str): The term to start parsing on. Defaults to the term from
                the constructor.

        Returns:
            any: Whatever the parse_* functions return.
        """
        if term is None:
            term = self.root_term
        ast = self.parser.parse(text, term)
        return self.parse_ast(ast)

    def parse_partial(self, text, term=None):
        """Parse a string with the traversal.

        Arguments:
            text (str): The text to parse.
            term (str): The term to start parsing on. Defaults to the term from
                the constructor.

        Returns:
            any: Whatever the parse_* functions return.
            int: The number of characters parsed
        """
        if term is None:
            term = self.root_term
        ast, parsed = self.parser.parse_partial(text, term)
        return self.parse_ast(ast), parsed

    def parse_ast(self, ast):
        """Parse an AST with the traversal.

        Arguments:
            ast (ASTNode): The AST to parse.

        Returns:
            any: Whatever the parse_* functions return.
        """
        if ast is None:
            return None
        else:
            return self._postorder_traversal(ast)[0]

    @staticmethod
    def term_in_definition(term, definition):
        """Determine if a definition could ever expand to include a term.

        Arguments:
            term (str): The term to find.
            definition (dict[str]): Dictionary of term definitions.

        Returns:
            bool: Whether the term could be in the definition.
        """
        return any(
            (term == element or (isinstance(element, tuple) and ASTWalker.term_in_definition(term, element)))
            for element in definition
        )


class EBNFWalker(ASTWalker):
    # pylint: disable=invalid-name,no-self-use,unused-argument
    """A traversal of the EBNF grammar to build up term definitions."""

    def __init__(self):
        """Initialize the traversal."""
        super().__init__(PEGParser(EBNF_DEFS), 'Syntax')

    def flatten(self, ast, results):
        """Convert results into a tuple.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results to include.

        Returns:
            tuple[any]: The results as a tuple.
        """
        return tuple((ast.term.upper(), *results))

    def _parse_Syntax(self, ast, results):
        """Parse a Syntax node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            dict[str]: Dictionary of term definitions.
        """
        return dict(results)

    def _parse_Definition(self, ast, results):
        """Parse a Definition node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            tuple[str]: A single term definition.
        """
        return tuple(results)

    def _parse_Disjunct(self, ast, results):
        """Parse a Disjunct node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            tuple[str]: A disjunct definition.
        """
        return self.flatten(ast, results)

    def _parse_Except(self, ast, results):
        """Parse an Except node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            tuple[str]: An except definition.
        """
        return self.flatten(ast, results)

    def _parse_Conjunct(self, ast, results):
        """Parse a Conjunct node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            tuple[str]: A conjunct definition.
        """
        return self.flatten(ast, results)

    def _parse_Repetition(self, ast, results):
        """Parse a Repetition node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            tuple[str]: A repetition definition.
        """
        return self.flatten(ast.first_descendant('*'), results)

    def _parse_Reserved(self, ast, results):
        """Parse a Reserved node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            str: A reserved (core) term.
        """
        return ast.match

    def _parse_Identifier(self, ast, results):
        """Parse an Identifier node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            str: A custom term.
        """
        return ast.match

    def _parse_Literal(self, ast, results):
        """Parse a Literal node.

        Arguments:
            ast (ASTNode): The AST term to head the tuple.
            results (list[any]): The results from descendants.

        Returns:
            str: The literal.
        """
        return ast.match


def main():
    """Parse a grammar."""
    from argparse import ArgumentParser
    from fileinput import input as fileinput
    arg_parser = ArgumentParser()
    arg_parser.add_argument('-e', dest='expression', help='starting expression; if omitted, first defined term is used')
    arg_parser.add_argument('-g', dest='grammar', default=EBNF_GRAMMAR, help='EBNF grammar file')
    arg_parser.add_argument('-v', dest='verbose', action='store_true', help='show what the parser is doing')
    arg_parser.add_argument('file', default='-', nargs='?', help='text file to be parsed')
    args = arg_parser.parse_args()
    if args.expression:
        term = args.expression
    else:
        grammar = ''
        with open(args.grammar, 'r') as fd:
            grammar = fd.read()
        term = PEGParser(EBNF_DEFS).parse(grammar, 'Syntax').first_descendant('Definition/Identifier').match
    contents = ''.join(fileinput(files=args.file))
    parser = create_parser_from_file(args.grammar, debug=args.verbose)
    parser.parse(contents, term).pretty_print()


if __name__ == '__main__':
    main()
