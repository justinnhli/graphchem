#!/usr/bin/env python3
"""A script to determine reactions necessary to synthesize a molecule."""

from collections import namedtuple, Counter
from decimal import Decimal
from itertools import product as cross_product
from os.path import realpath, join as join_path, dirname
from textwrap import dedent

from networkx import DiGraph

from pegparse import ASTWalker, create_parser_from_file

MoleculeCount = namedtuple('MoleculeCount', 'count molecule')
GroupCount = namedtuple('GroupCount', 'group count')


def least_common_mulitple(*ints):
    """Find the least common multiple.

    Arguments:
        ints (int): Integers

    Returns:
        int: The least common multiple
    """
    result = ints[0]
    for number in ints[1:]:
        if result > number:
            numerator, denominator = result, number
        else:
            numerator, denominator = number, result
        remainder = numerator % denominator
        while remainder != 0:
            numerator, denominator = denominator, remainder
            remainder = numerator % denominator
        gcd = denominator
        result = (result * number) // gcd
    return result


def molecule_key(string):
    """Key for comparing molecule strings.

    Arguments:
        string (str): The molecular formula.

    Returns:
        Tuple[int, str]: The key.
    """
    return len(string), string


def reaction_key(string):
    """Key for comparing reaction equations.

    Arguments:
        string (str): The reaction equation.

    Returns:
        Tuple[int, str]: The key.
    """
    return len(string), string


def pathway_key(pathway):
    """Key for comparing synthesis pathways.

    Arguments:
        pathway (Sequence[str]): A list of reaction equations.

    Returns:
        Tuple[int, int]: The key.
    """
    return len(pathway), sum(len(reaction) for reaction in pathway)


class ReactionWalker(ASTWalker):
    """A parser for chemical reactions."""

    # pylint: disable = invalid-name, unused-argument, no-self-use

    def __init__(self):
        """Initialize a ReactionWalker."""
        super().__init__(
            create_parser_from_file(join_path(
                dirname(realpath(__file__)),
                'reactions.peg',
            )),
            'reaction',
        )

    def _parse_reaction(self, ast, results):
        return Reaction(results[0], results[1], formula=ast.match)

    def _parse_molecule_list(self, ast, results):
        return results

    def _parse_molecule_count(self, ast, results):
        if len(results) == 1:
            return MoleculeCount(Decimal(1), results[0])
        elif len(results) == 2:
            return MoleculeCount(*results)
        else:
            assert False
            return None

    def _parse_molecule(self, ast, results):
        return Molecule(*results, formula=ast.match)

    def _parse_group_count(self, ast, results):
        if len(results) == 1:
            return GroupCount(results[0], 1)
        elif len(results) == 2:
            return GroupCount(*results)
        else:
            assert False
            return None

    def _parse_group(self, ast, results):
        return results

    def _parse_element(self, ast, results):
        return ast.match

    def _parse_number(self, ast, results):
        return Decimal(ast.match)

    def _parse_int(self, ast, results):
        return int(ast.match)


class Molecule:
    """A chemistry molecule."""

    def __init__(self, *components, name=None, formula=None):
        """Initialize the Molecule.

        Arguments:
            *components (GroupCount): The chemical formula components (element
                and group counts).
            name (str): The name of the molecule, if any. Optional.
        """
        self.name = name
        self.components = components
        self._formula = formula
        self._atoms = None

    @property
    def formula(self):
        """Get the chemical formula of the molecule.

        Returns:
            str: The chemical formula of the molecule.
        """
        if self._formula is None:
            self._formula = ''.join(
                self._build_formula(component)
                for component in self.components
            )
        return self._formula

    def _build_formula(self, component):
        if isinstance(component.group, str):
            result = component.group
        else:
            result = ''.join([
                '(',
                ''.join(
                    self._build_formula(subcomponent)
                    for subcomponent in component.group
                ),
                ')',
            ])
        if component.count == 1:
            return result
        else:
            return result + str(component.count)

    @property
    def atoms(self):
        """Get a count of the atoms in the molecule.

        Returns:
            Counter[str]: A count of the atoms in the molecule.
        """
        if self._atoms is None:
            self._atoms = Counter()
            for component in self.components:
                self._atoms.update(self._build_atoms(component))
        return self._atoms

    def _build_atoms(self, component):
        if isinstance(component.group, str):
            return Counter({component.group: component.count})
        else:
            result = Counter()
            for subcomponent in component.group:
                subresult = self._build_atoms(subcomponent)
                for _ in range(component.count):
                    result.update(subresult)
            return result

    def __eq__(self, other):
        return str(self) == str(other)

    def __lt__(self, other):
        return (len(str(self)), str(self)) < (len(str(other)), str(other))

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return self.formula


class Reaction:
    """A chemical reaction."""

    def __init__(self, reactants, products, formula=None, energy=None):
        """Initialize the Reaction.

        Arguments:
            reactants (List[MoleculeCount]): The reactants and their amounts.
            products (List[MoleculeCount]): The products and their amounts.
            energy (int): The Gibbs free energy of the reaction.
        """
        self.reactant_counts = reactants
        self.product_counts = products
        self.energy = energy
        self.formula = formula
        self._integerize()
        self._check_equality()

    @property
    def reactants(self):
        """Get the reactants.

        Returns:
            List[Molecule]: The reactants of this reaction.
        """
        return [reactant.molecule for reactant in self.reactant_counts]

    @property
    def products(self):
        """Get the products.

        Returns:
            List[Molecule]: The products of this reaction.
        """
        return [product.molecule for product in self.product_counts]

    @property
    def reactant_coefficients(self):
        """Get the coefficients of each reactant.

        Returns:
            List[int]: The coefficients of each reactant in this reaction.
        """
        return [reactant.count for reactant in self.reactant_counts]

    @property
    def product_coefficients(self):
        """Get the coefficients of each product.

        Returns:
            List[int]: The coefficients of each product in this reaction.
        """
        return [product.count for product in self.product_counts]

    def __str__(self):
        return self.formula

    def _integerize(self):
        """Convert reaction equation to integers.

        This function does not reduce the coefficients to the smallest values.
        """
        multiplier = least_common_mulitple(
            *(coeff.as_integer_ratio()[1] for coeff in self.reactant_coefficients),
            *(coeff.as_integer_ratio()[1] for coeff in self.product_coefficients),
        )
        self.reactant_counts = [
            MoleculeCount(int(multiplier * count), molecule)
            for count, molecule in self.reactant_counts
        ]
        self.product_counts = [
            MoleculeCount(int(multiplier * count), molecule)
            for count, molecule in self.product_counts
        ]

    def _check_equality(self):
        """Check that atoms are conserved in this reaction."""
        reactants_counter = Counter()
        for count, molecule in self.reactant_counts:
            for _ in range(count):
                reactants_counter.update(molecule.atoms)
        products_counter = Counter()
        for count, molecule in self.product_counts:
            for _ in range(count):
                products_counter.update(molecule.atoms)
        assert reactants_counter == products_counter, \
            '\n'.join([
                'Reaction is not balanced:',
                *(
                    ' '.join([
                        3 * ' ',
                        'Reactants contain',
                        f'{reactants_counter[element]} {element}',
                        'but products contain',
                        f'{products_counter[element]}'
                    ])
                    for element in (reactants_counter.keys() | products_counter.keys())
                    if reactants_counter[element] != products_counter[element]
                ),
            ])

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return ' '.join([
            ' + '.join(
                (
                    ('' if count == 1 else f'{count} ') +
                    str(molecule)
                )
                for count, molecule in self.reactant_counts
            ),
            '=',
            ' + '.join(
                (
                    ('' if count == 1 else f'{count} ') +
                    str(molecule)
                )
                for count, molecule in self.product_counts
            ),
        ])


class ReactionNetwork:
    r"""A network of chemical reactions.

    The actual network is a directed bipartite graph of molecules (alternately,
    reactants and products) and chemical reactions.

    The search algorithm in synthesize() is modified breadth first search, and
    consists of a loop with two phases:

    1. Based on the available reactants, check which reactions will trigger.
    2. Based on these reactions, check which new reactants are available.

    One insight that helps with efficiency is noticing that, at each step, we
    only care about which _new_ reactants are available. Any reaction that does
    not use those reactants either have already triggered, or could not be
    triggered at this step.

    Dynamic programming is used to collect all the reactions necessary to
    synthesize a product. Associated with each molecule is the set of pathways
    (each itself a set of reactions) that could synthesize it. This is most
    easily explained with a recurrence relation. First, initial reactants have
    the set of an empty pathway. This is to differentiate them from reactants
    that have yet to be synthesized, which would have no pathways (the empty
    set).

    Whenever a reaction is triggered, the products then get their pathways
    updated. Because each reactant of that reaction could have multiple
    pathways, each of those possible pathways must be taken into account. For
    example, for a reaction A + B + C = D + E, and A, B, and C has 2, 3, and 3
    synthesis pathways respectively, then D and E would have 2 * 3 * 3 = 18
    synthesis pathways - the cross-product of all pathways of the reactants,
    each with the addition of the current reaction.

    For formally, we can define a molecule m and its synthesis pathways
    paths(m). Similarly, we might write that a reaction r has reactants pred(r)
    and products succ(r). Then we initialize:

    paths(m) =
        {{}}    if m is an initial reactant
        {}      otherwise

    When a reaction r is triggered, then for all m in succ(r),

    paths(m) = union(
        paths(m),
        { union(P, {r}) where
            P is in the unordered product of paths(n) where n in pred(r) }
    )

    Or in LaTeX:

    paths(m) \cup \left\{
        P  \cup \{r\} : P \in \Pi_{n \in pred(r)} paths(n)
    \right\}
    """

    REACTION_PARSER = ReactionWalker()

    def __init__(self, reaction_strs):
        """Initialize a ReactionNetwork.

        Arguments:
            reaction_strs (List[str]): Reactions in the network.
        """
        self.graph = DiGraph()
        self._reactions = {} # str -> Reaction
        self._build_graph(reaction_strs)
        self._reset_synthesis()

    @property
    def molecules(self):
        yield from (
            (node, self.graph.nodes[node])
            for node in self.graph.nodes
            if self.graph.nodes[node]['type'] == 'molecule'
        )

    @property
    def reactions(self):
        yield from (
            (node, self.graph.nodes[node])
            for node in self.graph.nodes
            if self.graph.nodes[node]['type'] == 'reaction'
        )

    def __getitem__(self, key):
        """Returns the node for the given key.

        Arguments:
            key (str): Either a molecule or a reaction.

        Returns:
            dict: Data associated with that node.
        """
        return self.graph.nodes[key]

    def _build_graph(self, reaction_strs):
        for reaction_str in reaction_strs:
            reaction = self.REACTION_PARSER.parse(reaction_str)
            reaction_str = str(reaction)
            self._reactions[reaction_str] = reaction
            self.graph.add_node(reaction_str, type='reaction', triggered=False)
            for reactant in reaction.reactants:
                reactant_str = str(reactant)
                self.graph.add_node(reactant_str, type='molecule', pathways=set())
                self.graph.add_edge(reactant_str, reaction_str)
            for product in reaction.products:
                product_str = str(product)
                self.graph.add_node(product_str, type='molecule', pathways=set())
                self.graph.add_edge(reaction_str, product_str)

    def _reset_synthesis(self):
        """Reset the graph from the previous synthesis."""
        for node in self.graph.nodes:
            self.graph.nodes[node]['pathways'] = set()

    def _all_reactants_synthesized(self, reaction):
        """Check if a reaction has all its reactants synthesized.

        Arguments:
            reaction (str): The reaction to check.

        Returns:
            bool: True all the reaction's reactants are have been synthesized.
        """
        return (
            self.graph.nodes[reaction]['triggered'] or
            all(
                len(self.graph.nodes[reactant]['pathways']) > 0
                for reactant in self.graph.predecessors(reaction)
            )
        )

    def _trigger_new_reactions(self, reactant):
        """Identify reactions that were only missing a given reactant.

        Arguments:
            reactant (str): A reactant.

        Returns:
            Set[str]: The reactions that are newly triggered.
        """
        reactions = set()
        if reactant not in self.graph.nodes:
            return reactions
        for reaction in self.graph.successors(reactant):
            if self._all_reactants_synthesized(reaction):
                self.graph.nodes[reaction]['triggered'] = True
                reactions.add(reaction)
        return reactions

    def _synthesize_new_reactants(self, reaction):
        """Mark new reactants as synthesized.

        Arguments:
            reaction (str): The reaction to check.

        Returns:
            Set[str]: The new products that are now synthesizable.
        """
        return set(
            product for product in self.graph.successors(reaction)
            if not self.graph.nodes[product]['pathways']
        )

    def _propagate_synthesis_pathways(self, reaction):
        """Propagate synthesis pathways down the network.

        Arguments:
            reaction (str): The reaction to start the propagation form.
        """
        queue = [reaction]
        while queue:
            reaction = queue.pop(0)
            new_pathways = self._chain_synthesis_pathways(reaction)
            for product in self.graph.successors(reaction):
                product_node = self.graph.nodes[product]
                if product_node['pathways'] >= new_pathways:
                    continue
                product_node['pathways'] |= new_pathways
                product_node['pathways'] = self._simplify_pathways(product_node['pathways'])
                queue.extend(self.graph.successors(product))

    def _chain_synthesis_pathways(self, reaction):
        """Determine the synthesis pathway of a reaction.

        Arguments:
            reaction (str): The reaction to check.

        Returns:
            Set[FrozenSet[str]]: All the synthesis pathways of this reaction.
        """
        reactant_networks = [
            self.graph.nodes[reactant]['pathways']
            for reactant in self.graph.predecessors(reaction)
        ]
        pathways = set()
        for networks in cross_product(*reactant_networks):
            pathways.add(frozenset(set.union({reaction}, *networks)))
        return pathways

    def _simplify_pathways(self, pathways):
        """Remove redundant pathways. TODO

        The problem is that with a sufficiently completed network, there will
        be a combinatorial explosion in the number of synthesis pathways. For
        example, imagine that there molecule A and B could both be synthesized
        in two ways. If A and B could be used to synthesize C, there would be
        four different pathways for its synthesis. Clearly, this would not
        scale as the network includes more reactions.

        Arguments:
            pathways (FrozenSet[Set[str]]): Synthesis pathways

        Returns:
            (FrozenSet[Set[str]]): The same pathways, simplified.
        """
        return pathways

    def synthesize(self, product, *reactants):
        """Search for sub-networks that could synthesize the product.

        Arguments:
            product (str): The molecule to synthesize.
            *reactants (str): The initial reactants.

        Yields:
            frozenset[str]: Reactions that could synthesize the product.
        """
        self._reset_synthesis()
        for reactant in reactants:
            if reactant in self.graph.nodes:
                self.graph.nodes[reactant]['pathways'] = set([frozenset()])
        new_reactants = set(reactants)
        while new_reactants:
            reactions = set()
            for reactant in new_reactants:
                reactions = reactions.union(self._trigger_new_reactions(reactant))
            new_reactants = set()
            for reaction in reactions:
                new_reactants = new_reactants.union(self._synthesize_new_reactants(reaction))
                self._propagate_synthesis_pathways(reaction)
        yield from self.graph.nodes[product]['pathways']

    def to_dot(self, final_product='', *initial_reactants): # pylint: disable = keyword-arg-before-vararg
        """Create a Graphviz representation of the reaction network.

        Arguments:
            final_product (str): The molecule to synthesize. Optional
            *initial_reactants (str): The initial reactants.

        Returns:
            str: A Graphviz description of the reaction network.
        """
        # pylint: disable = too-many-branches
        # from Bokeh palettes Set1
        colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999']
        synthesis = final_product and initial_reactants
        # begin graphviz output
        dot = []
        dot.append('digraph {')
        dot.append('')
        # pre-calculate list of reactions
        reaction_strs = sorted(
            [reaction for reaction, _ in self.reactions],
            key=reaction_key,
        )
        # draw special nodes (all reactions, synthesis reactants and products)
        dot.append('    # PREAMBLE')
        dot.append('    node [fontsize=20, style=filled]')
        dot.append('')
        if synthesis:
            dot.append('    # initial reactants')
            dot.append('    node [shape=oval, fillcolor="#CCFFAA"]')
            for reactant in initial_reactants:
                dot.append(f'    "{reactant}"')
            dot.append('')
            dot.append('    # final product')
            dot.append('    node [shape=octagon, fillcolor="#FFBFCF"]')
            dot.append(f'    "{final_product}"')
            dot.append('')
        dot.append('    # all reactions')
        dot.append('    node [shape=box, fillcolor="#FFFFFF"]')
        for reaction in reaction_strs:
            dot.append(f'    "{reaction}"')
        dot.append('')
        # color code synthesis networks
        used_reactions = set()
        if synthesis:
            dot.append('    # SYNTHESIS PATHWAYS')
            dot.append('    node [shape=none, style=""]')
            dot.append('')
            pathways = set(self.synthesize(final_product, *initial_reactants))
            pathways = sorted(pathways, key=pathway_key)
            for index, pathway in enumerate(pathways, start=1):
                dot.append(f'    # pathway {index}:')
                color = colors[(index - 1) % len(colors)]
                for reaction_str in sorted(pathway, key=reaction_key):
                    dot.append(f'    #   {reaction_str}')
                dot.append(f'    edge [color="{color}"]')
                for reaction_str in sorted(pathway, key=reaction_key):
                    reaction = self._reactions[reaction_str]
                    for reactant in sorted(reaction.reactants):
                        dot.append(f'    "{reactant}" -> "{reaction}"')
                    for product in sorted(reaction.products):
                        dot.append(f'    "{reaction}" -> "{product}"')
                    used_reactions.add(reaction_str)
                dot.append('')
        # color code remainder of network
        dot.append('    # REACTIONS NOT USED FOR SYNTHESIS')
        dot.append('    edge [color="#808080", style=dashed]')
        dot.append('')
        for reaction_str in reaction_strs:
            if reaction_str in used_reactions:
                continue
            dot.append(f'    # {reaction_str}')
            for reactant_str in sorted(self.graph.predecessors(reaction_str), key=reaction_key):
                dot.append(f'    "{reactant_str}" -> "{reaction_str}"')
            for product_str in sorted(self.graph.successors(reaction_str), key=reaction_key):
                dot.append(f'    "{reaction_str}" -> "{product_str}"')
            dot.append('')
        dot.append('}')
        return '\n'.join(dot)


def test():
    """Test for ReactionNetwork."""
    def pretty_string(pathways):
        lines = []
        pathways = sorted(pathways, key=pathway_key)
        for index, pathway in enumerate(pathways, start=1):
            lines.append(str(index))
            for reaction in sorted(pathway, key=reaction_key):
                lines.append(reaction)
        return '\n'.join(lines)

    reactions = dedent('''
        CO2 + H2 = HCOOH
        CO + H2O = HCOOH
        HCOOH + H2 = CH2O + H2O
        CH2O + H2 = CH3OH
        CH3OH + H2 = CH4 + H2O
        CO + 2 H2 = CH3OH
    ''').strip().splitlines()
    reactions = [reaction for reaction in reactions if reaction]
    network = ReactionNetwork(reactions)
    actual = set(network.synthesize('CH3OH', 'CO', 'H2', 'H2O'))
    expected = set([
        frozenset({
            'CO + 2 H2 = CH3OH',
        }),
        frozenset({
            'CO + H2O = HCOOH',
            'CH2O + H2 = CH3OH',
            'HCOOH + H2 = CH2O + H2O',
        }),
        frozenset({
            'CO + H2O = HCOOH',
            'CH2O + H2 = CH3OH',
            'CH3OH + H2 = CH4 + H2O',
            'HCOOH + H2 = CH2O + H2O',
        }),
        frozenset({
            'CO + H2O = HCOOH',
            'CH2O + H2 = CH3OH',
            'CO + 2 H2 = CH3OH',
            'CH3OH + H2 = CH4 + H2O',
            'HCOOH + H2 = CH2O + H2O',
        }),
    ])
    assert actual == expected, '\n\n'.join([
        'Expected:',
        pretty_string(expected),
        'but got:',
        pretty_string(actual),
    ])


SMALL_REACTION_SET = [
    'CO2 + H2 = HCOOH',
    'CO + H2O = HCOOH',
    'HCOOH + H2 = CH2O + H2O',
    'CH2O + H2 = CH3OH',
    'CH3OH + H2 = CH4 + H2O',
    'CO + 2 H2 = CH3OH',
]

LARGE_REACTION_SET = [
    'CO2 + H2 = HCOOH',
    'CO + H2O = HCOOH',
    'HCOOH + H2 = CH2O + H2O',
    'CH2O + H2 = CH3OH',
    'CH3OH + H2 = CH4 + H2O',
    '2 CO2 + 2 H2 = 2 CH2O + O2',
    '2 CO2 + 4 H2 = 2 CH3OH + O2',
    'CO2 + 4 H2 = CH4 + 2 H2O',
    '2 HCOOH + 2 H2 = 2 CH3OH + O2',
    'HCOOH + H2 = CH4 + O2',
    'CH2O + 2 H2 = CH4 + H2O',
    '2 CO + O2 = 2 CO2',
    'CO + H2 = CH2O',
    'CO + 2 H2 = CH3OH',
    'CO + 3 H2 = CH4 + H2O',
]

def main():
    """Demonstrate the synthesis of methanol."""
    reactions = SMALL_REACTION_SET
    reactions = [reaction for reaction in reactions if reaction]
    network = ReactionNetwork(reactions)

    initial_reactants = ['H2', 'CO', 'H2O']
    final_product = 'CH3OH'

    print(network.to_dot(final_product, *initial_reactants))
    print()

    pathways = network.synthesize(final_product, *initial_reactants)
    pathways = sorted(pathways, key=pathway_key)
    print(f'{len(pathways)} synthesis pathways found:')
    print()
    for index, pathway in enumerate(pathways, start=1):
        print(f'pathway {index}:')
        for reaction in sorted(pathway, key=reaction_key):
            print(f'    {reaction}')
        print()


if __name__ == '__main__':
    test()
    main()
