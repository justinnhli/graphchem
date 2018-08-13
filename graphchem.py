#!/usr/bin/env python3
"""A script to determine reactions necessary to synthesize a molecule."""

from collections import namedtuple, Counter
from decimal import Decimal
from itertools import product as cross_product
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


class ReactionWalker(ASTWalker):
    """A parser for chemical reactions."""

    # pylint: disable = invalid-name, unused-argument, no-self-use

    def __init__(self):
        """Initialize a ReactionWalker."""
        super().__init__(create_parser_from_file('reactions.ebnf'), 'Reaction')

    def _parse_Reaction(self, ast, results):
        return Reaction(results[0], results[1])

    def _parse_MoleculeList(self, ast, results):
        return results

    def _parse_MoleculeCount(self, ast, results):
        if len(results) == 1:
            return MoleculeCount(Decimal(1), results[0])
        elif len(results) == 2:
            return MoleculeCount(*results)
        else:
            assert False
            return None

    def _parse_Molecule(self, ast, results):
        return Molecule(*results)

    def _parse_GroupCount(self, ast, results):
        if len(results) == 1:
            return GroupCount(results[0], 1)
        elif len(results) == 2:
            return GroupCount(*results)
        else:
            assert False
            return None

    def _parse_Group(self, ast, results):
        return results

    def _parse_Element(self, ast, results):
        return ast.match

    def _parse_Number(self, ast, results):
        return Decimal(ast.match)

    def _parse_Int(self, ast, results):
        return int(ast.match)


class Molecule:
    """A chemistry molecule."""

    def __init__(self, *components, name=None):
        """Initialize the Molecule.

        Arguments:
            *components (GroupCount): The chemical formula components (element
                and group counts).
            name (str): The name of the molecule, if any. Optional.
        """
        self.name = name
        self.components = components
        self._formula = None
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

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return self.formula


class Reaction:
    """A chemical reaction."""

    def __init__(self, reactants, products, energy=None):
        """Initialize the Reaction.

        Arguments:
            reactants (List[MoleculeCount]): The reactants and their amounts.
            products (List[MoleculeCount]): The products and their amounts.
            energy (int): The Gibbs free energy of the reaction.
        """
        self.reactant_counts = reactants
        self.product_counts = products
        self.energy = energy
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

    def _integerize(self):
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
    """A network of chemical reactions.

    The actual network is a directed bipartite graph of molecules (alternately,
    reactants and products) and chemical reactions.

    """

    REACTION_PARSER = ReactionWalker()

    def __init__(self, reaction_strs):
        """Initialize a ReactionNetwork.

        Arguments:
            reaction_strs (List[str]): Reactions in the network.
        """
        self.graph = DiGraph()
        self._build_graph(reaction_strs)
        self._reset_synthesis()

    def _build_graph(self, reaction_strs):
        for reaction_str in reaction_strs:
            reaction = self.REACTION_PARSER.parse(reaction_str)
            reaction_str = str(reaction)
            self.graph.add_node(reaction_str, type='reaction')
            for reactant in reaction.reactants:
                reactant_str = str(reactant)
                self.graph.add_node(
                    reactant_str,
                    type='molecule',
                    synthesis_networks=set(),
                )
                self.graph.add_edge(reactant_str, reaction_str)
            for product in reaction.products:
                product_str = str(product)
                self.graph.add_node(
                    product_str,
                    type='molecule',
                    synthesis_networks=set(),
                )
                self.graph.add_edge(reaction_str, product_str)

    def _reset_synthesis(self):
        for node in self.graph.nodes:
            self.graph.nodes[node]['synthesis_networks'] = set()

    def _all_reactants_synthesized(self, reaction):
        return all(
            len(self.graph.nodes[reactant]['synthesis_networks']) > 0
            for reactant, _ in self.graph.in_edges(reaction)
        )

    def _find_new_reactions(self, reactants):
        reactions = set()
        for reactant in reactants:
            for _, reaction in self.graph.out_edges(reactant):
                if self._all_reactants_synthesized(reaction):
                    reactions.add(reaction)
        return reactions

    def _calculate_synthesis_networks(self, reaction):
        reactant_networks = [
            self.graph.nodes[reactant]['synthesis_networks']
            for reactant in self.graph.predecessors(reaction)
        ]
        product_networks = set()
        for networks in cross_product(*reactant_networks):
            product_networks.add(frozenset(set.union({reaction}, *networks)))
        return product_networks

    def _synthesize_new_reactants(self, reactions):
        reactants = set()
        for reaction in reactions:
            new_product_networks = self._calculate_synthesis_networks(reaction)
            for _, product in self.graph.out_edges(reaction):
                product_node = self.graph.nodes[product]
                if not product_node['synthesis_networks']:
                    reactants.add(product)
                product_node['synthesis_networks'] |= new_product_networks
        return reactants

    def search(self, product, *reactants):
        """Search for sub-networks that could synthesize the product.

        Arguments:
            product (str): The molecule to synthesize.
            *reactants (str): The initial reactants.

        Yields:
            frozenset[str]: Reactions that could synthesize the product.
        """
        self._reset_synthesis()
        for reactant in reactants:
            reactant_node = self.graph.nodes[reactant]
            reactant_node['synthesis_networks'] = set([frozenset()])
        new_reactants = set(reactants)
        wave = 1
        while new_reactants:
            reactions = self._find_new_reactions(new_reactants)
            new_reactants = self._synthesize_new_reactants(reactions)
            # FIXME need to propagate synthesis networks to children
            wave += 1
        yield from self.graph.nodes[product]['synthesis_networks']

    def to_dot(self, final_product='', *initial_reactants): # pylint: disable = keyword-arg-before-vararg
        """Create a Graphviz representation of the reaction network.

        Arguments:
            final_product (str): The molecule to synthesize. Optional
            *initial_reactants (str): The initial reactants.

        Returns:
            str: A Graphviz description of the reaction network.
        """
        # colors from Bokeh palettes Set1
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']
        synthesis = final_product and initial_reactants
        # begin graphviz output
        dot = []
        dot.append('digraph {')
        dot.append('    node [fontsize=20]')
        dot.append('')
        # pre-calculate list of reactions
        reaction_strs = sorted(
            [
                node for node in self.graph.nodes
                if self.graph.nodes[node]['type'] == 'reaction'
            ],
            key=len,
        )
        # draw special nodes (all reactions, synthesis reactants and products)
        dot.append('    # PREAMBLE')
        dot.append('')
        if synthesis:
            dot.append('    # initial reactants')
            dot.append('    node [shape=oval, style=filled, color="#73D216"]')
            for reactant in initial_reactants:
                dot.append(f'    "{reactant}" [label=<<b>{reactant}</b>>]')
            dot.append('')
            dot.append('    # final product')
            dot.append('    node [shape=octagon, style=filled, color="#EF2929"]')
            dot.append(f'    "{final_product}" [label=<<b>{final_product}</b>>]')
            dot.append('')
        dot.append('    # all reactions')
        dot.append('    node [shape=box, style=solid, color="#000000"]')
        for reaction in reaction_strs:
            dot.append(f'    "{reaction}"')
        dot.append('')
        dot.append('    node [shape=none]')
        dot.append('')
        # color code synthesis networks
        used_reactions = set()
        if synthesis:
            dot.append('    # SYNTHESIS NETWORKS')
            dot.append('')
            syn_nets = self.search(final_product, *initial_reactants)
            syn_nets = sorted(syn_nets, key=(lambda net: (len(net), sum(len(eq) for eq in net))))
            for syn_index, syn_net in enumerate(syn_nets, start=1):
                dot.append(f'    # network {syn_index} reactions:')
                for reaction_str in syn_net:
                    dot.append(f'    #   {reaction_str}')
                for reaction_str in syn_net:
                    reaction = self.REACTION_PARSER.parse(reaction_str)
                    for reactant in reaction.reactants:
                        dot.append(f'    "{reactant}" -> "{reaction}" [color="{colors[(syn_index - 1) % len(colors)]}"]')
                    for product in reaction.products:
                        dot.append(f'    "{reaction}" -> "{product}" [color="{colors[(syn_index - 1) % len(colors)]}"]')
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
            for reactant_str in self.graph.predecessors(reaction_str):
                dot.append(f'    "{reactant_str}" -> "{reaction_str}"')
            for product_str in self.graph.successors(reaction_str):
                dot.append(f'    "{reaction_str}" -> "{product_str}"')
            dot.append('')
        dot.append('}')
        return '\n'.join(dot)


def main():
    """Example synthesis of methanol."""
    reactions = dedent("""
        CO2 + H2 = HCOOH
        CO + H2O = HCOOH
        HCOOH + H2 = CH2O + H2O
        CH2O + H2 = CH3OH
        CH3OH + H2 = CH4 + H2O
        CO + 2 H2 = CH3OH
    """).strip().splitlines()
    reactions = [reaction for reaction in reactions if reaction]
    network = ReactionNetwork(reactions)

    initial_reactants = ['H2', 'CO', 'H2O']
    final_product = 'CH3OH'

    print(network.to_dot(final_product, *initial_reactants))
    print()

    syn_nets = network.search(final_product, *initial_reactants)
    syn_nets = sorted(syn_nets, key=(lambda net: (len(net), sum(len(eq) for eq in net))))
    print(f'{len(syn_nets)} synthesis networks found:')
    print()
    for index, syn_net in enumerate(syn_nets, start=1):
        print(f'network {index}:')
        for reaction in sorted(syn_net, key=len):
            print(f'    {reaction}')
        print()


if __name__ == '__main__':
    main()
