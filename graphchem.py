#!/usr/bin/env python3
"""A script to determine reactions necessary to synthesize a molecule."""

from argparse import ArgumentParser
from collections import namedtuple, Counter, defaultdict
from decimal import Decimal
from heapq import heappush, heappop
from os.path import realpath, join as join_path, dirname

from pegparse import ASTWalker, create_parser_from_file

MoleculeCount = namedtuple('MoleculeCount', 'count molecule')
GroupCount = namedtuple('GroupCount', 'group count')


def least_common_mulitple(*ints):
    """Find the least common multiple.

    Parameters:
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
        super().__init__(
            create_parser_from_file(join_path(
                dirname(realpath(__file__)),
                'reactions.peg',
            )),
            'reaction',
        )

    def _parse_reaction(self, ast, results):
        return Reaction(results[0], results[1])

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

    def __init__(self, *components, formula=None, name=None):
        """Initialize the Molecule.

        Parameters:
            *components (GroupCount): The chemical formula components (element
                and group counts).
            formula (str): The formula for the molecule. Optional.
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

    def __init__(self, reactants, products, energy=None):
        """Initialize the Reaction.

        Parameters:
            reactants (List[MoleculeCount]): The reactants and their amounts.
            products (List[MoleculeCount]): The products and their amounts.
            energy (int): The Gibbs free energy of the reaction. Optional.
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


def num_atom_difference(source, target):
    """Calculate the number of atoms different between two molecules.

    This function implements a naive heuristic between two molecules: the number
    of atoms that need to be added/removed from the source molecule to be
    transformed to the target molecule.

    Parameters:
        source (Molecule): The source molecule.
        target (Molecule): The target molecule.

    Returns:
        int: The number of atoms needed for the transformation.
    """
    return sum(
        abs(source.atoms.get(element, 0) - target.atoms.get(element, 0))
        for element in (set(source.atoms) | set(target.atoms))
    )


def molecular_difference(source, target):
    """Calculate the "difference" between two molecules.

    Parameters:
        source (Molecule): The source molecule.
        target (Molecule): The target molecule.

    Returns:
        float: The difference between the molecules, by some distance metric.
    """
    return num_atom_difference(source, target) # TODO


def reaction_possible(reaction, produced, timeline):
    """Calculate the earliest time a reaction is favorable.

    Parameters:
        reaction (Reaction): The reaction to consider.
        produced (Mapping[Product, Tuple[Reaction, time]]):
            When different chemicals have been produced
        timeline (Any): The temperature and pressure timeline.

    Returns:
        int: The time at which the reaction is favorable, or -1 otherwise.
    """
    reactants_ready = max(produced[reactant][1] for reactant in reaction.reactants)
    return reactants_ready + 1 # TODO


def search(reactions, initial_reactants, final_product):
    """Greedy hill-climbing to synthesize the product.

    Parameters:
        reactions (Iterable[Reaction]): List of reactions to consider.
        initial_reactants (List[Molecule]): List of initial reactants.
        final_product (Molecule): The product to synthesize.

    Raises:
        Exception: If the final product could not be synthesized.
    """

    # data structures
    inputs = defaultdict(set) # type: Dict[Reactant, Set[Reaction]]
    reactants = defaultdict(set) # type: Dict[Reaction, Set[Reactant]]
    outputs = defaultdict(set) # type: Dict[Product, Set[Reaction]]

    # variables
    queue = [] # type: List[Tuple[Tuple[float, int, str], Reaction]]
    produced = {} # type: Dict[Product, Tuple[Reaction, int]]

    def produce(product, time, producer=None):
        produced[product] = (producer, time)
        for reaction in inputs[product]:
            reactants[reaction].remove(product)
            if not reactants[reaction]:
                # priority is made of three things, compared in order:
                # 1. the difference from the goal molecule, a heuristics of sorts
                # 2. the earliest time a reaction is possible
                # 3. the string representation of the reaction, as a tie-breaker
                priority = (
                    min(
                        molecular_difference(product, final_product)
                        for product in reaction.products
                    ),
                    max(produced[reactant][1] for reactant in reaction.reactants),
                    str(reaction),
                )
                heappush(queue, (priority, reaction))

    # organize the reactions into data structures
    for reaction in reactions:
        for reactant in reaction.reactants:
            inputs[reactant].add(reaction)
        reactants[reaction].update(reaction.reactants)
        for product in reaction.products:
            outputs[product].add(reaction)

    # initialize the variables
    for reactant in initial_reactants:
        produce(reactant, 0)

    # hill climb
    while queue and final_product not in produced:
        _, reaction = heappop(queue)
        earliest_time = reaction_possible(reaction, produced, None)
        if earliest_time == -1:
            continue
        for product in reaction.products:
            if product not in produced:
                produce(product, earliest_time, producer=reaction)

    # error if search failed
    if not queue:
        raise Exception(
            f'unable to synthesize {final_product} from: '
            + f'{", ".join(str(x) for x in initial_reactants)}'
        )

    return produced


def print_search_results(initial_reactants, final_product, produced):
    # re-trace synthesis steps
    priorities = {
        None: (0, -len(produced)),
    }
    steps = set()
    product_queue = [(final_product, 0)]
    while product_queue:
        product, distance = product_queue.pop()
        (reaction, time) = produced[product]
        steps.add((reaction, product))
        if reaction is not None:
            priorities[reaction] = (time, -distance)
            for reactant in reaction.reactants:
                product_queue.append((reactant, distance + 1))

    # print synthesis steps
    print(
        f'synthesizing {final_product} from: '
        + f'{", ".join(str(x) for x in initial_reactants)}'
    )
    print()
    for reaction, product in sorted(steps, key=(lambda step: priorities[step[0]])):
        time, _ = priorities[reaction]
        if reaction is None:
            print(f'{time}: {product} given')
        else:
            print(f'{time}: {product} produced by {reaction}')


def visualize_reactions(reactions):
    lines = []
    lines.append('digraph {')
    chemicals = set()
    for reaction in reactions:
        lines.append(f'    "{reaction}" [shape="box"]')
        for reactant in reaction.reactants:
            lines.append(f'    "{reactant}" -> "{reaction}"')
        for product in reaction.products:
            lines.append(f'    "{reaction}" -> "{product}"')
    lines.append('}')
    print('\n'.join(lines))


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
    arg_parser = ArgumentParser()
    arg_parser.add_argument(dest='action', nargs='?', choices=['demo', 'search', 'visualize'], default='search', help='The action to perform')
    arg_parser.add_argument('-i', '--input', action='append', help='An initial reactant')
    arg_parser.add_argument('-o', '--output', action='store', help='The final product')
    arg_parser.add_argument('--reactions', action='store', default='LARGE_REACTION_SET', help='The reactions to allow')
    args = arg_parser.parse_args()
    if args.action == 'demo':
        args.action = 'search'
        args.input = ['H2', 'CO']
        args.output = 'HCOOH'
        args.reactions = 'LARGE_REACTION_SET'

    if not args.input:
        arg_parser.error('no initial reactant given')
    if not args.output:
        arg_parser.error('no final product given')

    initial_reactants = args.input
    final_product = args.output
    if args.reactions == 'LARGE_REACTION_SET':
        reaction_set = LARGE_REACTION_SET
    elif args.reactions == 'SMALL_REACTION_SET':
        reaction_set = SMALL_REACTION_SET
    else:
        arg_parser.error(f'undefined reactions set: {args.reactions}')

    parser = ReactionWalker()
    reactions = [
        parser.parse(reaction)
        for reaction in reaction_set
    ]
    initial_reactants = [
        parser.parse(reactant, 'molecule')
        for reactant in initial_reactants
    ]
    final_product = parser.parse(final_product, 'molecule')

    if args.action == 'search':
        produced = search(reactions, initial_reactants, final_product)
        print_search_results(initial_reactants, final_product, produced)
    elif args.action == 'visualize':
        visualize_reactions(reactions)
    else:
        arg_parser.error(f'undefined action: {args.action}')


if __name__ == '__main__':
    main()
