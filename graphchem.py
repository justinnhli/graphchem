#!/usr/bin/env python3
"""A script to determine reactions necessary to synthesize a molecule."""

from argparse import ArgumentParser
from collections import namedtuple, Counter, defaultdict
from decimal import Decimal
from heapq import heappush, heappop
from os.path import realpath, join as join_path, dirname

from typing import Any, Optional, Generator, Iterable, Sequence, Mapping, Tuple, List, Set, Dict

from pegparse import ASTWalker, create_parser_from_file


MoleculeCount = namedtuple('MoleculeCount', 'count molecule')


def least_common_mulitple(*ints):
    # type: (int) -> int
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

    # pylint: disable = unused-argument, no-self-use

    def __init__(self):
        # type: (ReactionWalker) -> None
        """Initialize a ReactionWalker."""
        super().__init__(
            create_parser_from_file(join_path(
                dirname(realpath(__file__)),
                'reactions.peg',
            )),
            'reaction',
        )

    def _parse_reaction(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> Reaction
        return Reaction(results[0], results[1])

    def _parse_molecule_list(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> List[Any]
        return results

    def _parse_molecule_count(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> MoleculeCount
        assert len(results) in (1, 2)
        if len(results) == 1:
            return MoleculeCount(Decimal(1), results[0])
        else:
            return MoleculeCount(*results)

    def _parse_molecule(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> Molecule
        return Molecule(*results, formula=ast.match)

    def _parse_group_count(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> Dict[str, int]
        if len(results) == 1:
            return results[0]
        elif len(results) == 2:
            result = results[0]
            for key, count in result.items():
                result[key] = count * results[1]
            return result
        else:
            assert False
            return None

    def _parse_group(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> List[Any]
        return results

    def _parse_element(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> Dict[str, int]
        return Counter({ast.match})

    def _parse_number(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> Decimal
        return Decimal(ast.match)

    def _parse_int(self, ast, results):
        # type: (ReactionWalker, ASTNode, List[Any]) -> int
        return int(ast.match)


class Molecule:
    """A chemistry molecule."""

    def __init__(self, *components, formula="", name=""):
        # type: (Mapping[str, int], str, str) -> None
        """Initialize the Molecule.

        Parameters:
            *components (Mapping[str, int]): The chemical formula components
                (element and group counts).
            formula (str): The formula for the molecule. Optional.
            name (str): The name of the molecule, if any. Optional.
        """
        self.name = name
        if formula:
            self.formula = formula
        else:
            self.formula = ''.join(
                f'{element}{count}' if count > 1 else element
                for component in components
                for element, count in component.items()
            )
        self.atoms = Counter() # type: Dict[str, int]
        for component in components:
            for element, count in component.items():
                self.atoms[element] += count

    def __eq__(self, other):
        # type: (Molecule, Any) -> bool
        return isinstance(other, Molecule) and str(self) == str(other)

    def __hash__(self):
        # type: (Molecule) -> int
        return hash(str(self))

    def __str__(self):
        # type: (Molecule) -> str
        return self.formula


Reactant = Molecule
Product = Molecule


class Reaction:
    """A chemical reaction."""

    def __init__(self, reactants, products, energy=0):
        # type: (Reaction, Sequence[MoleculeCount], Sequence[MoleculeCount], float) -> None
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
        # type: (Reaction) -> List[Molecule]
        """Get the reactants.

        Returns:
            List[Molecule]: The reactants of this reaction.
        """
        return [reactant.molecule for reactant in self.reactant_counts]

    @property
    def products(self):
        # type: (Reaction) -> List[Molecule]
        """Get the products.

        Returns:
            List[Molecule]: The products of this reaction.
        """
        return [product.molecule for product in self.product_counts]

    def _integerize(self):
        # type: (Reaction) -> None
        """Convert reaction equation to integers.

        This function does not reduce the coefficients to the smallest values.
        """
        multiplier = least_common_mulitple(
            *(group.count.as_integer_ratio()[1] for group in self.reactant_counts),
            *(group.count.as_integer_ratio()[1] for group in self.product_counts),
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
        # type: (Reaction) -> None
        """Check that atoms are conserved in this reaction."""
        reactants_counter = Counter() # type: Dict[Molecule, int]
        for count, molecule in self.reactant_counts:
            for _ in range(count):
                reactants_counter.update(molecule.atoms)
        products_counter = Counter() # type: Dict[Molecule, int]
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
        # type: (Reaction, Any) -> bool
        return isinstance(other, Reaction) and str(self) == str(other)

    def __hash__(self):
        # type: (Reaction) -> int
        return hash(str(self))

    def __str__(self):
        # type: (Reaction) -> str
        return ' '.join([
            ' + '.join(
                ('' if count == 1 else f'{count} ') + str(molecule)
                for count, molecule in self.reactant_counts
            ),
            '=',
            ' + '.join(
                ('' if count == 1 else f'{count} ') + str(molecule)
                for count, molecule in self.product_counts
            ),
        ])


TempPres = namedtuple('TempPres', 'temperature, pressure') # in Celsius and kilopascals
Priority = namedtuple('Priority', 'heuristic, time, distance, string')
ProductionMetadata = namedtuple('ProductionMetadata', 'reaction, time, distance')


Molecules = Iterable[Molecule] # pylint: disable = unused-variable
Reactions = Iterable[Reaction] # pylint: disable = unused-variable
Timeline = Sequence[TempPres] # pylint: disable = unused-variable
SearchResult = Dict[Product, ProductionMetadata] # pylint: disable = unused-variable


def num_atom_difference(source, target):
    # type: (Molecule, Molecule) -> int
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
    # type: (Molecule, Molecule) -> float
    """Calculate the "difference" between two molecules.

    Parameters:
        source (Molecule): The source molecule.
        target (Molecule): The target molecule.

    Returns:
        float: The difference between the molecules, by some distance metric.
    """
    return num_atom_difference(source, target) # TODO


def reaction_possible(reaction, temp_pres):
    # type: (Reaction, TempPres) -> bool
    """Determine whether a reaction is possible at a temperature and pressure.

    Parameters:
        reaction (Reaction): The reaction to consider.
        temp_pres (TempPres): The temperature and pressure.

    Returns:
        bool: If the reaction is possible.
    """
    # pylint: disable = unused-argument
    return True # TODO


def reaction_first_possible(reaction, produced, timeline):
    # type: (Reaction, SearchResult, Timeline) -> int
    """Calculate the earliest time a reaction is possible.

    This function assumes all reactants are already produced.

    Parameters:
        reaction (Reaction): The reaction to consider.
        produced (SearchResult): When different chemicals have been produced
        timeline (Timeline): The temperature and pressure timeline.

    Returns:
        int: The time at which the reaction is possible, or -1 otherwise.
    """
    reactants_ready = max(produced[reactant][1] for reactant in reaction.reactants)
    for time, temp_pres in enumerate(timeline[reactants_ready:], start=reactants_ready):
        if reaction_possible(reaction, temp_pres):
            return time
    return -1


def search(reactions, initial_reactants, final_product, timeline):
    # type: (Reactions, Molecules, Molecule, Timeline) -> Generator[SearchResult, None, None]
    """Greedily hill-climbing to synthesize the product.

    This search algorithm is slightly weird because the optimization metric is
    not the same as the heuristic/fitness function. We would like to know the
    earliest time at which the final product could be synthesized, which we can
    determine via a topological traversal of the reaction network (with
    additional timeline constraints). This is computationally expensive,
    however, as it would require running all possible reactions at every
    timestep.

    Instead, we base our algorithm on hill-climbing, so that reactions with
    products that resemble the final product occur first. The disadvantage of
    this approach is that the final product will be synthesized at a later time
    than the minimum. To remedy this disadvantage, the algorithm is coded to
    have the anytime property: after it finds one pathway for synthesizing the
    final product, it will continue searching for pathways that take less time
    to do so. This is why the function is a generator: each successive result
    will be a synthesis pathway that completes earlier than the previous result.
    In the limit, this algorithm will find the same pathway (or one that
    completes in the same amount of time) as the topological approach, although
    it sacrifices both runtime (due to hill-climbing exploring later reactions)
    and memory (to allow for the anytime property) to do so.

    Algorithmically, we use a priority queue to keep track of which reactions
    do not have any missing reactants, sorted by four metrics:

    1. the difference from the goal molecule, a pseudo-heuristic
    2. a pessimistic estimate of the earliest time all reactants are synthesized
    3. the shortest distance to an initial reactant
    4. the string representation of the reaction, as a tie-breaker

    Each time we pop a reaction from the queue, its products are considered
    synthesized, and we record the time that the reaction occurred. To achieve
    the anytime property, products can be "re-synthesized" if the algorithm
    encounters a different reaction that can occur earlier. The reactions that
    consume that product are then added back into the queue, with a new estimate
    of when they could occur, thus propagating this earlier synthesis pathway.

    This narrative omits an optimization: when a reaction is popped from the
    queue, it is not always necessary to consider its effects:

    1. If the reaction is never possible given the temperature and pressure
    conditions, the reaction could obviously never occur and can be ignored.

    2. If the reaction has already occurred at an earlier time. This is possible
    if a reaction is added to the priority queue multiple times, due to a
    reactant having an improved synthesis time estimate. The later reaction will
    not enable any new products, nor lead to an earlier synthesis time, and thus
    can be ignored.

    3. If the reaction cannot occur until after the final product is
    synthesized. For a similar reason, this reaction could not lead to an
    earlier synthesis pathway and can be ignored.

    Similarly, when a molecule is synthesized, not all consuming reactions need
    to be added to the priority queue:

    1. If the reaction has other reactants that have not yet been synthesized,
    the reaction obviously cannot occur at this point.

    2. If the reaction has already occurred at an earlier time compared to when
    the molecule was produced, there is no need to add the reaction to the
    queue. This will be caught by the second case above regardless, but it
    doesn't hurt to reduce the queue size.

    3. If new production time does not allow the reaction could occur any
    earlier, and so do not have to be added to the queue.

    Parameters:
        reactions (Reactions): List of reactions to consider.
        initial_reactants (Molecules): List of initial reactants.
        final_product (Molecule): The product to synthesize.
        timeline (Timeline): The temperature and pressure timeline.

    Yields:
        SearchResult: The produced molecules and the reaction
            and time they were produced.

    Raises:
        Exception: If the final product could not be synthesized.
    """

    # data structures
    consumed_by = defaultdict(set) # type: Dict[Reactant, Set[Reaction]]
    produced_by = defaultdict(set) # type: Dict[Product, Set[Reaction]]

    # search variables
    missing_reactants = defaultdict(set) # type: Dict[Reaction, Set[Reactant]]
    queue = [] # type: List[Tuple[Priority, Reaction]]
    produced = {} # type: SearchResult
    reacted = {} # type: Dict[Reaction, int]

    # organize the reactions into data structures
    for reaction in reactions:
        for reactant in reaction.reactants:
            consumed_by[reactant].add(reaction)
        missing_reactants[reaction].update(reaction.reactants)
        for product in reaction.products:
            produced_by[product].add(reaction)

    def produce(product, time, distance, producer=None):
        # type: (Product, int, int, Optional[Reaction]) -> None
        """Add a molecule to the production record.

        Parameters:
            product (Product): The molecule that was produced.
            time (int): The time at which the molecule was produced.
            distance (int): The minimum distance to an initial reactant.
            producer (Optional[Reaction]): The reaction that produced the molecule.
        """
        prev_time = -1
        if product in produced:
            prev_time = produced[product].time
        produced[product] = ProductionMetadata(producer, time, distance)
        for reaction in consumed_by[product]:
            missing_reactants[reaction].discard(product)
            if missing_reactants[reaction]:
                continue
            if reacted.get(reaction, time + 1) <= time:
                continue
            reactants_ready = max(produced[reactant].time for reactant in reaction.reactants)
            if 0 <= prev_time <= reactants_ready:
                continue
            priority = Priority(
                min(
                    molecular_difference(product, final_product)
                    for product in reaction.products
                ),
                reactants_ready,
                distance,
                str(reaction),
            )
            heappush(queue, (priority, reaction))

    # initialize the variables
    for reactant in initial_reactants:
        produce(reactant, time=0, distance=0)

    # hill climb
    while queue:
        priority, reaction = heappop(queue)
        earliest_time = reaction_first_possible(reaction, produced, timeline)
        if earliest_time == -1:
            continue
        if reacted.get(reaction, earliest_time + 1) <= earliest_time:
            continue
        if final_product in produced and produced[final_product].time < earliest_time:
            continue
        reacted[reaction] = earliest_time
        for product in reaction.products:
            if product not in produced or earliest_time < produced[product].time:
                produce(product, earliest_time, priority.distance + 1, producer=reaction)
                if product == final_product:
                    yield produced


def print_search_results(initial_reactants, final_product, timeline, produced):
    # type: (Reactions, Molecule, Timeline, SearchResult) -> None
    """Print search results.

    Parameters:
        initial_reactants (Molecules): List of initial reactants.
        final_product (Molecule): The product to synthesize.
        timeline (Timeline): The temperature and pressure timeline.
        produced (SearchResult): When different chemicals have been produced
    """
    # re-trace synthesis steps
    priorities = {
        None: (0, -len(produced)),
    } # type: Dict[Optional[Reaction], Tuple[int, float]]
    steps = set()
    product_queue = [final_product]
    while product_queue:
        product = product_queue.pop()
        (reaction, time, distance) = produced[product]
        steps.add((reaction, product))
        if reaction is not None:
            priorities[reaction] = (time, distance)
            for reactant in reaction.reactants:
                product_queue.append(reactant)

    # print synthesis steps
    print(
        f'synthesizing {final_product} from: '
        + f'{", ".join(str(x) for x in initial_reactants)}'
    )
    print()
    for reaction, product in sorted(steps, key=(lambda step: priorities[step[0]])):
        time, _ = priorities[reaction]
        temperature, pressure = timeline[time]
        if reaction is None:
            print(f'{time} ({temperature}C, {pressure}kPa): {product} given')
        else:
            print(f'{time} ({temperature}C, {pressure}kPa): {product} produced by {reaction}')


def visualize_reactions(reactions):
    # type: (Reactions) -> None
    lines = []
    lines.append('digraph {')
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
    # type: () -> None
    arg_parser = ArgumentParser()
    arg_parser.add_argument(
        dest='action', nargs='?',
        choices=['demo', 'search', 'visualize'], default='search',
        help='The action to perform',
    )
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
    timeline = [TempPres(0, 100), TempPres(100, 100)]

    if args.action == 'search':
        results = search(reactions, initial_reactants, final_product, timeline)
        for result in results:
            print()
            print_search_results(
                initial_reactants,
                final_product,
                timeline,
                result,
            )
    elif args.action == 'visualize':
        visualize_reactions(reactions)
    else:
        arg_parser.error(f'undefined action: {args.action}')


if __name__ == '__main__':
    main()
