#!/usr/bin/env python3
"""A script to determine reactions necessary to synthesize a molecule."""

from collections import namedtuple, Counter
from decimal import Decimal
from itertools import product as cross_product
from os.path import realpath, join as join_path, dirname
from textwrap import dedent

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


if __name__ == '__main__':
    main()
