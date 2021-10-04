"""
Implements functons for printing various objects in different forms.
"""
import re


def string_reaction_network_to_ode_system(reaction_network, rate_names=None):
    if rate_names is None:
        if reaction_network.RATE_NAMES is None:
            rate_names = [f"k_{i}" for i in range(len(reaction_network.rates))]
        else:
            rate_names = reaction_network.RATE_NAMES
    equations = [r.lower() + "'=" for r in reaction_network.internal_reactants]
    for i, reaction in enumerate(reaction_network.reactions):
        term = rate_names[i]
        for reactant in reaction_network.internal_reactants:
            order = reaction.left.stoichiometry[reactant]
            if order == 0:
                continue
            term += reactant.lower()
            if order > 1:
                term += '^' + str(order)
        if term == rate_names[i]:
            raise ValueError(f"No reactants in term for reaction {reaction.__repr__()}")
        for j in range(len(equations)):
            multiplicity = reaction.stoichiometry[reaction_network.internal_reactants[j]]
            if multiplicity == 0:
                continue
            equations[j] += '+' if multiplicity > 0 else '-'
            if abs(multiplicity) != 1:
                equations[j] += str(abs(multiplicity))
            equations[j] += term
    output = '\n'.join(equations)
    output = re.sub(r'=\+', '=', output)
    return output


def latex_reaction_network_to_ode_system(reaction_network, rate_names=None):
    if rate_names is None:
        if reaction_network.RATE_NAMES is None:
            rate_names = [f"k_{i}" for i in range(len(reaction_network.rates))]
        else:
            rate_names = reaction_network.RATE_NAMES
    equations = ["\\dot " + r.lower() + "(t) &=" for r in reaction_network.internal_reactants]
    for i, reaction in enumerate(reaction_network.reactions):
        term = rate_names[i] + ' '
        for reactant in reaction_network.internal_reactants:
            order = reaction.left.stoichiometry[reactant]
            if order == 0:
                continue
            term += reactant.lower()
            if order > 1:
                term += '^{' + str(order) + '}'
        if term == rate_names[i]:
            raise ValueError(f"No reactants in term for reaction {reaction.__repr__()}")
        for j in range(len(equations)):
            multiplicity = reaction.stoichiometry[reaction_network.internal_reactants[j]]
            if multiplicity == 0:
                continue
            equations[j] += ' + ' if multiplicity > 0 else ' - '
            if abs(multiplicity) != 1:
                equations[j] += str(abs(multiplicity)) + ' '
            equations[j] += term
    output = r'\\ '.join(equations)
    output = re.sub(r'= \+', '=', output)
    return output


def latex_reaction_network(reaction_network):
    """
    Generates valid LaTeX for an align environment.
    """
    lines = []
    for reaction, rate_name in zip(reaction_network.reactions, reaction_network.rate_names):
        line = reaction.__repr__()
        line = line.replace('->', rf"&\xrightarrow{{{rate_name}}}")
        lines.append(line)
    return '\\\\\n'.join(lines)

