import unittest

from ode_system import ODESystem
from reaction_network.abstract import ReactionNetwork


class TestODEToReactionNetwork(unittest.TestCase):
    REACTION_RATES = {'k': 1, 'm': 1, 'n': 1}

    def testLotka(self):
        system = ODESystem("s'=-ksx+ns\nx'=ksx-mx", ['s', 'x'])
        expected_reaction_network = ReactionNetwork("X->Q\nS->2S\nS+X->2X", external_reactants={"Q"}, rates=['n', 'm', 'k'])
        self.assertEqual(expected_reaction_network, system.to_reaction_network(self.REACTION_RATES))

    def testModifiedGompertz(self):
        system = ODESystem("s'=-ns\nx'=ksx+ns", ['s', 'x'])
        expected_reaction_network = ReactionNetwork("S->X\nS+X->2X+S", external_reactants={"Q"}, rates=['n', 'k'])
        self.assertEqual(expected_reaction_network, system.to_reaction_network(self.REACTION_RATES))

    def testGompertz(self):
        system = ODESystem("s'=-ns\nx'=ksx", ['s', 'x'])
        expected_reaction_network = ReactionNetwork("S->Q\nS+X->2X+S", external_reactants={"Q"}, rates=['n', 'k'])
        self.assertEqual(expected_reaction_network, system.to_reaction_network(self.REACTION_RATES))

    def testLogistic(self):
        system = ODESystem("s'=-ksx\nx'=ksx", ['s', 'x'])
        expected_reaction_network = ReactionNetwork("S+X->2X", external_reactants={"Q"}, rates=['k'])
        self.assertEqual(expected_reaction_network, system.to_reaction_network(self.REACTION_RATES))

    def testModifiedGompertzBateman(self):
        system = ODESystem("s'=-ks\nx'=npx\np'=ks-npx", ['s', 'x', 'p'])
        expected_reaction_network = ReactionNetwork("P+X->2X\n"
                                                    "S->P", external_reactants={"Q"}, rates=['k', 'n'])
        self.assertEqual(expected_reaction_network, system.to_reaction_network(self.REACTION_RATES))
