import unittest

from hand_sides import HandSide


class TestHandSide(unittest.TestCase):
    def testFloatCoefficients(self):
        formula = '2.54A + 0.4B + C + 5D'
        hand_side = HandSide(formula)
        self.assertEqual(2.54, hand_side.stoichiometry['A'])
        self.assertEqual(0.4, hand_side.stoichiometry['B'])
        self.assertEqual(1, hand_side.stoichiometry['C'])
        self.assertEqual(5, hand_side.stoichiometry['D'])