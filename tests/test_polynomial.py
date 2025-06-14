"""Test file for polynomial's classe"""


import unittest

from doctest import testmod

from reedsolomon import Polynomial


class TestPolynomial(unittest.TestCase):
    """Tests to check all of Polynomial's functionality"""

    # Docstring test
    def test_docstring(self):
        result = testmod(__import__("reedsolomon.polynomial"), verbose=True)

        self.assertEqual(result.failed, 0, f"Doctests failed: {result.failed} out of {result.attempted}")

    # String test
    def test_str_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(str(f), "f(x) = -3x² + 23x - 67")

    def test_str_2(self):
        g = Polynomial(x15=-1, x56=14, x0=35, name="g")

        self.assertEqual(str(g), "g(x) = 14x^56 - x^15 + 35")

    def test_str_3(self):
        coefficients = (-1, -1, 0)
        h = Polynomial(coefficients, name="h")

        self.assertEqual(str(h), "h(x) = -x² - x")

    # Representation test
    def test_repr_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(repr(f), "Polynomial(x2=-3, x=23, -67)")

    def test_repr_2(self):
        g = Polynomial(x15=-1, x56=14, x0=35, name="g")

        self.assertEqual(repr(g), "Polynomial(x56=14, x15=-1, 35)")

    def test_repr_3(self):
        coefficients = (-1, -1, 0)
        h = Polynomial(coefficients, name="h")

        self.assertEqual(repr(h), "Polynomial(x2=-1, x=-1)")

    # Attribute test
    def test_attributes_alphabetic_degree_1(self):
        coefficients = {'x5':-4, 'x4':3, 'x3':-56, 'x2':1, 'x1':26, 'x0':3}
        f = Polynomial(**coefficients)

        self.assertTupleEqual((f.a, f.b, f.c, f.d, f.e, f.f), (-4, 3, -56, 1, 26, 3))

    def test_attributes_alphabetic_degree_2(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)

        self.assertTupleEqual((f.a, f.b, f.c), coefficients)

    def test_attributes_alphabetic_degree_3(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)
        f.a = 34

        self.assertEqual(repr(f), "Polynomial(x2=34, x=23, -67)")

    def test_attributes_alphabetic_degree_4(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)
        del f.c

        self.assertEqual(repr(f), "Polynomial(x2=-3, x=23)")

    def test_attributes_x_degree_1(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)

        self.assertTupleEqual((f.x2, f.x1, f.x0), coefficients)

    def test_attributes_x_degree_2(self):
        coefficients = (45, -67, 89, 18, 0)
        f = Polynomial(coefficients)

        self.assertTupleEqual((f.x5, f.x4, f.x3, f.x2, f.x1, f.x0), (0, 45, -67, 89, 18, 0))

    def test_attributes_x_degree_3(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)
        f.x2 = 34

        self.assertEqual(repr(f), "Polynomial(x2=34, x=23, -67)")

    def test_attributes_x_degree_4(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)
        del f.x0

        self.assertEqual(repr(f), "Polynomial(x2=-3, x=23)")

    def test_attributes_index_degree_1(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)

        self.assertTupleEqual((f[2], f[1], f[0]), coefficients)

    def test_attributes_index_degree_2(self):
        coefficients = (45, -67, 89, 18, 0)
        f = Polynomial(coefficients)

        self.assertTupleEqual((f[5], f[4], f[3], f[2], f[1], f[0]), (0, 45, -67, 89, 18, 0))

    def test_attributes_index_degree_3(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)
        f[2] = 34

        self.assertEqual(repr(f), "Polynomial(x2=34, x=23, -67)")

    def test_attributes_index_degree_4(self):
        coefficients = (-3, 23, -67)
        f = Polynomial(coefficients)
        del f[0]

        self.assertEqual(repr(f), "Polynomial(x2=-3, x=23)")

    def test_tips_zero_1(self):
        f = Polynomial(0, 0, -3, 23, -67, 0)

        self.assertTupleEqual((f.coefficients, f.sparse), ([-3, 23, -67, 0], {'x3': -3, 'x2': 23, 'x1': -67, 'x0': 0}))

    def test_tips_iter_1(self):
        f = Polynomial(-3, 23, -67)
        x2, x1, x0 = f

        self.assertTupleEqual((x2, x1, x0), (-3, 23, -67))

    def test_tips_len_1(self):
        f = Polynomial(0, 0, 0, -3, 23, -67, 0)

        self.assertEqual(len(f), 4)

    def test_tips_list_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(list(f), [-3, 23, -67])

    def test_tips_contains_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertIn('x2', f)

    def test_tips_contains_2(self):
        f = Polynomial(-3, 23, -67)

        self.assertNotIn('x36', f)

    def test_tips_call_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(f(5), -27)

    def test_tips_call_2(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(f(0), -67)

    def test_tips_list_1(self):
        f = Polynomial(0, -3, 23, -67, 0)

        self.assertListEqual(list(f), [-3, 23, -67, 0])

    def test_tips_reversed_1(self):
        f = Polynomial(0, -3, 23, -67, 0)

        self.assertListEqual(list(reversed(f)), [-67, 23, -3])

    def test_tips_dict_1(self):
        f = Polynomial(0, -3, 23, -67, 0)

        self.assertDictEqual(f.sparse, {"x3":-3, "x2":23, "x1":-67, "x0":0})

    def test_degree_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(f.degree, 2)

    def test_degree_2(self):
        f = Polynomial(x56=0, x35=-1)

        self.assertEqual(f.degree, 35)

    def test_items_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(f.items(), [('x2', -3), ('x1', 23), ('x0', -67)])
        self.assertEqual(dict(f.items()), f.sparse)

    def test_coefficients_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertListEqual(f.coefficients, [-3, 23, -67])

        f.coefficients = (3, 56, -1, 89)

        self.assertListEqual(f.coefficients, [3, 56, -1, 89])

    def test_sparce_1(self):
        f = Polynomial(x3=3, x2=56, x1=-1, x0=89)

        self.assertDictEqual(f.sparse, {'x3': 3, 'x2': 56, 'x1': -1, 'x0': 89})

        f.sparse = {"x2":-3, "x1":23, "x0":-67}

        self.assertDictEqual(f.sparse, {'x2': -3, 'x1': 23, 'x0': -67})

    def test_derive_1(self):
        f = Polynomial(-3, 23, -67)

        self.assertEqual(f.derive(), Polynomial(x1=-134, x0=23))

    def test_copy_1(self):
        h = Polynomial(9, -30, 25, name="h")
        sh = str(h)
        
        _h = h.copy()
        _sh = str(_h)

        self.assertEqual(sh, _sh)

    def test_solve_equation_1(self):
        f = Polynomial(2, -26)

        self.assertEqual(f.solve(), 13.0)

    def test_solve_2_roots_equation_1(self):
        f = Polynomial(1, 1, -12)

        self.assertEqual(f.delta, 49)
        self.assertTupleEqual(f.solve(), (-4.0, 3.0))

    def test_solve_1_roots_equation_1(self):
        h = Polynomial(9, -30, 25, name="h")

        self.assertEqual(h.delta, 0)
        self.assertTupleEqual(h.solve(), (1.6666666666666667,))

    def test_solve_equation_2(self):
        f = Polynomial(x5=3, x2=23, x0=14)

        self.assertTupleEqual(f.solve(), ((-2.061777647931321+0j),
                                          (-0.023547692027920532-0.7769582400529271j),
                                          (-0.023547692027920532+0.7769582400529271j),
                                          (1.054436515993581-1.6230188809835164j),
                                          (1.054436515993581+1.6230188809835164j)))

    def test_expression_developed_1(self):
        f = Polynomial(1, 1, -12)

        self.assertEqual(f.developed(), "x² + x - 12")

    def test_expression_canonic_1(self):
        f = Polynomial(1, 1, -12)

        self.assertTupleEqual((f.alpha, f.beta), (-0.5, -12.25))
        self.assertEqual(f.canonic(), "(x + 0.5)² - 12.25")

    def test_expression_factorise_1(self):
        f = Polynomial(1, 1, -12)

        self.assertEqual(f.factorised(), "(x + 4.0)(x - 3.0)")

    def test_expression_factorise_2(self):
        h = Polynomial(9, -30, 25, name="h")

        self.assertEqual(h.factorised(), "9(x - 1.667)²")

    def test_manipulation_equal_1(self):
        f = Polynomial(1, 1, -12)
        h = Polynomial(9, -30, 25, name="h")

        self.assertNotEqual(f, h)

    def test_manipulation_equal_2(self):
        f = Polynomial(9, -30, 25)
        h = Polynomial(9, -30, 25, name="h")

        self.assertEqual(f, h)

    def test_manipulation_negative_1(self):
        f = Polynomial(9, -30, 25)
        h = -f
        h.name = "h"

        self.assertEqual(str(h), "h(x) = -9x² + 30x - 25")

    def test_add_1(self):
        one = Polynomial([2,4,7,3])
        two = Polynomial([5,2,4,2])

        r = one + two

        self.assertEqual(list(r.coefficients), [7, 6, 11, 5])

    def test_add_2(self):
        one = Polynomial([2,4,7,3,5,2])
        two = Polynomial([5,2,4,2])

        r = one + two

        self.assertEqual(list(r.coefficients), [2,4,12,5,9,4])

    def test_add_3(self):
        one = Polynomial([7,3,5,2])
        two = Polynomial([6,8,5,2,4,2])

        r = one + two

        self.assertEqual(list(r.coefficients), [6,8,12,5,9,4])

    def test_add_4(self):
        f = Polynomial(x5=35, x2=-3, x0=1)
        g = Polynomial(21, 5, -1, name="g")

        h = f + g
        h.name = "h"

        self.assertEqual(str(h), "h(x) = 35x^5 + 18x² + 5x")

    def test_add_5(self):
        f = Polynomial(x5=35, x2=-3, x0=1)
        g = Polynomial(21, 5, -1, name="g")

        h = g + 1
        h.name = "h"

        self.assertEqual(str(h), "h(x) = 21x² + 5x")

    def test_sub_1(self):
        f = Polynomial(x5=35, x2=-3, x0=1)
        g = Polynomial(21, 5, -1, name="g")

        h = f - g
        h.name = "h"

        self.assertEqual(str(h), "h(x) = 35x^5 - 24x² - 5x + 2")

    def test_sub_2(self):
        f = Polynomial(x5=35, x2=-3, x0=1)
        g = Polynomial(21, 5, -1, name="g")

        h = f - 1
        h.name = "h"

        self.assertEqual(str(h), "h(x) = 35x^5 - 3x^2")

    def test_mul_1(self):
        f = Polynomial(x5=35, x2=-3, x0=1)
        g = Polynomial(21, 5, -1, name="g")

        h = f * 3
        h.name = "h"

        self.assertEqual(str(h), "h(x) = 105x^5 - 9x² + 3")

    def test_mul_2(self):
        f = Polynomial(x5=35, x2=-3, x0=1)
        g = Polynomial(21, 5, -1, name="g")

        h = f * g
        h.name = "h"

        self.assertEqual(str(h), "h(x) = 735x^7 + 175x^6 - 35x^5 - 63x^4 - 15x^3 + 24x² + 5x - 1")

    def test_mul_3(self):
        i = Polynomial(x5=1, x2=1, x0=1)
        j = Polynomial(x2=1, x1=1, x0=1)

        h = i * j
        h.name = "h"

        self.assertEqual(str(h), "h(x) = x^7 + x^6 + x^5 + x^4 + x^3 + 2x² + x + 1")

    def test_div_1(self):
        f = Polynomial(x5=35, x2=-3, x0=1)
        self._extracted_from_test_div_1(f, "h(x) = 17.5x^5 - 1.5x² + 0.5")

        f = Polynomial(2, 3, -1, 5)
        self._extracted_from_test_div_1(f, "h(x) = x^3 + 1.5x² - 0.5x + 2.5")

    def _extracted_from_test_div_1(self, f, arg1):
        result = f / 2

        return self._extracted_from_rename(result, arg1)

    def test_div_2(self):
        f = Polynomial(2, 3, -1, 5)
        g = Polynomial(1, 0, 1, name="g")

        h = f / g
        h.name = "h"

        self.assertEqual(str(h), "h(x) = 2x + 3")

    def test_div_3(self):
        one = Polynomial([1,4,0,3])
        two = Polynomial([1,0,1])

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, [1,4])
        self.assertEqual(r.coefficients, [-1,-1])

    def test_div_4(self):
        one = Polynomial([1,0,0,2,2,0,1,2,1])
        two = Polynomial([1,0,-1])

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, [1,0,1,2,3,2,4])
        self.assertEqual(r.coefficients, [4,5])

    def test_div_5(self):
        # 0 quotient
        one = Polynomial([1,0,-1])
        two = Polynomial([1,1,0,0,-1])

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, [0])
        self.assertEqual(r.coefficients, [1,0,-1])

    def test_div_6(self):
        one = Polynomial([1,0,0,2,2,0,1,-2,-4])
        two = Polynomial([1,0,-1])

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, [1,0,1,2,3,2,4])
        self.assertEqual(r.coefficients, [0])

    def test_floor_division_1(self):
        f = Polynomial(2, 3, -1, 5)
        self._extracted_from_test_floor_division_1(f, 2, "h(x) = x^3 + 1.5x² - 0.5x + 2.5")

        g = Polynomial(1, 0, 1, name="g")
        self._extracted_from_test_floor_division_1(f, g, "h(x) = 2x + 3")

    def _extracted_from_test_floor_division_1(self, f, arg1, arg2):
        result = f // arg1
        return self._extracted_from_rename(result, arg2)

    def test_modulo_1(self):
        f = Polynomial(2, 3, -1, 5)
        self._extracted_from_test_modulo_1(f, 2, "h(x) = 0")

        g = Polynomial(1, 0, 1, name="g")
        self._extracted_from_test_modulo_1(f, g, "h(x) = -3x + 2")

    def _extracted_from_test_modulo_1(self, f, arg1, arg2):
        result = f % arg1
        return self._extracted_from_rename(result, arg2)

    def test_divmod_1(self):
        f = Polynomial(2, 3, -1, 5)
        g = Polynomial(1, 0, 1, name="g")

        quotient, remainder = divmod(f, g)

        self.assertEqual(quotient, Polynomial(x1=2.0, x0=3.0))
        self.assertEqual(remainder, Polynomial(x1=-3.0, x0=2.0))

        h = Polynomial(2, 3) * Polynomial(1, 0, 1) + Polynomial(-3, 2)
        h.name = "h"

        self.assertEqual(f, h)

    def test_pow_1(self):
        f = Polynomial(5, -3, 1)

        h = self._extracted_from_test_pow_1(f, 2, "h(x) = 25x^4 - 30x^3 + 19x² - 6x + 1")
        h = self._extracted_from_test_pow_1(f, 3, "h(x) = 125x^6 - 225x^5 + 210x^4 - 117x^3 + 42x² - 9x + 1")

        f **= 3

        self.assertEqual(f, h)

    def _extracted_from_test_pow_1(self, f, arg1, arg2):
        result = f**arg1
        return self._extracted_from_rename(result, arg2)

    # TODO Rename this here and in `_extracted_from_test_div_1`, `_extracted_from_test_floor_division_1`, `_extracted_from_test_modulo_1` and `_extracted_from_test_pow_1`
    def _extracted_from_rename(self, result, arg1):
        result.name = "h"
        self.assertEqual(str(result), arg1)
        return result

    def test_bool_1(self):
        f = Polynomial(5, -3, 1)

        self.assertTrue(f)

        h = Polynomial()

        self.assertFalse(h)

    def test_int_1(self):
        f = Polynomial(5, -3, 1)

        self.assertEqual(int(f), 2)
        self.assertEqual(f.degree, 2)

if __name__ == '__main__':
    unittest.main()