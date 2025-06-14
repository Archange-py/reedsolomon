#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file contains the Polynomial class, the basic object for use in coordination with Galois fields,
and together with the Reed-Solomon corrector code. This class was inspired by the attributes and
initialization of the Polynomial class in this project's polynomial.py file:
https://github.com/lrq3000/unireedsolomon.

In addition to implementing the algorithms required for the Reed-Solomon code, it has been given
additional functionality taken from first-grade courses, so don't worry if it seems both simple
and complicated !
"""

from typing import Iterable, Iterator, Self, Any

try:
    sp = __import__("sympy")

except ImportError:
    sp = None

class PolynomialError(Exception):
    """Base class for handling polynomial errors"""
    pass


class Polynomial:
    """Basic class for implementing a Polynomial in object-oriented programming
    
    Args:
        *coefficients (tuple, optional): An iterable containing the various coefficients of the polynomial, or a list of coefficients from largest to smallest.
        name (str, optional): function name. Defaults to "f".
        **sparse (dict, optional): An association of x with its associated coefficient.

    Examples:

        * There are three ways to initialize a Polynomial object.

        Simply, creates a polynomial using the arguments
        as coefficients in order of decreasing power:
        >>> print(Polynomial(5, 4, 3, 2, 1, 0))
        f(x) = 5x^5 + 4x^4 + 3x^3 + 2x² + x

        With a list, tuple, or other iterable, creates a polynomial using
        the items as coefficients in order of decreasing power:
        >>> print(Polynomial([5, 0, 0, 0, 0, 0]))
        f(x) = 5x^5

        With keyword arguments such as for example x3=5, sets the
        coefficient of x^3 to be 5:
        >>> print(Polynomial(x32=5, x64=8))
        f(x) = 8x^64 + 5x^32

        >>> print(Polynomial(x5=5, x9=4, x0=2)) 
        f(x) = 4x^9 + 5x^5 + 2

        With no arguments, creates an empty polynomial:
        >>> Polynomial() == Polynomial(0) == Polynomial([0, 0, 0])
        True

        * You can get attributes using multiple syntaxes:

        X notation:

        >>> f = Polynomial(-3, 23, -67)
        >>> f.x2, f.x1, f.x0
        (-3, 23, -67)
        >>> f.x2 = 34
        >>> f
        Polynomial(x2=34, x=23, -67)

        Alphabetic notation (max is z):

        >>> f = Polynomial(-3, 23, -67)
        >>> f.a, f.b, f.c
        (-3, 23, -67)
        >>> f.a = -3
        >>> f
        Polynomial(x2=-3, x=23, -67)
        >>> f = Polynomial(x5=-4, x4=3, x3=-56, x2=1, x1=26, x0=3)
        >>> f.a, f.b, f.c, f.d, f.e, f.f
        (-4, 3, -56, 1, 26, 3)

        List notation:

        >>> f = Polynomial(-3, 23, -67)
        >>> f[2], f[1], f[0]
        (-3, 23, -67)
        >>> f[2] = 34
        >>> f
        Polynomial(x2=34, x=23, -67)

        * Operations:

        Addition:

        >>> f = Polynomial(x5=35, x2=-3, x0=1)
        >>> g = Polynomial(21, 5, -1, name="g")

        >>> h = f + g
        >>> h.name = "h"
        >>> print(h)
        h(x) = 35x^5 + 18x² + 5x

        >>> h = g + 1
        >>> h.name = "h"
        >>> print(h)
        h(x) = 5x + 1

        Subtraction:

        >>> f = Polynomial(x5=35, x2=-3, x0=1)
        >>> g = Polynomial(21, 5, -1, name="g")

        >>> h = f - g
        >>> h.name = "h"
        >>> print(h)
        h(x) = 35x^5 - 24x² - 5x + 2

        >>> h = f - 1
        >>> h.name = "h"
        >>> print(h)
        h(x) = 35x^5 - 3x^2

        Multiplication:

        >>> h = f * g
        >>> h.name = "h"
        >>> print(h)
        h(x) = 175x^6 - 15x^3 + 5x

        >>> h = f * 3
        >>> h.name = "h"
        >>> print(h)
        h(x) = 105x^5 - 9x² + 3

        Negative:

        >>> f = Polynomial(5, -3, 1)
        >>> h = -f
        >>> h.name = "h"
        >>> print(h)
        h(x) = -5x² + 3x - 1

        Power:

        >>> f = Polynomial(5, -3, 1)
        >>> print(f)
        f(x) = 5x² - 3x + 1
        >>> print(f * f)
        f(x) = 25x^4 - 30x^3 + 19x² - 6x + 1

        >>> h = f**2 # only positive and integers power
        >>> h.name = "h"
        >>> print(h)
        h(x) = 25x^4 - 30x^3 + 19x² - 6x + 1

        Division:

        >>> f = Polynomial(2, 3, -1, 5)
        >>> g = Polynomial(1, 0, 1, name="g")
        >>> print(f)
        f(x) = 2x^3 + 3x² - x + 5
        >>> print(g)
        g(x) = x² + 1

        >>> h = f / 2
        >>> h.name = "h"
        >>> print(h)
        h(x) = x^3 + 1.5x² - 0.5x + 2.5

        >>> h = f / g # # Warning: h is just the quotient, a remainder can exist !
        >>> h.name = "h"
        >>> print(h)
        h(x) = 2x + 3

        Floor Division:

        >>> h = f // 2
        >>> h.name = "h"
        >>> print(h)
        h(x) = x^3 + 1.5x² - 0.5x + 2.5

        >>> h = f // g # it's the same of f / g
        >>> h.name = "h"
        >>> print(h)
        h(x) = 2x + 3

        Modulo:

        >>> h = f % 2
        >>> h.name = "h"
        >>> print(h)
        h(x) = 0

        >>> h = f % g
        >>> h.name = "h"
        >>> print(h)
        h(x) = -3x + 2

        Divmod:

        >>> divmod(f, g)
        (Polynomial(x=2, 3), Polynomial(x=-3, 2))

        >>> h = Polynomial(2, 3) * Polynomial(1, 0, 1) + Polynomial(-3, 2)
        >>> h.name = "h"
        >>> print(h)
        h(x) = 2x^3 + 3x² - x + 5

        >>> h == f
        True

        * Tips:

        >>> f = Polynomial(-3, 23, -67)
        >>> x2, x1, x0 = f
        >>> x2, x1, x0
        (-3, 23, -67)

        >>> len(f)
        3

        >>> list(f)
        [-3, 23, -67]
        >>> list(reversed(f))
        [-67, 23, -3]

        >>> 'x2' in f
        True
        >>> 'x36' in f
        False

        >>> f(5)
        -27
        >>> f(0)
        -67

        >>> int(f) == f.degree
        True

        >>> bool(f)
        True
        >>> bool(Polynomial())
        False

        >>> f = Polynomial(1, 1, -12)
        >>> h = Polynomial(9, -30, 25, name="h")
        >>> f == h
        False
        >>> f != h
        True
        >>> f = Polynomial(9, -30, 25)
        >>> f == h
        True

        >>> f = Polynomial(-3, 23, -67)
        >>> [coef for coef in f]
        [-3, 23, -67]
        >>> [coef for coef in reversed(f)]
        [-67, 23, -3]
        >>> [(x, v) for x, v in f.items()]
        [('x2', -3), ('x1', 23), ('x0', -67)]

    Raises:
        PolynomialError: Specify coefficients list or keyword terms, not both.
        PolynomialError: Sparse must be follow this exact syntax: 'x' + int.
        PolynomialError: Coefficients must be (contains) a list of int or float.
    """
    _dictionary = "abcdefghijklmnopqrstuvwxyz" # to use the alphabetic notation: f.a, f.b, f.c

    def __init__(self, *coefficients, name: str = "f", **sparse):
        """Basic class for implementing a Polynomial in object-oriented programming
        (see the Polynomial's docstring for more informations)

        Args:
            *coefficients (tuple, optional): An iterable containing the various coefficients
            of the polynomial, or a list of coefficients from largest to smallest.
            name (str, optional): function name. Defaults to "f".
            **sparse (dict, optional): An association of x with its associated coefficient: x34=5, x0=8, ....

        Raises:
            PolynomialError: Specify coefficients list or keyword terms, not both.
            PolynomialError: Sparse must be follow this exact syntax: 'x' + int.
            PolynomialError: Coefficients must be (contains) a list of int or float.
        """
        self._init(coefficients, sparse)
        self.name = name

    def _init(self, coefficients: tuple = (), sparse: dict = None):
        """Simple method for initializing a polynomial object after instantiating it.

        Args:
            coefficients (tuple, optional): An iterable containing the various coefficients of the polynomial, or a list of coefficients from largest to smallest. Defaults to ().
            sparse (dict, optional): An association of x with its associated coefficient: x34=5, x0=8, .... Defaults to {}.

        Raises:
            PolynomialError: Specify coefficients list or keyword terms, not both.
            PolynomialError: Sparse must be follow this exact syntax: 'x' + int.
            PolynomialError: Coefficients must be (contains) a list of int or float.
        """
        if sparse is None:
            sparse = {}

        if coefficients and sparse:
            raise PolynomialError("Specify coefficients list or keyword terms, not both.")

        if not all(len(x) > 1 and x[0] == "x" and x[1:].isdigit() for x in sparse):
            raise PolynomialError("Sparse must be follow this exact syntax: 'x' + int.")

        if not sparse:
            if len(coefficients) == 1 and issubclass(type(coefficients[0]), Iterable):
                coefficients = list(coefficients[0])
            else:
                coefficients = list(coefficients)

            if not all(isinstance(x, (int, float)) for x in coefficients):
                raise PolynomialError("Coefficients must be (contains) a list of int or float.")

            coefficients.reverse()

            self._sparse = {f"x{str(x)}": v for x, v in enumerate(coefficients)}

        else:
            self._sparse = sparse

        self._checkfloatsparse()

    @property
    def coefficients(self) -> list[int | float]:
        """The list of coefficients, including null coefficients.

        Examples:
            >>> f = Polynomial(-3, 23, -67)
            >>> f.coefficients
            [-3, 23, -67]

            >>> f.coefficients = [45, 0, -3.8, 1]
            >>> print(f)
            f(x) = 45x^3 - 3.8x + 1

            >>> f = Polynomial(0)
            >>> f.coefficients
            [0]

        """
        x: list = [coef for _, coef in self.items()]
        x = x[self.__class__._position_end_zeros(x):]

        if not x:
            x.append(0)

        return x

    @coefficients.setter
    def coefficients(self, _coef: Iterable):
        self._init((_coef))

    @coefficients.deleter
    def coefficients(self):
        raise PolynomialError("You can't delete coefficients attribute.")

    @property
    def sparse(self) -> dict[str, int]:
        """A dictionary comprising an association of the
        x power number and its associated coefficient.

        Examples:
            >>> f = Polynomial(-3, 23, -67)
            >>> f.sparse
            {'x2': -3, 'x1': 23, 'x0': -67}

            >>> f.sparse = {'x3':45, 'x1':-3.8, 'x0':1}
            >>> print(f)
            f(x) = 45x^3 - 3.8x + 1
        """
        items = list(self._sparse.items())
        items.sort(key= lambda x: int(x[0][1:]), reverse=True)
        n = self.__class__._position_end_zeros([value for _, value in items])

        self._sparse = dict(items[n:])
        return self._sparse

    @sparse.setter
    def sparse(self, _sparse: dict):
        self._init(sparse=_sparse)

    @sparse.deleter
    def sparse(self):
        raise PolynomialError("You can't delete sparse attribute.")

    @property
    def degree(self) -> int:
        """The maximum degree of the polynomial.

        Examples:
            >>> f = Polynomial(-3, 23, -67)
            >>> f.degree
            2

            >>> f[56] = 5
            >>> f
            Polynomial(x56=5, x2=-3, x=23, -67)
            >>> f.degree
            56
        """
        return self._degree()

    @degree.setter
    def degree(self, _):
        raise PolynomialError("You can't change the degree of a polynomial object.")

    @degree.deleter
    def degree(self):
        raise PolynomialError("You can't delete degree attribute.")

    @property
    def delta(self) -> float:
        """Polynomial discriminant, only available for second-degree
        polynomials, otherwise returns an error. Result of b² - 4ac

        Examples:
            >>> f = Polynomial(1, 1, -12)
            >>> print(f)
            f(x) = x² + x - 12
            >>> f.delta
            49.0
        """
        if self.degree == 2 and sp:
            a, b, c = sp.symbols("a b c")

            result = sp.simplify(b**2 - 4 * a * c)
            result = result.subs({a: self.a, b: self.b, c: self.c})

            return float(result.evalf())

        elif self.degree == 2:
            return self.b**2 - 4 * self.a * self.c

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd degree equations.")

    @delta.setter
    def delta(self, _):
        raise PolynomialError("You can't change the discriminant.")

    @delta.deleter
    def delta(self):
        raise PolynomialError("You can't delete delta attribute.")

    @property
    def alpha(self) -> int | float:
        """Alpha of the polynomial, only available for second-degree polynomials,
        otherwise returns an error. Result of -b / 2a.

        Examples:
            >>> f = Polynomial(1, 1, -12)
            >>> f.alpha, f.beta
            (-0.5, -12.25)
            >>> f.canonic()
            '(x + 0.5)² - 12.25'
        """
        if self._degree() == 2:
            return -self.b / (2 * self.a)

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd degree equations.")

    @alpha.setter
    def alpha(self, _):
        raise PolynomialError("You can't change it.")

    @alpha.deleter
    def alpha(self):
        raise PolynomialError("You can't delete alpha attribute.")

    @property
    def beta(self) -> int | float:
        """Beta of the polynomial, only available for second-degree polynomials,
        otherwise returns an error. Result of 4ac - b² / 4a.

        Examples:
            >>> f = Polynomial(1, 1, -12)
            >>> f.alpha, f.beta
            (-0.5, -12.25)
            >>> f.canonic()
            '(x + 0.5)² - 12.25'
        """
        if self._degree() == 2:
            return (4 * self.a * self.c - self.b**2) / 4 * self.a

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd degree equations.")

    @beta.setter
    def beta(self, _):
        raise PolynomialError("You can't change it.")

    @beta.deleter
    def beta(self):
        raise PolynomialError("You can't delete beta attribute.")

    def __getattr__(self, x: str) -> int | float:
        """Method for retrieving attributes using multiple syntaxes.

        Args:
            x (str): The attribute in question.

        Raises:
            PolynomialError: x must be follow this exact syntax: 'x' + int.
            PolynomialError: x must be in this list: a, b, c, d, ...
            PolynomialError: Polynomial must has coefficients in this order: self.degree, ..., 2, 1, 0.
            AttributeError: x was not found.

        Returns:
            int | float: Returns the associated coefficient.
        """
        if len(x) > 1 and x[0] == "x":
            if not x[1:].isdigit():
                raise PolynomialError(f"{x} must be follow this exact syntax: 'x' + int.")

            return 0 if x not in self.sparse.keys() else self.sparse[x]

        elif x in self.__class__._dictionary and self._degree() <= 27:
            dico = self.__class__._dictionary[:self._degree() + 1]

            if x not in dico:
                raise PolynomialError(f"{x} must be in this list: {", ".join(list(dico))}")

            try:
                list_x: list[int] = [int(x[1:]) for x, _ in self.items()]
                list_x.reverse()

                for i, _x in zip(range(self._degree() + 1), list_x):
                    assert i == _x

            except AssertionError as e:
                raise PolynomialError(
                    "Polynomial must has coefficients in this order: self.degree, ..., 2, 1, 0."
                ) from e

            return self.sparse[f"x{str(len(dico) - dico.index(x) - 1)}"]

        else:
            raise AttributeError(f"{x} was not found.")

    def __setattr__(self, x: str, value: Any):
        """Method for updating attributes.

        Args:
            x (str): The attribute in question.
            value (Any): The vew value.
        """
        if isinstance(value, (int, float)):
            if len(x) > 1 and x[0] == "x" and x[1:].isdigit():
                self.sparse[x] = value
                self._checkfloatsparse()

            elif x in self.__class__._dictionary and self._degree() <= 27:
                dico: str = self.__class__._dictionary[:self._degree() + 1]
                self.sparse[f"x{str(len(dico) - dico.index(x) - 1)}"] = value
                self._checkfloatsparse()

        else:
            super().__setattr__(x, value)

    def __delattr__(self, x: str):
        """Deletes the associated coefficient.

        Args:
            x (str): The attribute in question.
        """
        if len(x) > 1 and x[0] == "x" and x[1:].isdigit():
            del self.sparse[x]

        elif x in self.__class__._dictionary and self._degree() <= 27:
            dico: str = self.__class__._dictionary[:self._degree() + 1]
            del self.sparse[f"x{str(len(dico) - dico.index(x) - 1)}"]

        else:
            super().__delattr__(x)

    def __getitem__(self, x: int) -> int:
        """In the manner of a list, returns the coefficient associated with the input degree.

        Args:
            x (int): The degree that you want, like the index in a list.

        Returns:
            int: Returns the associated coefficient
        """
        try:
            return self.sparse[f"x{x}"]

        except KeyError:
            return 0

    def __setitem__(self, x: int, value: int):
        """In the manner of a list, set the new coefficient associated with the input degree.

        Args:
            x (int): The degree that you want, like the index in a list.
            value (int): The new value.
        """
        self.sparse[f"x{x}"] = value
        self._checkfloatsparse()

    def __delitem__(self, x: int):
        """In the manner of a list, del the input degree.

        Args:
            x (int): The degree that you want, like the index in a list.
        """
        del self.sparse[f"x{x}"]

    def __len__(self) -> int:
        """Returns the “length” of the polynomial, i.e. its degree + 1"""
        return len(self.sparse.keys())

    def __iter__(self) -> Iterator:
        """An iterator on the coefficients of the polynomial"""
        return iter(self.coefficients)

    def __reversed__(self) -> list[int | float]:
        """Reverses the order of coefficients in the polynomial"""
        coefficients = self.coefficients
        coefficients.reverse()

        return coefficients[self.__class__._position_end_zeros(coefficients):]

    def __contains__(self, x: str) -> bool:
        """Checks whether the input degree is less than the maximum degree of the polynomial

        Args:
            x (str): The x power.

        Raises:
            PolynomialError: x must be follow this exact syntax: 'x' + int.

        Returns:
            bool: Return this statement as a bool.
        """
        if not (len(x) > 1 and x[0] == "x" and x[1:].isdigit()):
            raise PolynomialError(f"{x} must be follow this exact syntax: 'x' + int.")

        return x in self.sparse.keys()

    def __bool__(self) -> bool:
        """Determinate if the polynomial is not null"""
        return not (len(self.coefficients) == 1 and not self.coefficients[0])

    def __call__(self, x: int | float) -> int | float:
        """Method used to give a mathematical notation when calling a function, such as f(6), or f(-1).

        Args:
            x (int | float): The x-axis.

        Returns:
            int | float: The y-axis.
        """
        string = " + ".join(
            f"{str(self.sparse[x])} * x**{x[1:]}"
            for x, _ in self.items()
            if self.sparse[x] != 0
        )
        string = string.replace("+ -", "- ").replace("x^1 ", "x").replace("x^0", "").replace("1x", "x")

        return eval(string)

    def __int__(self) -> int:
        """Return the degree of the polynomial"""
        return self.degree

    def __str__(self) -> str:
        """Returns the expanded form of the function"""
        return f"{self.name}(x) = {self.developed()}"

    def __repr__(self) -> str:
        """Returns the formal form of the function"""
        string: str = ", ".join(f"{x}={str(v)}" for x, v in self.items() if v != 0).replace("x1=", "x=").replace("x0=", "") or '0'

        return f'{self.__class__.__qualname__}({string})'

    def __eq__(self, other: Self) -> bool:
        """Tests whether two polynomials are equal by comparing their internal sparse dictionary"""
        if isinstance(other, self.__class__):
            return self.sparse == other.sparse

        else:
            raise NotImplementedError

    def __ne__(self, other: Self) -> bool:
        """Tests whether two polynomials are not equal by comparing their internal sparse dictionary."""
        if isinstance(other, self.__class__):
            return self.sparse != other.sparse

        else:
            raise NotImplementedError

    def __add__(self, other: Self | (int | float)) -> Self:
        """Adds two polynomials or a polynomial and a number"""
        if isinstance(other, self.__class__):
            P = self.__class__(name=self.name)

            for x, value in self.items():
                if other.sparse.get(x) is None:
                    P.sparse[x] = value

                else:
                    P.sparse[x] = value + other.sparse[x]
                    del other.sparse[x]

            for x, value in other.items():
                P.sparse[x] = value

            for x, value in P.items():
                if not value:
                    del P.sparse[x]

            return P

        elif isinstance(other, (int, float)):
            P = self.copy()

            if P.sparse.get("x0") is None:
                P[0] = 0

            P[0] += other

            return P

        else:
            raise NotImplementedError

    __radd__ = __iadd__ = __add__

    def __sub__(self, other: Self | (int | float)) -> Self:
        """Subs two polynomials or a polynomial and a number."""
        if isinstance(other, self.__class__):
            P = self.__class__(name=self.name)

            for x, value in self.items():
                if other.sparse.get(x) is None:
                    P.sparse[x] = value

                else:
                    P.sparse[x] = value - other.sparse[x]
                    del other.sparse[x]

            for x, value in other.items():
                P.sparse[x] = -value

            for x, value in P.items():
                if not value:
                    del P.sparse[x]

            return P

        elif isinstance(other, (int, float)):
            P = self.copy()

            if P.sparse.get("x0") is None:
                P[0] = 0

            P[0] -= other

            return P

        else:
            raise NotImplementedError

    __rsub__ = __isub__ = __sub__

    def __mul__(self, other: Self | (int | float)) -> Self:
        """Multiply two polynomials or a polynomial and a number."""
        if isinstance(other, self.__class__):
            liste, n = [], 0

            for i in range(self.degree, -1, -1):
                if not self[i]:
                    continue

                liste.append(self.__class__())
                for j in range(other.degree, -1, -1):
                    liste[n][i + j] = self[i] * other[j]

                n += 1

            return sum(liste)

        elif isinstance(other, (int, float)):
            return self.__class__(**{x:value * other for x, value in self.items()})

        else:
            raise NotImplementedError

    __rmul__ = __imul__ = __mul__

    def __truediv__(self, other: Self | (int | float)) -> Self:
        """Divs it's self and a number or a polynomial."""
        if isinstance(other, (int, float)):
            if not other:
                raise ZeroDivisionError("You can't divide a polynomial by zero.")

            return self.__mul__(other**-1)

        elif isinstance(other, self.__class__):
            return self.__floordiv__(other)

        else:
            raise NotImplementedError

    __rtruediv__ = __itruediv__ = __truediv__

    def __floordiv__(self, other: Self | (int | float)) -> Self:
        """Return the quotient if the division is between two polynomial, or the result of a normal division."""
        if isinstance(other, self.__class__):
            return divmod(self, other)[0]

        elif isinstance(other, (int, float)):
            return self.__truediv__(other)

        else:
            raise NotImplementedError

    __ifloordiv__ = __floordiv__

    def __mod__(self, other: Self | (int | float)) -> Self:
        """Return the remainder if the division is between two polynomial, or a null polynomial."""
        if isinstance(other, self.__class__):
            return divmod(self, other)[1]

        elif isinstance(other, (int, float)):
            return self.__class__()

        else:
            raise NotImplementedError

    def __divmod__(self, other: Self | (int | float)) -> Self:
        """Return the quotient and the remainder."""
        if isinstance(other, self.__class__):
            if not other:
                raise ZeroDivisionError("You can't divide a polynomial by an other null polynomial.")

            elif self.degree < other.degree:
                return (self.__class__(), self.copy())

            else:
                msg_out = list(self.copy())
                other = list(other)
                normalizer = other[0]

                for i in range(len(self) - (len(other)-1)):
                    msg_out[i] /= normalizer
                    coef = msg_out[i]
                    if coef != 0:
                        for j in range(1, len(other)):
                            if other[j] != 0:
                                msg_out[i + j] += -other[j] * coef

                separator = -(len(other)-1)
                return self.__class__(msg_out[:separator]), self.__class__(msg_out[separator:])

        elif isinstance(other, (int, float)):
            return self / other

        else:
            raise NotImplementedError

    def __pow__(self, other: int) -> Self:
        """Implemented the power of a polynomial, with a positive integer"""
        if isinstance(other, int) and other >= 0:
            P = self.copy()

            for _ in range(other - 1):
                P *= self

            return P

        else:
            raise NotImplementedError("The power must be a positive and an integer number.")

    __ipow__ = __pow__

    def __invert__(self) -> Self:
        raise NotImplementedError

    def __neg__(self) -> Self:
        """Multiply the polynomial by -1"""
        P = self.copy()
        P *= -1

        return P

    def __pos__(self) -> Self:
        """Return it self"""
        return self

    def __hash__(self) -> int:
        """Return the hash from the dictionary"""
        return hash(self.sparse)

    def copy(self) -> Self:
        """Copy the polynomial
        
        Examples:
            >>> h = Polynomial(9, -30, 25, name="h")
            >>> print(h)
            h(x) = 9x² - 30x + 25

            >>> _h = h.copy()
            >>> print(_h)
            h(x) = 9x² - 30x + 25
        """
        P = self.__class__(name=self.name)
        P._init(sparse=self.sparse)

        return P

    def _degree(self) -> int:
        """Returns the updated degree of the polynomial

        Examples:
            >>> f = Polynomial(-3, 23, -67)
            >>> f.degree
            2

            >>> f[56] = 5
            >>> f
            Polynomial(x56=5, x2=-3, x=23, -67)
            >>> f.degree
            56
        """
        if len(self.coefficients) == 1 and not self.coefficients[0]:
            return 0

        else:
            return max(int(x[1:]) for x in self.sparse.keys())

    @staticmethod
    def _position_end_zeros(iterable: Iterable[int | float]) -> int:
        """From the list of coefficients, return the position furthest to
        the right from the left where the zeros cancel the powers of x.

        Args:
            iterable (Iterable): The ordered list of coefficients in descending order.

        Examples:
            >>> l = [0 for _ in range(5)] + [6, 0, -3, 1.5]
            >>> l
            [0, 0, 0, 0, 0, 6, 0, -3, 1.5]
            >>> n = Polynomial._position_end_zeros(l)
            >>> n
            5
            >>> l[n:]
            [6, 0, -3, 1.5]

        Returns:
            int: The position from left.
        """
        n = 0
        for x in iterable:
            if not x:
                n += 1

            else:
                break

        return n

    @staticmethod
    def _checkfloat(iterable: Iterable[int | float]) -> Iterable[int | float]:
        """A function to change float coefficient in integers, when it's possible.

        Args:
            iterable (Iterable[int  |  float]): The iterator containing int or float numbers.

        Example:
            >>> l = [5.0, -3.67, 1.0, 98.1]
            >>> Polynomial._checkfloat(l)
            [5, -3.67, 1, 98.1]

        Returns:
            Iterable[int | float]: The checked iterator.
        """
        _iter: list[int | float] = []

        for coef in iterable:
            if coef:
                n = str(coef).find('.')

                if (len(str(coef)[n:]) == 2 and not int(str(coef)[-1])) or n == -1:
                   _iter.append(int(coef))

                else:
                    _iter.append(coef)

            else:
                _iter.append(0)

        return _iter

    def _checkfloatsparse(self):
        """A function to update coefficients attributes. By default call in the _init, and sets method.

        Examples:
            >>> f = Polynomial(5.0, -3.67, 1.0, 98.1)
            >>> print(f)
            f(x) = 5x^3 - 3.67x² + x + 98.1
        """
        iterator: list[int | float] = reversed([self[i] for i in range(self.degree, -1, -1)])

        for i, coef in enumerate(self.__class__._checkfloat(iterator)):
            if coef:
                self._sparse[f"x{str(i)}"] = coef

    def items(self) -> list[tuple]:
        """Returns sparse but ordered dictionary items in descending order
        
        Examples:
            >>> f = Polynomial(-3, 23, -67)
            >>> f.items() # just ordered
            [('x2', -3), ('x1', 23), ('x0', -67)]
            >>> dict(f.items()) == f.sparse
            True
        """
        items = list(self.sparse.items())
        items.sort(key= lambda x: int(x[0][1:]), reverse=True)

        return items

    def solve(self) -> tuple[float]:
        """Solve the equation f(x) = 0, and can use SymPy to get roots.
        
        Examples:
        
            * Solve 1st degrees equations
            
            >>> f = Polynomial(2, -26)
            >>> f.solve()
            13.0

            * Solve 2nd degrees equations

            With 2 roots:

            >>> f = Polynomial(1, 1, -12)
            >>> print(f)
            f(x) = x² + x - 12
            >>> f.delta
            49.0
            >>> f.solve()
            (-4.0, 3.0)

            With 1 root:
            
            >>> h = Polynomial(9, -30, 25, name="h")
            >>> print(h)
            h(x) = 9x² - 30x + 25
            >>> h.delta
            0.0
            >>> h.solve()
            (1.6666666666666667,)

            * Solve any equations with SymPy

            >>> f = Polynomial(x5=3, x2=23, x0=14)
            >>> print(f)
            f(x) = 3x^5 + 23x² + 14
            >>> solutions = f.solve()
            >>> solutions
            ((-2.061777647931321+0j), (-0.023547692027920532-0.7769582400529271j), (-0.023547692027920532+0.7769582400529271j), (1.054436515993581-1.6230188809835164j), (1.054436515993581+1.6230188809835164j))

        """
        delta = self.delta if self.degree == 2 else None

        if self.degree == 1:
            return (-self.b / self.a)

        elif self.degree == 2 and delta > 0:
            solutions = [(-self.b - delta**0.5) / (2 * self.a), (-self.b + delta**0.5) / (2 * self.a)]

            if solutions[0] > solutions[1]:
                solutions[0], solutions[1] = solutions[1], solutions[0]

            return tuple(solutions)

        elif self.degree == 2 and delta == 0:
            return (-self.b / (2 * self.a), )

        elif sp:
            x = sp.symbols("x")
            expression = eval(
                " + ".join(
                    f"{str(value)} * x**{_x[1:]}" for _x, value in self.items()
                ),
                {"x": x},
            )

            return tuple(complex(sol.evalf()) for sol in sp.solve(expression))

        else:
            raise PolynomialError("Solve an equation greater than the 2nd degree isn't implemented. Please install the 'sympy' library if you want to solve them.")

    def developed(self) -> str:
        """Returns the function's expanded form

        Examples:
            >>> f = Polynomial(-3, 23, -67)
            >>> f.developed()
            '-3x² + 23x - 67'
        """
        string = " + ".join(
            f"{str(self.sparse[x])}x^{x[1:]}"
            for x, _ in self.items()
            if self.sparse[x] != 0
        )
        string = string.replace("+ -", "- ").replace("x^1 ", "x ").replace("x^0", "").replace(" 1x", " x").replace("-1x", "-x").replace("x^2 ", "x² ")

        if string.endswith("^1"):
            string = string[:len(string) - 2]

        if string.startswith("1x"):
            string = string[1:]

        if not string:
            string = "0"

        return string

    def canonic(self, decimal: int = 3) -> str:
        """Returns the canonical form, only if the polynomial is of the second degree

        Args:
            decimal (int, optional): The number of decimal if coefficients, alpha or beta are floating numbers. Defaults to 3.

        Examples:
            >>> f = Polynomial(1, 1, -12)
            >>> f.alpha, f.beta
            (-0.5, -12.25)
            >>> f.canonic()
            '(x + 0.5)² - 12.25'

        Raises:
            NotImplementedError: The polynomial is not of the second degree.

        Returns:
            str: A string of the canonical form.
        """
        if self.degree != 2:
            raise NotImplementedError("For the moment, only implemented for 2nd degree equations.")

        string = f"{self.a}(x - {round(self.alpha, decimal)})² + {round(self.beta, decimal)}"
        string = string.replace("- -", "+ ").replace("+ -", "- ").replace("-1(", "-(").replace("² + 0", "").replace("² - 0", "")

        if string.startswith("1("):
            string = string[1:]

        return string

    def factorised(self, decimal: int = 3) -> str:
        """Returns the factorized form, only if the polynomial is of the second degree

        Args:
            decimal (int, optional): The number of decimal if coefficients are floating numbers. Defaults to 3.

        Examples:
            >>> f = Polynomial(1, 1, -12)
            >>> f.factorised()
            '(x + 4.0)(x - 3.0)'

            >>> h = Polynomial(9, -30, 25, name="h")
            >>> h.factorised()
            '9(x - 1.667)²'

        Raises:
            PolynomialError: Causes an error if the polynomial does not pass on the x-axis.
            NotImplementedError: The polynomial is not of the second degree.

        Returns:
            str: _description_
        """
        if self.degree != 2:
            raise NotImplementedError("For the moment, only implemented for 2nd degree equations.")

        solutions = self.solve()

        if solutions is None:
            raise PolynomialError("Factorised expression doesn't exist.")

        elif len(solutions) == 1:
            string = f"{self.a}(x - {round(solutions[0], decimal)})²"

        elif len(solutions) == 2:
            string = f"{self.a}(x - {round(solutions[0], decimal)})(x - {round(solutions[1], decimal)})"

        string = string.replace("- -", "+ ").replace("-1(", "-(")

        if string.startswith("1("):
            string = string[1:]

        return string

    def derive(self) -> Self:
        """Returns the derivative of the polynomial
        
        Examples:
            >>> f = Polynomial(-3, 23, -67)
            >>> f
            Polynomial(x2=-3, x=23, -67)

            >>> h = f.derive()
            >>> h
            Polynomial(x=-134, 23)

        """
        L = len(self)-1
        return self.__class__([(L-i) * self[i] for i in range(L)])
