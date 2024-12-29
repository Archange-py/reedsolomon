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

import sympy as sp


class PolynomialError(Exception):
    """Base class for handling polynomial errors"""
    pass


class Polynomial:
    _dictionnary = "abcdefghijklmnopqrstuvwxyz" # to use the notation f.a, f.b, f.c for example

    def __init__(self, *coefficients, name: str = "f", **sparse):
        """Basic class for implementing a Polynomial in object-oriented programming

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

    def _init(self, coefficients: tuple = (), sparse: dict = {}):
        """Simple method for initializing a polynomial object after instantiating it.

        Args:
            coefficients (tuple, optional): An iterable containing the various coefficients
            of the polynomial, or a list of coefficients from largest to smallest. Defaults to ().
            sparse (dict, optional): An association of x with its associated coefficient: x34=5, x0=8, .... Defaults to {}.

        Raises:
            PolynomialError: Specify coefficients list or keyword terms, not both.
            PolynomialError: Sparse must be follow this exact syntax: 'x' + int.
            PolynomialError: Coefficients must be (contains) a list of int or float.
        """
        if coefficients and sparse:
            raise PolynomialError("Specify coefficients list or keyword terms, not both.")

        if not all([len(x) > 1 and x[0] == "x" and x[1:].isdigit() for x in sparse.keys()]):
            raise PolynomialError(f"Sparse must be follow this exact syntax: 'x' + int.")

        if not sparse:
            if len(coefficients) == 1 and issubclass(type(coefficients[0]), Iterable):
                coefficients = list(coefficients[0])
            else:
                coefficients = list(coefficients)

            if not all([isinstance(x, (int, float)) for x in coefficients]):
                raise PolynomialError(f"Coefficients must be (contains) a list of int or float.")

            coefficients.reverse()

            self._sparse = {"x" + str(x):v for x, v in enumerate(coefficients)}

        else:
            self._sparse = sparse

    @property
    def coefficients(self) -> list[int | float]:
        """The list of coefficients, including null coefficients"""
        x = [coef for _, coef in self.items()]

        return x[self._position_end_zeros(x):]

    @coefficients.setter
    def coefficients(self, _coef: Iterable):
        self._init((_coef))

    @coefficients.deleter
    def coefficients(self):
        raise PolynomialError("You can't delete coefficients attribute.")

    @property
    def sparse(self) -> dict[str, int]:
        """A dictionary comprising an association of the
        x power number and its associated coefficient."""
        items = list(self._sparse.items())
        items.sort(key= lambda x: int(x[0][1:]), reverse=True)
        n = self._position_end_zeros([value for _, value in items])

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
        """The maximum degree of the polynomial."""
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
        polynomials, otherwise returns an error. Result of b² - 4ac"""
        if self._degree() == 2:
            a, b, c = sp.symbols("a b c")

            result = sp.simplify(b**2 - 4 * a * c)
            result = result.subs({a: self.a, b: self.b, c: self.c})

            return float(result.evalf())

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
        otherwise returns an error. Result of -b / 2a."""
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
        otherwise returns an error. Result of 4ac - b² / 4a."""
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
            if x[1:].isdigit():
                if not x in self.sparse.keys():
                    return 0
 
                return self.sparse[x]

            else:
                raise PolynomialError(f"{x} must be follow this exact syntax: 'x' + int.")

        elif x in Polynomial._dictionnary and self._degree() <= 27:
            dico = Polynomial._dictionnary[:self._degree() + 1]

            if x not in dico:
                raise PolynomialError(f"{x} must be in this list: {", ".join(list(dico))}")

            try:
                list_x: list[int] = [int(x[1:]) for x, _ in self.items()]
                list_x.reverse()

                for i, _x in zip(range(self._degree() + 1), list_x):
                    assert i == _x

            except AssertionError:
                raise PolynomialError(f"Polynomial must has coefficients in this order: self.degree, ..., 2, 1, 0.")

            return self.sparse["x" + str(len(dico) - dico.index(x) - 1)]

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

            elif x in Polynomial._dictionnary and self._degree() <= 27:
                dico: str = Polynomial._dictionnary[:self._degree() + 1]
                self.sparse["x" + str(len(dico) - dico.index(x) - 1)] = value

        else:
            super().__setattr__(x, value)

    def __delattr__(self, x: str):
        """Deletes the associated coefficient.

        Args:
            x (str): The attribute in question.
        """
        if len(x) > 1 and x[0] == "x" and x[1:].isdigit():
            del self.sparse[x]

        elif x in Polynomial._dictionnary and self._degree() <= 27:
            dico: str = Polynomial._dictionnary[:self._degree() + 1]
            del self.sparse["x" + str(len(dico) - dico.index(x) - 1)]

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
            return self.sparse["x" + str(x)]

        except KeyError:
            return 0

    def __setitem__(self, x: int, value: int):
        """In the manner of a list, set the new coefficient associated with the input degree.

        Args:
            x (int): The degree that you want, like the index in a list.
            value (int): The new value.
        """
        self.sparse["x" + str(x)] = value

    def __delitem__(self, x: int):
        """In the manner of a list, del the input degree.

        Args:
            x (int): The degree that you want, like the index in a list.
        """
        del self.sparse["x" + str(x)]

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

        return coefficients[self._position_end_zeros(coefficients):]

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

    def __call__(self, x: int | float) -> int | float:
        """Method used to give a mathematical notation when calling a function, such as f(6), or f(-1).

        Args:
            x (int | float): The x-axis.

        Returns:
            int | float: The y-axis.
        """
        string = " + ".join(str(self.sparse[x]) + " * x**" + x[1:] for x, _ in self.items() if self.sparse[x] != 0)
        string = string.replace("+ -", "- ").replace("x^1 ", "x").replace("x^0", "").replace("1x", "x")

        return eval(string)

    def __str__(self) -> str:
        """Returns the expanded form of the function"""
        return f"{self.name}(x) = {self.developped()}"

    def __repr__(self) -> str:
        """Returns the formal form of the function"""
        string = ", ".join(x + "=" + str(v) for x, v in self.items())
        string = string.replace("x1=", "x=").replace("x0=", "").replace(", 0", "")

        if not string:
            string = "0"

        return f'{self.__class__.__qualname__}({string})'

    def __eq__(self, other: Self) -> bool:
        """Tests whether two polynomials are equal by comparing their internal sparse dictionary"""
        if isinstance(other, Polynomial):
            return self.sparse == other.sparse

        else:
            raise NotImplementedError

    def __ne__(self, other: Self) -> bool:
        """Tests whether two polynomials are not equal by comparing their internal sparse dictionary."""
        if isinstance(other, Polynomial):
            return self.sparse != other.sparse

        else:
            raise NotImplementedError

    def __add__(self, other: Self | (int | float)) -> Self:
        """Adds two polynomials or a polynomial and a number"""
        if isinstance(other, Polynomial):
            P = Polynomial(name=self.name)

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
        if isinstance(other, Polynomial):
            P = Polynomial(name=self.name)

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
        """Muls two polynomials or a polynomial and a number."""
        if isinstance(other, Polynomial):
            liste, n = [], 0

            for i in range(self.degree, -1, -1):
                if self[i]:
                    liste.append(Polynomial())
                    for j in range(other.degree, -1, -1):
                        liste[n][i + j] = self[i] * other[j]

                else:
                    continue

                n += 1

            return sum(liste)

        elif isinstance(other, (int, float)):
            return Polynomial(**{x:value * other for x, value in self.items()})

        else:
            raise NotImplementedError

    __rmul__ = __imul__ = __mul__

    def __truediv__(self, other: Self | (int | float)) -> Self:
        """Divs it's self and a number."""
        if isinstance(other, (int, float)):
            return self.__mul__(other**-1)

        else:
            raise NotImplementedError

    __rtruediv__ = __itruediv__ = __truediv__

    def __neg__(self) -> Self:
        """Multiply the polynomial by -1"""
        P = self.copy()
        P *= -1

        return P

    def __pos__(self) -> Self:
        """Return it self"""
        return self

    def __hash__(self) -> int:
        """Return the hash from the dictionary."""
        return hash(self.sparse)

    def copy(self) -> Self:
        """Copy the polynomial"""
        P = Polynomial(name=self.name)
        P._init(sparse=self.sparse)

        return P

    def _degree(self) -> int:
        """Returns the updated degree of the polynomial"""
        return max([int(x[1:]) for x in self.sparse.keys()])

    def _position_end_zeros(self, iterable: Iterable) -> int:
        """From the list of coefficients, return the position furthest to
        the right from the left where the zeros cancel the powers of x.

        Args:
            iterable (Iterable): The oredered list of coefficients in descending order.

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

    def items(self) -> list[tuple]:
        """Returns sparse but ordoned dictionary items in descending order"""
        items = list(self.sparse.items())
        items.sort(key= lambda x: int(x[0][1:]), reverse=True)

        return items

    def solve(self) -> tuple[float]:
        """Solve the equation f(x) = 0"""
        delta = self.delta if self.degree == 2 else None

        if self.degree == 1:
            return (-self.b / self.a)

        elif self.degree == 2 and delta > 0:
            solutions = [(-self.b - delta**0.5) / (2 * self.a), (-self.b + delta**0.5) / (2 * self.a)]

            if solutions[0] > solutions[1]:
                solutions[0], solutions[1] = solutions[1], solutions[0]

            return tuple(solutions)

        elif self.degree == 2 and delta == 0:
            return tuple([-self.b / (2 * self.a)])

        else:
            x = sp.symbols("x")

            expression = eval(" + ".join(str(value) + " * x**" + _x[1:] for _x, value in self.items()), {"x":x})
            solutions = tuple(complex(sol.evalf()) for sol in sp.solve(expression))

            return solutions

    def developped(self) -> str:
        """Returns the function's expanded form"""
        string = " + ".join(str(self.sparse[x]) + "x^" + x[1:] for x, _ in self.items() if self.sparse[x] != 0)
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

        Raises:
            NotImplementedError: The polynomial is not of the second degree.

        Returns:
            str: A string of the canonical form.
        """
        if self.degree == 2:
            string = f"{self.a}(x - {round(self.alpha, decimal)})² + {round(self.beta, decimal)}"
            string = string.replace("- -", "+ ").replace("+ -", "- ").replace("-1(", "-(").replace("² + 0", "").replace("² - 0", "")

            if string.startswith("1("):
                string = string[1:]

            return string

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd degree equations.")

    def factorised(self, decimal: int = 3) -> str:
        """Returns the factorized form, only if the polynomial is of the second degree

        Args:
            decimal (int, optional): The number of decimal if coefficients are floating numbers. Defaults to 3.

        Raises:
            PolynomialError: Causes an error if the polynomial does not pass on the x-axis.
            NotImplementedError: The polynomial is not of the second degree.

        Returns:
            str: _description_
        """
        if self.degree == 2:
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

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd degree equations.")

    def derive(self) -> Self:
        """Returns the derivative of the polynomial"""
        L = len(self)-1
        return Polynomial([(L-i) * self[i] for i in range(L)])
