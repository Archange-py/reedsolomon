#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

"""

from typing import Iterable, Iterator, Self, Any


class PolynomialError(Exception):
    pass


class Polynomial:
    _dictionnary = "abcdefghijklmnopqrstuvwxyz"

    def __init__(self, *coefficients, name: str = "f", **sparse):
        self._init(coefficients, sparse)
        self.name = name

    def _init(self, coefficients: tuple = (), sparse: dict = {}):
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
    def coefficients(self) -> list[int]:
        x = [coef for _, coef in self.items()]

        n = 0
        for coef in x:
            if not coef:
                n += 1

            else:
                break

        return x[n:]

    @coefficients.setter
    def coefficients(self, _coef: Iterable):
        self._init((_coef))

    @coefficients.deleter
    def coefficients(self):
        raise PolynomialError("You can't delete coefficients attribute.")

    @property
    def sparse(self) -> dict[str, int]:
        items = list(self._sparse.items())
        items.sort(key= lambda x: int(x[0][1:]), reverse=True)

        n = 0
        for _, value in items:
            if not value:
                n += 1

            else:
                break

        self._sparse = dict(items[n:])
        return self._sparse

    @sparse.setter
    def sparse(self, _sparse: dict):
        self._init(sparse=_sparse)

    @sparse.deleter
    def sparse(self):
        raise PolynomialError("You can't delete sparse attribute.")

    @property
    def degree(self) -> dict[str, int]:
        return self._degree()

    @degree.setter
    def degree(self, _):
        raise PolynomialError("You can't change the degree of a polynomial object.")

    @degree.deleter
    def degree(self):
        raise PolynomialError("You can't delete degree attribute.")

    @property
    def delta(self) -> int | float:
        if self._degree() == 2:
            return self.b**2 - 4 * self.a * self.c

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd _degree equations.")

    @delta.setter
    def delta(self, _):
        raise PolynomialError("You can't change the discriminant.")

    @delta.deleter
    def delta(self):
        raise PolynomialError("You can't delete delta attribute.")


    @property
    def alpha(self) -> int | float:
        if self._degree() == 2:
            return -self.b / (2 * self.a)

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd _degree equations.")

    @alpha.setter
    def alpha(self, _):
        raise PolynomialError("You can't change it.")

    @alpha.deleter
    def alpha(self):
        raise PolynomialError("You can't delete alpha attribute.")

    @property
    def beta(self) -> int | float:
        if self._degree() == 2:
            return (4 * self.a * self.c - self.b**2) / 4 * self.a

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd _degree equations.")

    @beta.setter
    def beta(self, _):
        raise PolynomialError("You can't change it.")

    @beta.deleter
    def beta(self):
        raise PolynomialError("You can't delete beta attribute.")

    def __getattr__(self, x: str) -> int | float:
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
                raise PolynomialError(f"Polynomial must has coefficients in this order: self._degree(), ..., 2, 1, 0.")

            return self.sparse["x" + str(len(dico) - dico.index(x) - 1)]

        else:
            raise AttributeError(f"{x} was not found.")

    def __setattr__(self, x: str, value: Any):
        if isinstance(value, (int, float)):
            if len(x) > 1 and x[0] == "x" and x[1:].isdigit():
                self.sparse[x] = value

            elif x in Polynomial._dictionnary and self._degree() <= 27:
                dico: str = Polynomial._dictionnary[:self._degree() + 1]
                self.sparse["x" + str(len(dico) - dico.index(x) - 1)] = value

        else:
            super().__setattr__(x, value)

    def __delattr__(self, x: str):
        if len(x) > 1 and x[0] == "x" and x[1:].isdigit():
            del self.sparse[x]

        elif x in Polynomial._dictionnary and self._degree() <= 27:
            dico: str = Polynomial._dictionnary[:self._degree() + 1]
            del self.sparse["x" + str(len(dico) - dico.index(x) - 1)]

        else:
            super().__delattr__(x)

    def __getitem__(self, x: int) -> int:
        try:
            return self.sparse["x" + str(x)]

        except KeyError:
            return 0

    def __setitem__(self, x: int, value: int):
        self.sparse["x" + str(x)] = value

    def __delitem__(self, x: int):
        del self.sparse["x" + str(x)]

    def __len__(self) -> int:
        return len(self.sparse.keys())

    def __iter__(self) -> Iterator:
        return iter(self.coefficients)

    def __contains__(self, x: str) -> bool:
        if not (len(x) > 1 and x[0] == "x" and x[1:].isdigit()):
            raise PolynomialError(f"{x} must be follow this exact syntax: 'x' + int.")

        return x in self.sparse.keys()

    def __call__(self, x: int | float) -> int | float:
        string = " + ".join(str(self.sparse[x]) + " * x**" + x[1:] for x, _ in self.items() if self.sparse[x] != 0)
        string = string.replace("+ -", "- ").replace("x^1 ", "x").replace("x^0", "").replace("1x", "x")

        return eval(string)

    def __str__(self) -> str:
        return f"{self.name}(x) = {self.devellopped()}"

    def __repr__(self) -> str:
        string = ", ".join(x + "=" + str(v) for x, v in self.items())
        string = string.replace("x1=", "x=").replace("x0=", "").replace(", 0", "")

        if not string:
            string = "0"

        return f'{self.__class__.__qualname__}({string})'

    def __eq__(self, other: Self) -> bool:
        if isinstance(other, Polynomial):
            return self.sparse == other.sparse

        else:
            raise NotImplementedError

    def __ne__(self, other: Self) -> bool:
        if isinstance(other, Polynomial):
            return self.sparse != other.sparse

        else:
            raise NotImplementedError

    def __add__(self, other: Self | (int | float)) -> Self:
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
        if isinstance(other, (int, float)):
            return Polynomial(**{x:value * other for x, value in self.items()})

        else:
            raise NotImplementedError

    __rmul__ = __imul__ = __mul__

    def __truediv__(self, other: Self | (int | float)) -> Self:
        if isinstance(other, (int, float)):
            return self.__mul__(other**-1)

        else:
            raise NotImplementedError

    __rtruediv__ = __itruediv__ = __truediv__

    def __neg__(self) -> Self:
        P = self.copy()
        P *= -1

        return P

    def __pos__(self) -> Self:
        return self

    def __hash__(self) -> int:
        return hash(self.sparse)

    def copy(self) -> Self:
        P = Polynomial(name=self.name)
        P._init(sparse=self.sparse)

        return P

    def _degree(self) -> int:
        return max([int(x[1:]) for x in self.sparse.keys()])

    def items(self) -> list[tuple]:
        items = list(self.sparse.items())
        items.sort(key= lambda x: int(x[0][1:]), reverse=True)

        return items

    def solve(self) -> tuple[int |float] | None:
        delta = self.delta

        if delta > 0:
            solutions = [(-self.b - delta**0.5) / (2 * self.a), (-self.b + delta**0.5) / (2 * self.a)]

            if solutions[0] > solutions[1]:
                solutions[0], solutions[1] = solutions[1], solutions[0]

            return tuple(solutions)

        elif delta == 0:
            return tuple([-self.b / (2 * self.a)])

        else:
            return None

    def devellopped(self) -> str:
        string = " + ".join(str(self.sparse[x]) + "x^" + x[1:] for x, _ in self.items() if self.sparse[x] != 0)
        string = string.replace("+ -", "- ").replace("x^1 ", "x ").replace("x^0", "").replace(" 1x", " x").replace("-1x", "-x").replace("-1x^", "-x^").replace("x^2 ", "x² ")

        if string.endswith("^1"):
            string = string[:len(string) - 2]

        if not string:
            string = "0"

        return string

    def canonic(self, decimal: int = 3) -> str:
        if self._degree() == 2:
            string = f"{self.a}(x - {round(self.alpha, decimal)})² + {round(self.beta, decimal)}"
            string = string.replace("- -", "+ ").replace("+ -", "- ").replace("1(", "(").replace("-1(", "-(").replace("² + 0", "").replace("² - 0", "")

            return string

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd _degree equations.")

    def factorised(self, decimal: int = 3) -> str:
        if self._degree() == 2:
            solutions = self.solve()

            if solutions is None:
                raise PolynomialError("Factorised expression doesn't exist.")

            elif len(solutions) == 1:
                string = f"{self.a}(x - {round(solutions[0], decimal)})²"
                string = string.replace("- -", "+ ").replace("1(", "(").replace("-1(", "-(")

            elif len(solutions) == 2:
                string = f"{self.a}(x - {round(solutions[0], decimal)})(x - {round(solutions[1], decimal)})"
                string = string.replace("- -", "+ ").replace("1(", "(").replace("-1(", "-(")

            return string

        else:
            raise NotImplementedError("For the moment, only implemented for 2nd _degree equations.")

    def derive(self) -> Self:
        L = len(self)-1
        return Polynomial([(L-i) * self[i] for i in range(L)])
