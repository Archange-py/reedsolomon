
class GF(set):
    def __init__(self, length: int):
        self.length = length
        super().__init__(range(self.length))

    def inv(self, nbr: int) -> int:
        pass

    def add(self, nbr1: int, nbr2: int) -> int:
        return (nbr1 + nbr2) % self.length

    def sub(self, nbr1: int, nbr2: int) -> int:
        return (nbr1 - nbr2) % self.length

    def mul(self, nbr1: int, nbr2: int) -> int:
        return (nbr1 * nbr2) % self.length

    def div(self, nbr1: int, nbr2: int) -> int:
        return self.mul(nbr1, self.inv(nbr2))
