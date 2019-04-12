class Matrix(list):
    def __init__(self, array):
            list.__init__(self, array)

    def __str__(self) -> str:
        # fix print so lines with bigger numbers dont stick out
        first_line = self.matricize(0, ('\u23A1', '\u23A4'))

        all_lines = first_line
        for i in range(1, self.rows - 1):
            middle_line = self.matricize(i)
            all_lines += middle_line
        if self.rows != 1:
            last_line = self.matricize(-1, ('\u23A3', '\u23A6'))
            all_lines += last_line

        return all_lines

    def __mul__(self, other):
        try:
            return self.dot_product(other)
        except (ValueError, TypeError):
            try:
                return self.scalar_product(other)
            except ValueError:
                pass

    # add .from_document() generation method

    # Could expand into seperate class that continues generation
    @classmethod
    def from_alg(cls, alg, rows, cols=None):
        """Takes a function alg(i, j) and uses it to generate values"""
        if cols is None:
            cols = rows

        matrix = cls.generate_null(rows, cols)

        for i in range(rows):
            for j in range(cols):
                matrix[i][j] = alg(i, j)

        return matrix

    @classmethod
    def from_matrix(cls, matrix):
        from copy import deepcopy
        return deepcopy(matrix)

    @classmethod
    def generate_null(cls, rows, cols=None):
        if cols is None:
            cols = rows

        return cls([[0 for j in range(cols)]
                    for i in range(rows)])

    # needs a better name
    def matricize(self, row: int, brackets: tuple = ('\u23A2', '\u23A5')):
        row_str = str(self[row]).strip('[').strip(']').replace(',', '')
        line_str = f"\n {brackets[0]} {row_str} {brackets[1]}"
        return line_str

    def scalar_product(self, scalar):
        """Multiplies a matrix by a scalar"""

        return Matrix([list(map(lambda value: value * scalar, row))
                      for row in self])

    def dot_product(self, other):
        """Returns the dot product of two matrices"""

        matrix1 = self
        matrix2 = other

        matrix1_length = len(matrix1[0])
        matrix1_height = len(matrix1)

        matrix2_length = len(matrix2[0])
        matrix2_height = len(matrix2)

        if matrix1_length != matrix2_height:
            raise Exception("Incompatible matrices")  # needs more description

        output = []
        for row in range(matrix1_height):
            new_row = []
            for col in range(matrix2_length):
                total = 0
                for index in range(matrix1_length):
                    total += (matrix1[row][index]
                              * matrix2.transpose[col][index])
                new_row.append(total)
            output.append(new_row)
        return Matrix(output)

    def is_square(self):
        return self.rows == self.cols

    @property
    def rows(self):
        return len(self)

    @property
    def cols(self):
        return len(self[0])

    @property
    def trace(self):
        """Returns the sum of the values on the current matrix's diagonal"""
        sum(self.diag)

    @property
    def inverse(self):
        if self.determinant == 0:
            raise Exception('Matrix determinant is undefined')
        if not self.is_square():
            raise Exception('Only square matrices have inverses')

    @property
    def diag(self) -> list:
        """Returns the set of values in the current matrix
        where i = j
        """
        return [self[i][i] for i in range(self.rows)]

    @property
    def determinant(self):
        """Returns the determinant of the current matrix based off
        of its LU decomposition
        """
        if not self.is_square():
            raise Exception('Only square matrices have determinants')

        import operator
        import functools

        def prod(factors):
            return functools.reduce(operator.mul, factors, 1)

        return prod(self.LUdecomposition['upper'].diag)

    @property
    def LUdecomposition(self):
        """Uses Doolittle's Algorithm to find the decomposition
        of the current matrix A into a pair of upper and lower
        triangular matrices such that A = LU
        """
        lower = Matrix.generate_null(self.rows, self.cols)
        upper = Matrix.from_matrix(lower)
        # TODO add matrix .clone() functionality

        for i in range(self.rows):
            for k in range(i, self.rows):
                summation = sum([lower[i][j] * upper[j][k] for j in range(i)])
                upper[i][k] = self[i][k] - summation

            for k in range(i, self.rows):
                if i == k:
                    lower[i][i] = 1
                else:
                    summation = sum([lower[k][j] * upper[j][i]
                                    for j in range(i)])
                    lower[k][i] = (self[k][i] - summation) // upper[i][i]
        return {
            'lower': lower,
            'upper': upper
        }

    @property
    def transpose(self):
        """Returns the transpose of the current matrix"""

        return Matrix(map(list, zip(*self)))


class IdentityMatrix(Matrix):
    def __init__(self, rows, cols=None):

        if cols is None:
            cols = rows

        if rows != cols:
            raise Exception('Identity matrices must be square')

        array = [[1 if i == j else 0
                 for j in range(cols)]
                 for i in range(rows)]

        Matrix.__init__(self, array)

    @property
    def determinant(self):
        return 1

    @property
    def trace(self):
        return self.rows


class ZeroMatrix(Matrix):
    def __init__(self, rows, cols=None):

        if cols is None:
            cols = rows

        array = [[0 for j in range(cols)]
                 for i in range(rows)]

        Matrix.__init__(self, array)

    @property
    def trace(self):
        return 0


class PascalMatrix(Matrix):
    def __init__(self, rows, cols=None):

        if cols is None:
            cols = rows

        if rows != cols:
            raise Exception('Pascal matrices must be square')

        from math import factorial as fact

        def pascal(i, j): return fact(i + j) // (fact(i) * fact(j))

        array = [[pascal(i, j)
                 for j in range(cols)]
                 for i in range(rows)]

        Matrix.__init__(self, array)


if __name__ == "__main__":
    m1 = Matrix([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]
    ])

    def f(i, j):
        if 0 in (i, j) or i == j:
            return 1
        else:
            return 0

    m2 = Matrix.from_alg(f, 20)

    print(m2)
