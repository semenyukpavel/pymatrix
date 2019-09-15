class Matrix:
    def __init__(self, nrows, ncols):
        self.n = nrows
        self.m = ncols
        self.rows = [[0 for _ in range(ncols)] for _ in range(nrows)]

    def __str__(self):
        s = str()
        for i in range(self.n):
            s += '\n'
            s += ' '.join(str(item) for item in self.rows[i])
        s += '\n'
        return s

    def __getitem__(self, index):
        return self.rows[index]

    @classmethod
    def create_id(cls, m):
        res = Matrix(m, m)
        for i in range(m):
            res[i][i] = 1
        return res

    @classmethod
    def create_diag(cls, arr):
        res = Matrix(len(arr), len(arr))
        for i in range(len(arr)):
            res[i][i] = arr[i]
        return res

    @classmethod
    def create_values(cls, vals):
        n = len(vals)
        m = len(vals[0])
        res = Matrix(n, m)
        for i in range(n):
            for j in range(m):
                res[i][j] = vals[i][j]
        return res

    @classmethod
    def create_vector_column(cls, val):
        n = len(val)
        res = Matrix(n, 1)
        for i in range(n):
            res[i][0] = val[i]
        return res

    @classmethod
    def create_vector_row(cls, val):
        m = len(val)
        res = Matrix(1, m)
        for i in range(m):
            res[0][i] = val[i]
        return res

    def __add__(self, right):
        if self.n != right.n or self.m != right.m:
            print('Cannot add two matrices with different sizes')
            exit(1)
        res = Matrix(self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = self[i][j] + right[i][j]
        return res

    def __sub__(self, right):
        if self.n != right.n or self.m != right.m:
            print('Cannot subtract two matrices with different sizes')
            exit(1)
        res = Matrix(self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = self[i][j] - right[i][j]
        return res

    def __eq__(self, right):
        if self.n != right.n or self.m != right.m:
            return False
        for i in range(self.n):
            for j in range(self.m):
                if self[i][j] != right[i][j]:
                    return False
        return True

    def __iadd__(self, right):
        if self.n != right.n or self.m != right.m:
            print('Cannot add two matrices with different sizes')
            exit(1)
        for i in range(self.n):
            for j in range(self.m):
                self[i][j] += right[i][j]
        return self

    def __isub__(self, right):
        if self.n != right.n or self.m != right.m:
            print('Cannot subtract two matrices with different sizes')
            exit(1)
        for i in range(self.n):
            for j in range(self.m):
                self[i][j] -= right[i][j]
        return self

    def __mul__(self, right):
        if self.m != right.n:
            print('Cannot multiply matrices with given sizes')
            exit(2)
        res = Matrix(self.n, right.m)
        for i in range(self.n):
            for j in range(right.m):
                for k in range(len(self[i])):
                    res[i][j] += self[i][k]*right[k][j]
        return res

    def __imul__(self, right):
        if self.m != right.n:
            print('Cannot multiply matrices with given sizes')
            exit(2)
        res = Matrix(self.n, right.m)
        for i in range(self.n):
            for j in range(right.m):
                for k in range(len(self[i])):
                    res[i][j] += self[i][k]*right[k][j]
        for i in range(self.n):
            for j in range(self.m):
                self[i][j] = res[i][j]
        return self

    def mult_by_scalar(self, scal):
        res = Matrix(self.n, self.m)
        for i in range(self.n):
            for j in range(self.m):
                res[i][j] = scal*self[i][j]
        return res

    def transpose(self):
        res = Matrix(self.m, self.n)
        for i in range(self.m):
            for j in range(self.n):
                res[i][j] = self[j][i]
        return res

    def trace(self):
        if self.n != self.m:
            print('Trace can be evaluated only for square matrices')
            exit(3)
        res = 0
        for i in range(self.n):
            res += self[i][i]
        return res

    def is_identity(self):
        if self.n != self.m:
            return False
        identity = Matrix.create_id(self.n)
        return identity == self

    def is_diagonal(self):
        if self.n != self.m:
            return False
        for i in range(self.n):
            for j in range(self.n):
                if i != j and self[i][j] != 0:
                    return False
        return True

    def is_null(self):
        for i in range(self.n):
            for j in range(self.m):
                if self[i][j] != 0:
                    return False
        return True

    def transform_1(self, i, j):
        for k in range(len(self[i])):
            self[i][k], self[j][k] = self[j][k], self[i][k]

    def transform_2(self, i, j, scal):
        for k in range(len(self[i])):
            self[j][k] += scal*self[i][k]

    def transform_3(self, i, scal):
        for k in range(len(self[i])):
            self[i][k] *= scal

    def gauss_form(self):
        if self.n != self.m:
            print('Determinant can be evaluated only for square matrices')
            exit(4)
        for i in range(self.n-1):
            scal = self[i][i]
            for j in range(i+1, self.n):
                self.transform_2(i, j, -self[j][i]/scal)

    def det(self):
        matrix = Matrix.create_values(self.rows)
        matrix.gauss_form()
        res = 1
        for i in range(matrix.n):
            res *= matrix[i][i]
        return res

    def rank(self):
        matrix = Matrix.create_values(self.rows)
        matrix.gauss_form()
        rank = matrix.n
        for i in range(matrix.n-1, -1, -1):
            if matrix[i][i] == 0:
                rank -= 1
            else:
                return rank
        return 0

    def inverse(self):
        matrix = Matrix.create_values(self.rows)
        if matrix.n != matrix.m:
            print('Inverse matrix exists only for square matrices')
            exit(5)
        if matrix.det() == 0:
            print('Inverse matrix does not exist for matrix with det = 0')
            exit(6)
        inv = Matrix.create_id(matrix.n)
        for i in range(matrix.n-1):
            scal = matrix[i][i]
            for j in range(i+1, matrix.n):
                coef = -matrix[j][i] / scal
                matrix.transform_2(i, j, coef)
                inv.transform_2(i, j, coef)
        for i in range(matrix.n-1, 0, -1):
            scal = matrix[i][i]
            for j in range(i-1, -1, -1):
                coef = -matrix[j][i] / scal
                matrix.transform_2(i, j, coef)
                inv.transform_2(i, j, coef)
        for i in range(matrix.n):
            scal = matrix[i][i]
            matrix.transform_3(i, 1/scal)
            inv.transform_3(i, 1/scal)
        return inv

    def minor(self, rows, columns):
        vals = []
        for i in rows:
            row = []
            for j in columns:
                row.append(self[i][j])
            vals.append(row)
        m = Matrix.create_values(vals)
        return m


if __name__ == '__main__':
    vals_a = [[1/(i+j+2) for j in range(3)] for i in range(3)]
    a = Matrix.create_values(vals_a)

    print(a)

    a.gauss_form()
    print(a)
