#include "hw1.h"

#include <random>
#include <iomanip>

namespace algebra {
    /**
     * Creates a matrix of zeros with the specified dimensions.
     *
     * @param n the number of rows in the matrix
     * @param m the number of columns in the matrix
     *
     * @return a matrix of zeros with dimensions n x m
     *
     * @throws std::invalid_argument if n or m is less than or equal to 0
     */
    Matrix zeros(size_t n, size_t m)
    {
        if (n <= 0 || m <= 0) 
            throw std::invalid_argument("n and m should be greater than 0");
        return Matrix(n, std::vector<double>(m, 0.0));
    }

    /**
     * Generates a matrix filled with ones.
     *
     * @param n the number of rows in the matrix
     * @param m the number of columns in the matrix
     *
     * @return a matrix of size n x m with all elements set to 1.0
     *
     * @throws std::invalid_argument if n or m is less than or equal to 0
     */
    Matrix ones(size_t n, size_t m)
    {
        if (n <= 0 || m <= 0) 
            throw std::invalid_argument("n and m should be greater than 0");
        return Matrix(n, std::vector<double>(m, 1.0));
    }

    /**
     * Generates a random matrix of size n x m with values between min and max.
     *
     * @param n the number of rows in the matrix
     * @param m the number of columns in the matrix
     * @param min the minimum value of the random number
     * @param max the maximum value of the random number
     *
     * @return the randomly generated matrix
     *
     * @throws std::invalid_argument if n or m is less than or equal to 0
     * @throws std::logic_error if min is greater than max
     */
    Matrix random(size_t n, size_t m, double min, double max)
    {
        if (n <= 0 || m <= 0) 
            throw std::invalid_argument("n and m should be greater than 0");
        if (min > max) {
            throw std::logic_error("min cannot be greater than max");
        }
        Matrix matrix = zeros(n, m);

        std::default_random_engine engine;
        std::uniform_real_distribution<double> distribution(min, max);
        for (auto& row : matrix) {
            for (auto& elem : row) {
                elem = distribution(engine);
            }
        }
        return matrix;
    }

    /**
     * Displays the elements of a Matrix object in a formatted manner.
     *
     * @param matrix the Matrix object to be displayed
     *
     * @throws None
     */
    void show(const Matrix& matrix)
    {
        std::cout << std::noshowpoint << std::setprecision(3);
        for (auto& row : matrix) {
            for (auto& elem : row) {
                std::cout << elem << "\t";
            }
            std::cout << std::endl;
        }
    }

    /**
     * Multiply a matrix by a scalar value.
     *
     * @param matrix the matrix to be multiplied
     * @param c the scalar value to multiply the matrix by
     *
     * @return the resulting matrix after multiplication
     *
     * @throws None
     */
    Matrix multiply(const Matrix& matrix, double c)
    {
        int n = matrix.size();
        int m = matrix[0].size();
        Matrix mat = zeros(n, m);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                mat[i][j] = matrix[i][j] * c;
            }
        }
        return mat;
    }

    /**
     * Multiplies two matrices.
     *
     * @param matrix1 the first matrix to be multiplied
     * @param matrix2 the second matrix to be multiplied
     *
     * @return the resulting matrix after multiplication
     *
     * @throws std::logic_error if the dimensions of the matrices are incompatible for multiplication
     */
    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
    {
        if (matrix1.size() == 0 || matrix2.size() == 0) return {};

        if (matrix1[0].size() != matrix2.size()) {
            throw std::logic_error("The matrix with wrong dimensions cannot be multiplied");
        }
        Matrix mat = zeros(matrix1.size(), matrix2[0].size());
        for (int i=0; i<mat.size(); i++) {
            for (int j=0; j<mat[i].size(); j++) {
                for (int k=0; k<matrix1[0].size(); k++) {
                    mat[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return mat;
    }

    /**
     * Calculates the sum of a matrix and a constant value.
     *
     * @param matrix The matrix to be summed.
     * @param c The constant value to be added to each element of the matrix.
     *
     * @return The resulting matrix after the sum operation.
     *
     * @throws None
     */
    Matrix sum(const Matrix& matrix, double c)
    {
        if (matrix.size() == 0) return {};

        int n = matrix.size();
        int m = matrix[0].size();
        Matrix mat = zeros(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                mat[i][j] = matrix[i][j] + c;
            }
        }
        return mat;
    }

    /**
     * Sums two matrices element-wise and returns the result.
     *
     * @param matrix1 the first matrix to be summed
     * @param matrix2 the second matrix to be summed
     *
     * @return the resulting matrix after element-wise summation
     *
     * @throws std::logic_error if the matrices have different dimensions
     */
    Matrix sum(const Matrix& matrix1, const Matrix& matrix2)
    {
        if (matrix1.size() == 0 && matrix2.size() == 0) return {};

        if (matrix1.size() == 0 || matrix2.size() == 0) {
            throw std::logic_error("The matrix with wrong dimensions cannot be summed");
        }

        if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
            throw std::logic_error("The matrix with wrong dimensions cannot be summed");
        }

        Matrix mat = zeros(matrix1.size(), matrix1[0].size());
        for (int i=0; i<mat.size(); i++) {
            for (int j=0; j<mat[i].size(); j++) {
                mat[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return mat;
    }

    /**
     * Transposes a given matrix.
     *
     * @param matrix the matrix to be transposed
     *
     * @return the transposed matrix
     *
     * @throws None
     */
    Matrix transpose(const Matrix& matrix)
    {
        if (matrix.size() == 0) return {};

        int n = matrix.size();
        int m = matrix[0].size();
        Matrix mat = zeros(m, n);

        for (int i=0; i<n; i++) {
            for (int j=0; j<m; j++) {
                mat[j][i] = matrix[i][j];
            }
        }
        return mat;
    }

    /**
     * Generates the minor matrix of the given matrix at the specified position.
     *
     * @param matrix The input matrix.
     * @param n The row index of the element to exclude.
     * @param m The column index of the element to exclude.
     *
     * @return The minor matrix.
     *
     * @throws std::logic_error if the matrix dimensions are incorrect.
     */
    Matrix minor(const Matrix& matrix, size_t n, size_t m)
    {
        if (matrix.size() == 0) return {};

        if (matrix.size() <= n || matrix[0].size() <= m) {
            throw std::logic_error("The matrix with wrong dimensions cannot be minor");
        }
        int skipRow = n, skipCol = m;
        Matrix mat = zeros(matrix.size()-1, matrix[0].size()-1);
        for (int i=0; i<mat.size(); i++) {
            for (int j=0; j<mat[i].size(); j++) {
                int ii = (i<skipRow) ? i : i+1;
                int jj = (j<skipCol) ? j : j+1;
                mat[i][j] = matrix[ii][jj];
            }
        }
        return mat;
    }
    /**
     * Calculates the determinant of a square matrix.
     *
     * @param matrix the input matrix
     *
     * @return the determinant of the matrix
     *
     * @throws std::logic_error if the matrix has wrong dimensions
     */
    double determinant(const Matrix& matrix)
    {
        int n = matrix.size();
        if (n == 0) return 1;

        if (matrix.size() != matrix[0].size()) {
            throw std::logic_error("The matrix with wrong dimensions cannot be determinant");
        }
        if (n == 1) return matrix[0][0];
        if (n == 2) return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
        if (n == 3) {
            return  matrix[0][0] * matrix[1][1] * matrix[2][2] +
                    matrix[0][1] * matrix[1][2] * matrix[2][0] +
                    matrix[0][2] * matrix[1][0] * matrix[2][1] -
                    matrix[0][2] * matrix[1][1] * matrix[2][0] -
                    matrix[0][0] * matrix[1][2] * matrix[2][1] -
                    matrix[0][1] * matrix[1][0] * matrix[2][2];
        }
        double ret = 0;
        for (int i=0; i<n; i++) {
           Matrix m = minor(matrix, 0, i);
           int sign = pow(-1, 1+(i+1));
           ret += matrix[0][i] * sign *  determinant(m);
        }
        return ret;
    }
    /**
     * Generates the adjoint matrix of the given matrix.
     *
     * @param matrix the input matrix
     *
     * @return the adjoint matrix of the input matrix
     *
     * @throws None
     */
    static Matrix adj(const Matrix& matrix) {
        int n = matrix.size();
        Matrix mat = zeros(n, n);
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                mat[i][j] = pow(-1, i+j+2) * determinant(minor(matrix, j, i));
            }
        }
        return mat;
    }
    /**
     * Calculates the inverse of a square matrix.
     *
     * @param matrix the input matrix
     *
     * @return the inverse of the input matrix
     *
     * @throws std::logic_error if the input matrix has wrong dimensions
     * @throws std::logic_error if the input matrix is not invertible
     */
    Matrix inverse(const Matrix& matrix)
    {
        int n = matrix.size();
        if (n == 0) return matrix;
        if (matrix.size() != matrix[0].size()) {
            throw std::logic_error("The matrix with wrong dimensions cannot be determinant");
        }
        double det = determinant(matrix);
        if (det == 0) {
            throw std::logic_error("The matrix cannot be inverse");
        }
        Matrix mat = adj(matrix);
        return multiply(mat, 1/det);
    }
    /**
     * Concatenates two matrices along a specified axis.
     *
     * @param matrix1 The first matrix to concatenate.
     * @param matrix2 The second matrix to concatenate.
     * @param axis The axis along which to concatenate the matrices. Default is 0.
     *
     * @return The concatenated matrix.
     *
     * @throws std::logic_error if the axis is not 0 or 1, or if the matrices have incompatible dimensions.
     */
    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0)
    {
        if (axis != 0 && axis != 1) {
            throw std::logic_error("The axis should be 0 or 1");
        }
        int n1 = matrix1.size();
        int m1 = matrix1[0].size();
        int n2 = matrix2.size();
        int m2 = matrix2[0].size();
        if (axis == 0) {
            if (m1 != m2) {
                throw std::logic_error("The matrix with wrong dimensions cannot be concatenated");
            }
            Matrix mat = zeros(n1+n2, m1);
            for (int i=0; i<n1; i++) {
                for (int j=0; j<m1; j++) {
                    mat[i][j] = matrix1[i][j];
                }
            }
            for (int i=0; i<n2; i++) {
                for (int j=0; j<m2; j++) {
                    mat[i+n1][j] = matrix2[i][j];
                }
            }
            return mat;
        } else {
            if (n1 != n2) {
                throw std::logic_error("The matrix with wrong dimensions cannot be concatenated");
            }
            Matrix mat = zeros(n1, m1+m2);
            for (int i=0; i<n1; i++) {
                for (int j=0; j<m1; j++) {
                    mat[i][j] = matrix1[i][j];
                }
                for (int j=0; j<m2; j++) {
                    mat[i][m1+j] = matrix2[i][j];
                }
            }
            return mat;
        }
    }

    Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2)
    {
        return {};
    }

    Matrix ero_multiply(const Matrix& matrix, size_t r, double c)
    {
        return {};
    }

    Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2)
    {
        return {};
    }

    Matrix upper_triangular(const Matrix& matrix)
    {
        return {};
    }
}