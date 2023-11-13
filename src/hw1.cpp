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

    Matrix sum(const Matrix& matrix1, const Matrix& matrix2)
    {
        return {};
    }

    Matrix transpose(const Matrix& matrix)
    {
        return {};
    }

    Matrix minor(const Matrix& matrix)
    {
        return {};
    }

    double determinant(const Matrix& matrix)
    {
        return {};
    }

    Matrix inverse(const Matrix& matrix)
    {
        return {};
    }

    Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0)
    {
        return {};
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
