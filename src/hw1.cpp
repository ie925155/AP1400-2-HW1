#include "hw1.h"

#include <random>

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

    void show(const Matrix& matrix)
    {

    }

    Matrix multiply(const Matrix& matrix, double c)
    {
        return {};
    }

    Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
    {
        return {};
    }

    Matrix sum(const Matrix& matrix, double c)
    {
        return {};
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
