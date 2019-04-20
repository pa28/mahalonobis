//
// Created by richard on 2019-04-19.
//

#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace matrix_math {

    template<typename T>
    class matrix {
    public:
        using size_type = typename std::vector<T>::size_type;

        size_type n_rows, n_cols;
        std::vector<std::vector<T>> data;

        matrix() = delete;

        matrix(size_type rows, size_type cols, T initial_value = 0)
                : n_rows{rows},
                  n_cols{cols},
                  data{rows} {
            for (size_type i = 0; i < n_rows; ++i) {
                for (size_type j = 0; j < n_cols; ++j) {
                    data[i].push_back(initial_value);
                }
            }
        }

        explicit matrix(std::vector<T> row)
                : n_rows{1},
                  n_cols{row.size()},
                  data{1} {
            for (auto value : row)
                data[0].push_back(value);
        }

        matrix(matrix &&m) noexcept : n_rows{m.n_rows}, n_cols{m.n_cols}, data{std::move(m.data)} {}

        matrix &operator=(matrix const &m) = default;

        matrix &operator=(matrix &&) noexcept = default;

        matrix(std::initializer_list<std::initializer_list<T>> &&initializerList) : data{} {
            for (auto ir : initializerList) {
                data.push_back(std::vector<T>{ir});
            }
            n_rows = data.size();
            n_cols = data[0].size();
        }

        auto begin() {
            return data.begin();
        }

        auto end() {
            return data.end();
        }

        auto cbegin() const {
            return data.cbegin();
        }

        auto cend() const {
            return data.cend();
        }

        std::ostream &print_on(std::ostream &os) const {
            for (size_type i = 0; i < n_rows; ++i) {
                for (size_type j = 0; j < n_cols; ++j) {
                    os << std::setw(5) << data[i][j];
                }
                os << '\n';
            }
            return os;
        }

        std::vector<T> &operator[](size_type i) { return data[i]; }

        std::vector<T> const &operator[](size_type i) const { return data[i]; }

        matrix operator*(T t) const {
            matrix result{n_rows, n_cols};
            for (size_type i = 0; i < n_rows; ++i)
                for (size_type j = 0; j < n_cols; ++j)
                    result.data[i][j] = t * data[i][j];
            return result;
        }

        matrix operator*(matrix const &m) const {
            if (n_cols != m.n_rows)
                throw std::logic_error("Matrix size error for multiplication.");

            matrix result{n_rows, m.n_cols};
            for (size_type i = 0; i < result.n_rows; ++i)
                for (size_type j = 0; j < result.n_cols; ++j) {
                    result.data[i][j] = 0;
                    for (size_type k = 0; k < n_cols; ++k) {
                        result.data[i][j] += data[i][k] * m.data[k][j];
                    }
                }

            return result;
        }

        matrix mean_data() const {
            matrix result{1, n_cols};

            for (size_type j = 0; j < result.n_cols; ++j) {
                for (size_type k = 0; k < n_rows; ++k) {
                    result.data[0][j] += data[k][j];
                }
                result.data[0][j] /= n_rows;
            }

            return result;
        }

        matrix normalize_data() const {
            auto u = mean_data();
            matrix result{n_rows, n_cols};
            for (size_type i = 0; i < result.n_rows; ++i)
                for (size_type j = 0; j < result.n_cols; ++j) {
                    result.data[i][j] = data[i][j] - u[0][j];
                }
            return result;
        }

        matrix operator/(T c) const {
            return operator*(1 / c);
        }

        matrix operator+(matrix const &m) const {
            if (n_cols != m.n_cols || n_rows != m.n_rows)
                throw std::logic_error("Matrix size error for addition.");

            matrix result{n_rows, n_cols};
            for (size_type i = 0; i < result.n_rows; ++i)
                for (size_type j = 0; j < result.n_cols; ++j) {
                    result.data[i][j] = data[i][j] + m.data[i][j];
                }
            return result;
        }

        matrix operator-(matrix const &m) const {
            if (n_cols != m.n_cols || n_rows != m.n_rows)
                throw std::logic_error("Matrix size error for subtraction.");

            matrix result{n_rows, n_cols};
            for (size_type i = 0; i < result.n_rows; ++i)
                for (size_type j = 0; j < result.n_cols; ++j) {
                    result.data[i][j] = data[i][j] - m.data[i][j];
                }
            return result;
        }

        matrix operator~() const {
            matrix result{n_cols, n_rows};
            for (size_type i = 0; i < n_rows; ++i)
                for (size_type j = 0; j < n_cols; ++j)
                    result.data[j][i] = data[i][j];
            return result;
        }

        matrix covariance() const {
            auto a = normalize_data();
            auto r = (~a * a) / 5;
            return r;
        }

        // Adapted from https://www.geeksforgeeks.org/adjoint-inverse-matrix/
        // Function to get cofactor of A[p][q] in temp[][]. n is current
        // dimension of A[][]
        matrix getCofactor(size_type p, size_type q, size_type n) {
            size_type i = 0, j = 0;

            matrix result{n, n};

            // Looping for each element of the matrix
            for (size_type row = 0; row < n; row++) {
                for (size_type col = 0; col < n; col++) {
                    //  Copying into temporary matrix only those element
                    //  which are not in given row and column
                    if (row != p && col != q) {
                        result.data[i][j++] = data[row][col];

                        // Row is filled, so increase row index and
                        // reset col index
                        if (j == n - 1) {
                            j = 0;
                            i++;
                        }
                    }
                }
            }

            return result;
        }

        /* Recursive function for finding determinant of matrix.
           n is current dimension of A[][]. */
        T determinant(size_type n) {
            T D = 0; // Initialize result

            //  Base case : if matrix contains single element
            if (n == 1)
                return data[0][0];

            int sign = 1;  // To store sign multiplier

            // Iterate for each element of first row
            for (size_type f = 0; f < n; f++) {
                // Getting Cofactor of A[0][f]
                auto temp = getCofactor(0, f, n);
                D += sign * data[0][f] * temp.determinant(n - 1);

                // terms are to be added with alternate sign
                sign = -sign;
            }

            return D;
        }

        // Function to get adjoint of A[N][N] in adj[N][N].
        matrix adjoint() {
            matrix adj{n_rows, n_cols};
            if (n_rows == 1) {
                adj[0][0] = 1;
                return adj;
            }

            int sign = 1;

            for (int i = 0; i < n_rows; i++) {
                for (int j = 0; j < n_rows; j++) {
                    // Get cofactor of A[i][j]
                    auto temp = getCofactor(i, j, n_rows);

                    // sign of adj[j][i] positive if sum of row
                    // and column indexes is even.
                    sign = ((i + j) % 2 == 0) ? 1 : -1;

                    // Interchanging rows and columns to get the
                    // transpose of the cofactor matrix
                    adj[j][i] = (sign) * (temp.determinant(n_rows - 1));
                }
            }
            return adj;
        }

        // Function to calculate and store inverse, returns false if
        // matrix is singular
        matrix inverse() {
            // Find determinant of this
            T det = determinant(n_rows);
            if (det == 0)
                throw std::domain_error("Singular matrix, can't find its inverse");

            // Find adjoint
            auto adj = adjoint();

            // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
            return adj / det;
        }

        T mahlaonabis(std::vector<T> const &delta) const {
            T result;
            matrix matrix_delta{delta};
            result = (matrix_delta * operator*(~matrix_delta))[0][0];
            return result;
        }
    };

    template<typename T>
    std::ostream &operator<<(std::ostream &ostream, matrix<T> const &m) {
        return m.print_on(ostream);
    }

} // namespace matrix

template<typename T>
std::vector<T> operator+(std::vector<T> const &v0, std::vector<T> const &v1) {
    if (v0.size() == v1.size()) {
        std::vector<T> r(v0.size());
        for (size_t i = 0; i < v0.size(); ++i) {
            r[i] = v0[i] + v1[i];
        }
        return r;
    } else
        throw std::domain_error("Vectors must be the same length for addition.");
}

template<typename T>
std::vector<T> operator-(std::vector<T> const &v0, std::vector<T> const &v1) {
    if (v0.size() == v1.size()) {
        std::vector<T> r(v0.size());
        for (size_t i = 0; i < v0.size(); ++i) {
            r[i] = v0[i] - v1[i];
        }
        return r;
    } else
        throw std::domain_error("Vectors must be the same length for addition.");
}

