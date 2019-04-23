//
// Created by richard on 2019-04-20.
//

#pragma once

#include <vector>
#include <iostream>

namespace mahalonobis {

    template<typename T>
    class matrix;

    template<typename T>
    class vector {
    public:
        using size_type = size_t;
        using iterator = T *;
        using const_iterator = const T *;

    protected:
        T *data{nullptr};
        size_type size{0};
        bool transposed{false};

        friend class matrix<T>;

    public:

        vector() = default;

        vector(vector const &other) = default;

        vector(vector &&other) noexcept = default;

        vector& operator=(vector const &other) = default;

        vector& operator=(vector &&other) noexcept = default;

        ~vector() = default;

        explicit vector(size_type n)
                : size{n},
                  data{new T[n]} {
            std::fill(begin(), end(), 0);
        }

        explicit vector(std::vector<T> const &v)
                : size{v.size()},
                  data{new T[v.size()]} {
            std::copy(v.cbegin(), v.cend(), begin());
        }

        vector(std::initializer_list<T> initializerList)
                : size{initializerList.size()},
                  data{new T[initializerList.size()]} {
            std::copy(initializerList.begin(), initializerList.end(), begin());
        }

        iterator begin() { return data; }

        iterator end() { return data + size; }

        const iterator begin() const { return data; }

        const iterator end() const { return data + size; }

        const_iterator cbegin() const { return data; }

        const_iterator cend() const { return data + size; }

        T &operator[](size_type idx) {
            return data[idx];
        }

        const T &operator[](size_type idx) const {
            return data[idx];
        }

        vector operator+(vector const o) const {
            if (size == o.size && transposed == o.transposed) {
                vector r(o.size);
                const_iterator ti = cbegin();
                const_iterator oi = o.cbegin();
                iterator ri = r.begin();
                while (ri != r.end()) {
                    *ri++ = *ti++ + *oi++;
                }
                return r;
            } else
                throw std::logic_error("vectors are not congruent for addition.");
        }

        vector operator-(vector const o) const {
            if (size == o.size && transposed == o.transposed) {
                vector r(o.size);
                const_iterator ti = cbegin();
                const_iterator oi = o.cbegin();
                iterator ri = r.begin();
                while (ri != r.end()) {
                    *ri++ = *ti++ - *oi++;
                }
                return r;
            } else
                throw std::logic_error("vectors are not congruent for subtraction.");
        }

        T operator*(vector const o) const {
            if (size == o.size) {
                const_iterator ti = cbegin();
                const_iterator oi = o.cbegin();
                T sum = 0;
                while (ti != cend())
                    sum += *ti++ * *oi++;
                return sum;
            } else
                throw std::logic_error("vectors are not congruent for dot product.");
        }

        std::ostream &print_on(std::ostream &ostream) const {
            for (auto d : *this)
                ostream << d << ' ';
            return ostream;
        }
    };

    template<typename T>
    class matrix {
    public:
        using size_type = size_t;
        using iterator = vector<T> *;
        using const_iterator = const vector<T> *;

    protected:
        size_type rows{0};
        vector<T> *data{nullptr};

    public:

        matrix() = default;
        
        matrix(matrix const &other) = default;

        matrix(matrix &&other) noexcept = default;

        matrix& operator=(matrix const &other) = default;

        matrix& operator=(matrix &&other) noexcept = default;

        ~matrix() = default;

        matrix(std::initializer_list<vector<T>> initializerList)
                : rows{initializerList.size()},
                  data{new vector<T>[initializerList.size()]} {
            std::copy(initializerList.begin(), initializerList.end(), begin());
        }

        matrix(size_type rows, size_type cols)
                : rows{rows},
                  data{new vector<T>[rows]} {
            for (size_type idx = 0; idx < rows; ++idx)
                data[idx] = vector<T>(cols);
        }

        /**
         * @brief Begin iterator.
         */
        iterator begin() { return data; }

        /**
         * @brief End iterator.
         * @return
         */
        iterator end() { return data + rows; }

        /**
         * @brief Constant begin iterator.
         */
        const iterator begin() const { return data; }

        /**
         * @brief Constant end iterator.
         */
        const iterator end() const { return data + rows; }

        /**
         * @brief Constant cbegin iterator.
         */
        const_iterator cbegin() const { return data; }

        /**
         * @brief Constant cend iterator.
         */
        const_iterator cend() const { return data + rows; }

        /**
         * @brief Indexing operator, returns a selected row.
         * @param idx The index of the row requested
         * @return A vector object containing the requested row.
         */
        vector<T> &operator[](size_type idx) {
            return data[idx];
        }

        /**
         * @brief Constant indexing operator, returns a selected row.
         * @param idx The index of the row requested
         * @return A constant vector object containing the requested row.
         */
        const vector<T> &operator[](size_type idx) const {
            return data[idx];
        }

        /**
         * @brief Return the result of multiplying this matix by m.
         * @param m the multiplier.
         * @return The matrix product.
         */
        matrix operator*(matrix const &m) const {
            if (data[0].size != m.rows)
                throw std::logic_error("Matrix size error for multiplication.");

            matrix result{rows, m.data[0].size};
            for (size_type i = 0; i < result.rows; ++i)
                for (size_type j = 0; j < result.data[i].size; ++j) {
                    result.data[i][j] = 0;
                    for (size_type k = 0; k < m.rows; ++k) {
                        result.data[i][j] += data[i][k] * m.data[k][j];
                    }
                }

            return result;
        }

        vector<T> operator*(vector<T> const &v) const {
            if (data[0].size != v.size)
                throw std::logic_error("Matrix/vector size error for multiplication.");

            vector<T> result{rows};
            for (size_type i = 0; i < rows; ++i)
                result[i] = data[i] * v;
            return result;
        }

        /**
         * @brief Multiply this matrix by a scalar.
         * @param t the scalar multiplier.
         * @return the product.
         */
        matrix operator*(T t) const {
            matrix result(rows, data[0].size);
            for (size_type i = 0; i < rows; ++i)
                for (size_type j = 0; j < data[i].size; ++j)
                    result[i][j] = data[i][j] * t;
            return result;
        }

        /**
         * @brief Divide this matrix by a scalar.
         * @param t the scalar divisor.
         * @return the result.
         */
        matrix operator/(T t) const {
            matrix result(rows, data[0].size);
            for (size_type i = 0; i < rows; ++i)
                for (size_type j = 0; j < data[i].size; ++j)
                    result[i][j] = data[i][j] / t;
            return result;
        }

        /**
         * @brief Compute the mean of the data.
         * @return A vector mean
         * @details Treats the matrix as a set of multivariate values, one per row of the matrix, returns a
         * vector with the mean value.
         */
        vector<T> mean_data() const {
            vector<T> result(data[0].size);

            for (auto v : *this) {
                for (size_type j = 0; j < v.size; ++j) {
                    result[j] += v[j];
                }
            }

            for (auto &r : result)
                r /= rows;

            return result;
        }

        /**
         * @brief Normalize the matrix
         * @return the normalized data
         * @details Treats the matrix as a set of multivariate values, one per row of the matrix, returns a
         * matrix with the mean subtracted from each value.
         */
        matrix normalize_data() const {
            auto u = mean_data();
            matrix result(rows, data[0].size);
            auto ri = result.begin();
            for (auto &v : *this) {
                *ri++ = v - u;
            }
            return result;
        }

        /**
         * @brief Transpose the matrix.
         */
        matrix operator~() const {
            matrix result(data[0].size, rows);
            for (size_type i = 0; i < rows; ++i)
                for (size_type j = 0; j < data[0].size; ++j)
                    result[j][i] = data[i][j];
            return result;
        }

        /**
         * @brief Compute the covariance matrix of this matrix.
         * @return
         */
        matrix covariance() const {
            auto a = normalize_data();
            auto r = (~a * a) / rows;
            return r;
        }

        // Adapted from https://www.geeksforgeeks.org/adjoint-inverse-matrix/
        // Function to return the cofactor of A[p][q]. n is current effective
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

        // Function to return adjoint of A[N][N].
        matrix adjoint() {
            matrix adj{rows, rows};
            if (rows == 1) {
                adj[0][0] = 1;
                return adj;
            }

            int sign = 1;

            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < rows; j++) {
                    // Get cofactor of A[i][j]
                    auto temp = getCofactor(i, j, rows);

                    // sign of adj[j][i] positive if sum of row
                    // and column indexes is even.
                    sign = ((i + j) % 2 == 0) ? 1 : -1;

                    // Interchanging rows and columns to get the
                    // transpose of the cofactor matrix
                    adj[j][i] = (sign) * (temp.determinant(rows - 1));
                }
            }
            return adj;
        }

        // Function to calculate and return the inverse, returns false if
        // matrix is singular
        matrix inverse() {
            // Find determinant of this
            T det = determinant(rows);
            if (det == 0)
                throw std::domain_error("Singular matrix, can't find its inverse");

            // Find adjoint
            auto adj = adjoint();

            // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
            return adj / det;
        }

        std::ostream &print_on(std::ostream &ostream) const {
            for (auto v : *this)
                ostream << v << '\n';
            return ostream;
        }

        /**
         * @brief Compute the mahalanobis distance for vector delta
         * @param delta
         * @return The distance
         * @details Delta can be the result of the difference between two multivariate values,
         * or between a multivariate value and the mean value.
         */
        T mahalanobis(vector<T> const &delta) const {
            T result;
            matrix matrix_delta{delta};
            result = (matrix_delta * operator*(~matrix_delta))[0][0];
            return result;
        }
    };

    /**
     * @brief Allow the use of the comutative property of multiplication.
     * @tparam T
     * @param t a scalar multiplier
     * @param m a matrix
     * @return the result of m * t;
     */
    template<typename T>
    matrix<T> operator*(T t, matrix<T> const &m) {
        return m.operator*(t);
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &ostream, vector<T> const &v) {
        return v.print_on(ostream);
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &ostream, matrix<T> const &m) {
        return m.print_on(ostream);
    }
}