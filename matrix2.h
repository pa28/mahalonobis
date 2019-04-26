//
// Created by richard on 2019-04-25.
//

#pragma once

#include <cstring>

namespace matrix_math {

    template<typename T>
    class matrix;

    template<typename T>
    class vector {
    public:
        using value_type = T;
        using size_type = size_t;
        using iterator = T *;
        using const_iterator = const T *;

        friend class matrix<T>;

    private:
        size_t size{0};
        T *data{nullptr};

    public:
        explicit vector(size_t s)
                : size{s} {
            data = new T[size];
            for (size_t idx = 0; idx < size; ++idx)
                data[idx] = 0;
        }

        ~vector() {
            delete[] data;
        }

        vector(vector const &other)
                : size{other.size} {
            data = new T[size];
            if (other.data)
                memcpy(data, other.data, size * sizeof(T));
        }

        vector(vector &&other) noexcept
                : size{other.size} {
            data = other.data;
            other.data = nullptr;
            other.size = 0;
        }

        vector(std::initializer_list<T> initializerList)
                : size{initializerList.size()} {
            data = new T[initializerList.size()];
            std::copy(initializerList.begin(), initializerList.end(), begin());
        }

        vector &operator=(vector const &other) {
            if (size != other.size) {
                delete[] data;
                size = other.size;
                data = new T[size];
            }
            memcpy(data, other.data, size * sizeof(T));
            return *this;
        }

        vector &operator=(vector &&other) noexcept {
            delete[] data;
            size = other.size;
            data = other.data;
            other.size = 0;
            other.data = nullptr;
            return *this;
        }

        T &operator[](size_type idx) {
            return data[idx];
        }

        const T &operator[](size_type idx) const {
            return data[idx];
        }

        iterator begin() { return data; }

        iterator end() { return data + size; }

        const iterator begin() const { return data; }

        const iterator end() const { return data + size; }

        const_iterator cbegin() const { return data; }

        const_iterator cend() const { return data + size; }

        vector operator*(T x) const {
            vector v(size);
            auto a = v.begin();
            auto b = begin();
            while (a != v.end())
                *a++ = *b++ * x;
            return v;
        }

        std::ostream &print_on(std::ostream &ostream) const {
            for (auto v : *this)
                ostream << v << ' ';
            return ostream;
        }
    };

    template<typename T>
    vector<T> operator*(T x, vector<T> const &v) {
        return v.operator*(x);
    }

    template<typename T>
    class matrix {
    public:
        using value_type = T;
        using vector_type = vector<T>;
        using size_type = size_t;

        class iterator {
        private:
            friend class matrix;
            matrix &m;
            size_t offset{0};

            explicit iterator(matrix &m, size_t offset = 0) :m{m}, offset{offset} {}

        public:
            iterator() = delete;

            iterator& operator++() { ++offset; return *this; }

            const iterator& operator++() const { ++offset; return *this; }

            bool operator!=(const iterator &other) const { return offset != other.offset; }

            vector<T>& operator*() { return *m.data[offset]; }

            const vector<T>& operator*() const { return *m.data[offset]; }

            vector<T>* operator->() { return m.data[offset]; }

            const vector<T>* operator->() const { return m.data[offset]; }
        };

    private:
        size_t rows{0}, cols{0};
        vector<T> **data{nullptr};

    public:
        matrix(size_t r, size_t c)
                : rows{r},
                  cols{c} {
            data = new vector<T>*[rows];
            for (size_t idx = 0; idx < rows; ++idx)
                data[idx] = new vector<T>(cols);
        }

        ~matrix() {
            for (size_t r = 0; r < rows; ++r)
                delete data[r];
            delete[] data;
        }

        matrix(std::initializer_list<vector<T>> initializerList)
                : rows{initializerList.size()} {
            data = new vector<T>*[initializerList.size()];
            auto il = initializerList.begin();
            for (size_t idx = 0; idx < rows; ++idx, ++il) {
                data[idx] = new vector<T>(il->size);
                *data[idx] = *il;
            }
        }

        iterator begin() { return iterator(*this); }

        iterator end() { return iterator(*this,rows); }

        const iterator begin() const { return iterator(const_cast<matrix&>(*this)); }

        const iterator end() const { return iterator(const_cast<matrix&>(*this),rows); }

        const iterator cbegin() const { return iterator(const_cast<matrix&>(*this)); }

        const iterator cend() const { return iterator(const_cast<matrix&>(*this),rows); }

        vector<T>&operator[](size_t idx) { return *data[idx]; }

        const vector<T>&operator[](size_t idx) const { return *data[idx]; }

        std::ostream &print_on(std::ostream &ostream) const {
            for (auto idx = begin(); idx != end(); ++idx)
                idx->print_on(ostream) << std::endl;
            return ostream;
        }
    };
}
