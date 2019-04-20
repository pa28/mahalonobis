#include <iostream>
#include "matrix.h"

using matrix_math::matrix;

int main() {
    matrix<float> A{{90, 60, 90},
                    {90, 90, 30},
                    {60, 60, 60},
                    {60, 60, 90},
                    {30, 30, 30}};
    matrix<float> one{5, 5, 1};

    std::cout << A << '\n'
              << one << '\n';

    auto b = one * A / 5;
    auto a = A - b;
    matrix u{b[0]};
    std::cout << a << '\n'
              << A.normalize_data() << '\n'
              << A.mean_data() << '\n';

    auto S = A.normalize_data().covariance();
    auto Si = S.inverse();

    std::cout << S << '\n'
              << Si << '\n'
              << Si * S << '\n';

    for ( auto const &d : A)
        std::cout << Si.mahlaonabis(d - A[2]) << '\n';

    return 0;
}