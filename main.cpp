#include <iostream>
#include "matrix.h"
#include "mahalonobis.h"

using mahalonobis::vector;
using mahalonobis::matrix;

int main() {
    vector mv1{1, 2, 3};
    vector mv2{4, 5, 6};
    std::cout << mv1 << '\n'
              << mv2 << '\n'
              << mv1 + mv2 << '\n'
              << mv2 - mv1 << '\n'
              << mv1 * mv2 << '\n';

    matrix mm1{vector{1, 2, 3},
               vector{4, 5, 6}};

    std::cout << mm1 << '\n';

    matrix<long double> mm2(3, 3);

    mm2[0][0] = 5;

    std::cout << mm2 << '\n';

    matrix<long double> A{vector<long double>{90, 60, 90},
                    vector<long double>{90, 90, 30},
                    vector<long double>{60, 60, 60},
                    vector<long double>{60, 60, 90},
                    vector<long double>{30, 30, 30}};

    auto a = A.normalize_data();
    std::cout << A << '\n'
              << A.mean_data() << '\n'
              << a << '\n'
              << ~a << '\n'
              << (~a * a) / 5 << '\n';

    auto S = A.normalize_data().covariance();
    auto Si = S.inverse();

    std::cout << S << '\n'
              << Si << '\n'
              << S * Si << '\n';

/*
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
*/

    return 0;
}