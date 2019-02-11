#include <check.h>
#include <stdlib.h>

#include "vector.h"
#include "matrix.h"
#include "helpers.h"
#include "model.h"
#include "simplex.h"

/*
 * vector tests
 */
START_TEST(test_vector_new)
{
    vector *v = vector_new(3);
    ck_assert_ptr_nonnull(v);
    ck_assert_ptr_nonnull(v->data);
    ck_assert_uint_eq(v->size, 3);
    vector_free(v);
}
END_TEST

START_TEST(test_vector_elements)
{
    double tol = 1e-15;
    vector *v = vector_new(3);
    v->data[0] = .3;
    v->data[1] = .5;
    v->data[2] = .2;

    ck_assert_double_eq_tol(v->data[0], .3, tol);
    ck_assert_double_eq_tol(v->data[1], .5, tol);
    ck_assert_double_eq_tol(v->data[2], .2, tol);

    vector_free(v);
}
END_TEST

/*
 * matrix tests
 */
START_TEST(test_matrix_new)
{
    matrix *m = matrix_new(3, 2);
    ck_assert_ptr_nonnull(m);
    ck_assert_ptr_nonnull(m->data);
    ck_assert_uint_eq(m->nrow, 3);
    ck_assert_uint_eq(m->ncol, 2);
    matrix_free(m);
}
END_TEST

START_TEST(test_matrix_elements)
{
    double tol = 1e-15;
    size_t nrow = 3;
    size_t ncol = 2;

    matrix *m = matrix_new(nrow, ncol);
    double (*mm)[ncol] = (double (*)[ncol]) m->data;
    mm[0][1] = .3;

    ck_assert_double_eq_tol(mm[0][1], .3, tol);

    matrix_free(m);
}
END_TEST

/*
 * simplex tests
 */
START_TEST(test_simplex_transform)
{
    double tol = 1e-15;
    const size_t K=3;
    double th[K];
    double th_bar[K];
    double y[K-1];

    for(size_t k=0; k<K; k++)
        th[k] = 1. / K;

    transform(K, th, y);
    inverse_transform(K, y, th_bar);

    for(size_t k=0; k<K-1; k++)
        ck_assert_double_eq_tol(y[k], 0., tol);

    for(size_t k=0; k<K; k++)
        ck_assert_double_eq_tol(th[k], th_bar[k], tol);
}
END_TEST


START_TEST(test_sum)
{
    double tol = 1e-15;
    double x[3] = {1.1, 2.1, 0.07};
    double s = sum(x, 3);
    ck_assert_double_eq_tol(s, 3.27, tol);
}
END_TEST

START_TEST(test_log_sum_exp)
{
    double tol = 1e-15;
    double x[3] = { 1.0, 2.0, 3.0 };
    double expected = log(exp(1.0) + exp(2.0) + exp(3.0));
    double f = log_sum_exp(x, 3);
    ck_assert_double_eq_tol(expected, f, tol);
}
END_TEST

START_TEST(test_is_in)
{
    size_t x[2] = { 1, 2 };
    ck_assert_uint_eq(1, is_in(1, x, 2));
    ck_assert_uint_eq(0, is_in(3, x, 2));
}
END_TEST

START_TEST(test_loglik)
{
    double tol = 1e-15;
    double theta[3] = { 0.5, 0.3, 0.2 };
    size_t x[4][2] = { {0, 1}, {0, 1}, {0, 1}, {0, 2} };
    size_t y[4] = { 0, 0, 0, 2 };

    double ll_expected = 3 * (log(0.5) - log(0.5 + 0.3)) + (log(0.2) - log(0.5 + 0.2));
    double ll = loglik(4, 2, 3, theta, x, y);
    ck_assert_double_eq_tol(ll_expected, ll, tol);
}
END_TEST

START_TEST(test_full_conditionals_last_player_didnt_play)
{
    double tol = 1e-15;
    double theta[3] = { .5, .3, .2 };
    size_t x[1][2] = { {0, 1} };
    size_t y[1] = { 0 };

    double ll0_expected = log(theta[0]) - log(theta[0] + theta[1]);
    double ll0 = fullcond(1, 2, 3, 0, theta, x, y);
    ck_assert_double_eq_tol(ll0_expected, ll0, tol);

    double ll1_expected = -log(theta[0] + theta[1]);
    double ll1 = fullcond(1, 2, 3, 1, theta, x, y);
    ck_assert_double_eq_tol(ll1_expected, ll1, tol);
}
END_TEST

START_TEST(test_full_conditionals_last_player_has_played)
{
    double tol = 1e-15;
    double theta[3] = { .5, .3, .2 };
    size_t x[2][2] = { {0, 1}, {0, 2} };
    size_t y[2] = { 0, 0 };

    double ll0_expected = log(theta[0]) - log(theta[0] + theta[1])
                        + log(theta[0]) - log(theta[0] + 1 - theta[0] - theta[1]);
    double ll0 = fullcond(2, 2, 3, 0, theta, x, y);
    ck_assert_double_eq_tol(ll0_expected, ll0, tol);

    double ll1_expected = -log(theta[0] + theta[1])
                         - log(theta[0] + 1 - theta[0] - theta[1]);
    double ll1 = fullcond(2, 2, 3, 1, theta, x, y);
    ck_assert_double_eq_tol(ll1_expected, ll1, tol);
}
END_TEST

START_TEST(test_full_conditionals_last_player_has_played_2)
{
    double tol = 1e-15;
    double theta[3] = { .5, .3, .2 };
    size_t x[3][2] = { {0, 1}, {0, 1}, {0, 2} };
    size_t y[3] = { 0, 1, 0 };

    double ll0_expected = 2 * log(theta[0]) - 2 * log(theta[0] + theta[1]) - log(theta[0] + theta[2]);
    double ll0 = fullcond(3, 2, 3, 0, theta, x, y);
    ck_assert_double_eq_tol(ll0_expected, ll0, tol);

    double ll1_expected = log(theta[1]) - 2 * log(theta[0] + theta[1]) - log(theta[0] + theta[2]);
    double ll1 = fullcond(3, 2, 3, 1, theta, x, y);
    ck_assert_double_eq_tol(ll1_expected, ll1, tol);
}
END_TEST

START_TEST(test_full_conditionals_big_list)
{
    double tol = 1e-12;
    double theta[3] = { .5, .3, .2 };
    size_t x[10][2] = { {0, 1}, {0, 1},
                        {0, 2}, {0, 2}, {0, 2}, {0, 2}, {0, 2},
                        {1, 2}, {1, 2}, {1, 2} };
    size_t y[10] = { 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};

    double ll0_expected = 2 * log(theta[0]) + 4 * log(1-theta[0]-theta[1])
                         -2 * log(theta[0] + theta[1]) - 3 * log(1 - theta[0]);
    double ll0 = fullcond(10, 2, 3, 0, theta, x, y);
    ck_assert_double_eq_tol(ll0_expected, ll0, tol);

    double ll1_expected = 4 * log(theta[1]) + 4 * log(1 - theta[0] - theta[1])
                         -2 * log(theta[0] + theta[1]) - 5 * log(1-theta[1]);
    double ll1 = fullcond(10, 2, 3, 1, theta, x, y);
    ck_assert_double_eq_tol(ll1_expected, ll1, tol);
}
END_TEST

Suite * test_suite(void)
{
    Suite *s;
    TCase *tc_vector;
    TCase *tc_matrix;
    TCase *tc_core;
    TCase *tc_model;
    TCase *tc_simplex;

    s = suite_create("Tests");

    /* Vector test case */
    tc_vector = tcase_create("Vector");
    tcase_add_test(tc_vector, test_vector_new);
    tcase_add_test(tc_vector, test_vector_elements);
    suite_add_tcase(s, tc_vector);

    /* Matrix test case */
    tc_matrix = tcase_create("Matrix");
    tcase_add_test(tc_matrix, test_matrix_new);
    tcase_add_test(tc_matrix, test_matrix_elements);
    suite_add_tcase(s, tc_matrix);

    /* Core test case */
    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_sum);
    tcase_add_test(tc_core, test_log_sum_exp);
    tcase_add_test(tc_core, test_is_in);
    suite_add_tcase(s, tc_core);

    /* Model test case */
    tc_model = tcase_create("Model");
    tcase_add_test(tc_model, test_loglik);
    /* tcase_add_test(tc_model, test_full_conditionals_last_player_didnt_play); */
    /* tcase_add_test(tc_model, test_full_conditionals_last_player_has_played); */
    /* tcase_add_test(tc_model, test_full_conditionals_last_player_has_played_2); */
    tcase_add_test(tc_model, test_full_conditionals_big_list);
    suite_add_tcase(s, tc_model);

    tc_simplex = tcase_create("Simplex");
    tcase_add_test(tc_simplex, test_simplex_transform);
    suite_add_tcase(s, tc_simplex);

    return s;
}

int main(void) {
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = test_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);


    /* from the libcheck documentation:
     * (https://libcheck.github.io/check/doc/doxygen/html/check_8h.html)
     *
     * srunner_free()
     * Frees a suite runner, including all contained suite and test cases.
     *
     * This call is responsible for freeing all resources related to a suite runner
     * and all contained suites and test cases. Suite and test cases need not be
     * freed individually, as this call handles that.
     *
     */
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
