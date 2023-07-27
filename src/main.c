#include <math.h>
#include <stdio.h>
#include "uav_math/uav_math.h"

int main() {
    struct Matrix vec1;
    uav_matrix_init(&vec1, 3, 1);
    struct Matrix vec2;
    uav_matrix_init(&vec2, 3, 1);

    vec1.rows[0][0] = 1.0f;
    vec1.rows[1][0] = 1.0f;
    vec1.rows[2][0] = 0.0f;

    vec2.rows[0][0] = 0.0f;
    vec2.rows[1][0] = 1.0f;
    vec2.rows[2][0] = 1.0f;

    float angle = uav_vec_angle_between(&vec1, &vec2);

	printf("%.3f\n", angle);

	return 0;
}

