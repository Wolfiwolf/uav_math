#include <math.h>
#include <stdio.h>
#include "uav_math/uav_math.h"

int main() {
    struct Matrix vec1;
    uav_matrix_init(&vec1, 3, 1);

    vec1.rows[0][0] = 0.0f;
    vec1.rows[1][0] = 0.0f;
    vec1.rows[2][0] = 1.0f;

    struct Matrix euler_angles;
    uav_matrix_init(&euler_angles, 3, 1);

    euler_angles.rows[0][0] = 0.0f * 0.01745329f;
    euler_angles.rows[1][0] = 45.0f * 0.01745329f;
    euler_angles.rows[2][0] = 0.0f;

    uav_rotation_body_to_inertial(&vec1, &euler_angles);

    printf("x: %.3f\n", vec1.rows[0][0]);
    printf("y: %.3f\n", vec1.rows[1][0]);
    printf("z: %.3f\n", vec1.rows[2][0]);

    uav_matrix_destroy(&vec1);
    uav_matrix_destroy(&euler_angles);

	return 0;
}

