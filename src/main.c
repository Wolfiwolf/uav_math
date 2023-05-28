#include <stdio.h>
#include "uav_math/uav_math.h"

int main() {

	struct Matrix mat;
	uav_matrix_init(&mat, 3, 3);

	struct Matrix I;
	uav_matrix_init(&I, 3, 3);

	struct Matrix res;
	uav_matrix_init(&res, 3, 3);


	for (uint8_t i = 0; i < 3; ++i) {
		for (uint8_t j = 0; j < 3; ++j) {
			mat.rows[i][j] = i * 3 + j;
		}
	}

	for (uint8_t i = 0; i < 3; ++i) {
		for (uint8_t j = 0; j < 3; ++j) {
			I.rows[i][j] = i * 3 + j;
		}
	}


	uav_matrix_transpose(&mat, &res);


	for (uint8_t i = 0; i < 3; ++i) {
		for (uint8_t j = 0; j < 3; ++j) {
			printf("mat %.3f\n", res.rows[i][j]);
		}
	}

	return 0;
}

