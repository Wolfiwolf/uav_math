#include "uav_math.h"
#include "stdlib.h"
#include <alloca.h>
#include <math.h>

void uav_matrix_init(struct Matrix *mat, uint8_t M, uint8_t N) {
	mat->M = M;
	mat->N = N;
	mat->rows = (float**)malloc(M * sizeof(uint32_t*));

	for (uint8_t i = 0; i < M; ++i) {
		mat->rows[i] = malloc(N * sizeof(float));
	}
}

void uav_matrix_destroy(struct Matrix *mat) {
	for (uint8_t i = 0; i < mat->M; ++i) {
		free(mat->rows[i]);
	}

	free(mat->rows);
}

void uav_matrix_add_to(struct Matrix *mat1, struct Matrix *mat2) {
	for (uint8_t i = 0; i < mat1->M; ++i) {
		for (uint8_t j = 0; j < mat1->N; ++j) 
			mat1->rows[i][j] = mat1->rows[i][j] + mat2->rows[i][j];
	}
} 

void uav_matrix_add(struct Matrix *mat1, struct Matrix *mat2, struct Matrix *res) {
	for (uint8_t i = 0; i < mat1->M; ++i) {
		for (uint8_t j = 0; j < mat1->N; ++j) 
			res->rows[i][j] = mat1->rows[i][j] + mat2->rows[i][j];
	}
}

void uav_matrix_sub_to(struct Matrix *mat1, struct Matrix *mat2) {
	for (uint8_t i = 0; i < mat1->M; ++i) {
		for (uint8_t j = 0; j < mat1->N; ++j) 
			mat1->rows[i][j] = mat1->rows[i][j] - mat2->rows[i][j];
	}
} 

void uav_matrix_sub(struct Matrix *mat1, struct Matrix *mat2, struct Matrix *res) {
	for (uint8_t i = 0; i < mat1->M; ++i) {
		for (uint8_t j = 0; j < mat1->N; ++j) 
			res->rows[i][j] = mat1->rows[i][j] - mat2->rows[i][j];
	}
} 

void uav_matrix_multiply(struct Matrix *mat1, struct Matrix *mat2, struct Matrix *res) {
	for (uint8_t i = 0; i < mat1->M; ++i) {
		for (uint8_t j = 0; j < mat1->N; ++j) {
			float sum = 0.0f;

			for (uint8_t k = 0; k < mat1->N; ++k) {
				sum += mat1->rows[i][k] * mat2->rows[k][j];
			}

			res->rows[i][j] = sum;
		}
	}
} 

void uav_matrix_scalar_multiply(struct Matrix *mat1, float scalar) {
	for (uint8_t i = 0; i < mat1->M; ++i) {
		for (uint8_t j = 0; j < mat1->N; ++j) 
			mat1->rows[i][j] *= scalar;
	}
} 

void uav_matrix_transpose(struct Matrix *mat, struct Matrix *res) {
	for (uint8_t i = 0; i < mat->M; ++i) {
		for (uint8_t j = 0; j < mat->N; ++j) 
			res->rows[i][j] = mat->rows[j][i];
	}
} 

float uav_matrix_max(struct Matrix *mat) {
	float max = mat->rows[0][0];

	for (uint8_t i = 0; i < mat->M; ++i) {
		for (uint8_t j = 0; j < mat->N; ++j)
			if (max < mat->rows[i][j]) max = mat->rows[i][j];
	}

	return max;
} 

float uav_matrix_min(struct Matrix *mat) {
	float min = mat->rows[0][0];

	for (uint8_t i = 0; i < mat->M; ++i) {
		for (uint8_t j = 0; j < mat->N; ++j)
			if (min < mat->rows[i][j]) min = mat->rows[i][j];
	}

	return min;
} 

// VECTOR MATH

float uav_vec_magnitude(struct Matrix *vec) {
	float sum = 0.0f;

	for (uint8_t i = 0; i < vec->M; ++i) {
		sum += vec->rows[i][0] * vec->rows[i][0];
	}

	return sqrtf(sum);
} 

float uav_vec_dot(struct Matrix *vec1, struct Matrix *vec2) {
	float res = 0.0f;

	for (uint8_t i = 0; i < vec1->M; ++i) {
		res += vec1->rows[i][0] * vec2->rows[i][0];
	}

	return res;
} 

void uav_vec_cross(struct Matrix *vec1, struct Matrix *vec2, struct Matrix *res) {
	res->rows[0][0] = (vec1->rows[1][0] * vec2->rows[2][0]) - (vec1->rows[2][0] * vec2->rows[1][0]);
	res->rows[1][0] = (vec1->rows[2][0] * vec2->rows[0][0]) - (vec1->rows[0][0] * vec2->rows[2][0]);
	res->rows[2][0] = (vec1->rows[0][0] * vec2->rows[1][0]) - (vec1->rows[1][0] * vec2->rows[0][0]);
} 

float uav_vec_angle_between(struct Matrix *vec1, struct Matrix *vec2) {
	return acosf(uav_vec_dot(vec1, vec2) / (uav_vec_magnitude(vec1) * uav_vec_magnitude(vec1)));
} 

void uav_vec_differential(struct Matrix *vec, struct Matrix *w, struct Matrix *res, float delta_t_sec) {
	uav_vec_cross(w, vec, res);
	uav_matrix_scalar_multiply(res, delta_t_sec);
} 

// ORIENTATION MATH

void uav_orient_euler_to_q(struct Matrix *euler_angles, struct Matrix *q) {
	float phi = euler_angles->rows[0][0];
	float theta = euler_angles->rows[1][0];
	float psi = euler_angles->rows[2][0];

	float C11 = cosf(theta) * cosf(psi);
	float C12 = cosf(theta) * sinf(psi);
	float C13 = -sinf(theta);

	float C21 = cosf(psi) * sinf(theta) * sinf(phi) - cosf(theta) * sin(psi);
	float C22 = cosf(phi) * cosf(psi) + sinf(theta) * sinf(phi) * sinf(psi);
	float C23 = cosf(theta) * sinf(phi);

	float C31 = cosf(phi) * cosf(psi) * sinf(theta) + sinf(phi) * sinf(psi);
	float C32 = -cosf(psi) * sinf(phi) + cosf(phi) * sinf(theta) * sinf(psi);
	float C33 = cosf(theta) * cosf(phi);

	struct Matrix q_h;
	uav_matrix_init(&q_h, 4, 1);

	q_h.rows[0][0] = sqrtf(0.25f * (1 + C11 + C22 + C33));
	q_h.rows[1][0] = sqrtf(0.25f * (1 + C11 - C22 - C33));
	q_h.rows[2][0] = sqrtf(0.25f * (1 - C11 + C22 - C33));
	q_h.rows[3][0] = sqrtf(0.25f * (1 - C11 - C22 + C33));

	float q_h_max = uav_matrix_max(&q_h);

	if (q_h_max == q_h.rows[0][0]) {
		q->rows[0][0] = q_h_max;
		q->rows[1][0] = (C23 - C32) / (4 * q_h_max);
		q->rows[2][0] = (C31 - C13) / (4 * q_h_max);
		q->rows[3][0] = (C12 - C21) / (4 * q_h_max);
	} else if (q_h_max == q_h.rows[1][0]) {
		q->rows[0][0] = (C23 - C32) / (4 * q_h_max);
		q->rows[1][0] = q_h_max;
		q->rows[2][0] = (C12 - C21) / (4 * q_h_max);
		q->rows[3][0] = (C31 - C13) / (4 * q_h_max);
	} else if (q_h_max == q_h.rows[2][0]) {
		q->rows[0][0] = (C31 - C13) / (4 * q_h_max);
		q->rows[1][0] = (C12 + C21) / (4 * q_h_max);
		q->rows[2][0] = q_h_max;
		q->rows[3][0] = (C23 + C32) / (4 * q_h_max);
	} else if (q_h_max == q_h.rows[3][0]) {
		q->rows[0][0] = (C12 - C21) / (4 * q_h_max);
		q->rows[1][0] = (C31 + C13) / (4 * q_h_max);
		q->rows[2][0] = (C23 + C32) / (4 * q_h_max);
		q->rows[3][0] = q_h_max;
	}

	uav_matrix_destroy(&q_h);
}

void uav_orient_q_to_euler(struct Matrix *q, struct Matrix *euler_angles) {
	float qs = q->rows[0][0];
	float qx = q->rows[1][0];
	float qy = q->rows[2][0];
	float qz = q->rows[3][0];

	float C23 = 2.0f * (qy * qz + qx * qs);
	float C33 = qs * qs - qx * qx + qy * qy - qz * qz;
	float C13 = 2.0f * (qx * qz - qy * qs);
	float C12 = 2.0f * (qx * qy + qz * qs);
	float C11 = qs * qs + qx * qx - qy * qy - qz * qz;

	euler_angles->rows[0][0] = atan2f(C23, C33);
	euler_angles->rows[1][0] = -asinf(C13);
	euler_angles->rows[2][0] = atan2f(C12, C11);
}

void uav_orient_q_dot(struct Matrix *q, struct Matrix *w, struct Matrix *res, float delta_t_sec) {
	struct Matrix omega;
	uav_matrix_init(&omega, 4, 4);

	float wx = w->rows[0][0];
	float wy = w->rows[1][0];
	float wz = w->rows[2][0];

	omega.rows[0][0] = 0.0f;
	omega.rows[0][1] = -wx;
	omega.rows[0][2] = -wy;
	omega.rows[0][3] = -wz;

	omega.rows[1][0] = wx;
	omega.rows[1][1] = 0.0f;
	omega.rows[1][2] = wz;
	omega.rows[1][3] = -wy;

	omega.rows[2][0] = wy;
	omega.rows[2][1] = -wz;
	omega.rows[2][2] = 0.0f;
	omega.rows[2][3] = wx;

	omega.rows[3][0] = wz;
	omega.rows[3][1] = wy;
	omega.rows[3][2] = -wx;
	omega.rows[3][3] = 0.0f;


	uav_matrix_multiply(&omega, q, res);
	uav_matrix_scalar_multiply(res, 0.5f * delta_t_sec);

	uav_matrix_destroy(&omega);
}
