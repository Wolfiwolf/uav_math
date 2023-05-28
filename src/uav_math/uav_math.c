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
	
}

void uav_orient_q_to_euler(struct Matrix *q, struct Matrix *euler_angles) {

}

void uav_orient_q_dot(struct Matrix *q, struct Matrix *w, struct Matrix *res, float delta_t_sec) {

}
