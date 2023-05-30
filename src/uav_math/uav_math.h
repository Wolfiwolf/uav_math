#ifndef UAV_MATH_H
#define UAV_MATH_H

#include <stdint.h>

struct Matrix {
	uint8_t M;
	uint8_t N;
	float **rows;
};

void uav_matrix_init(struct Matrix *mat, uint8_t M, uint8_t N); 

void uav_matrix_destroy(struct Matrix *mat);

// MATRIX MATH

void uav_matrix_add_to(struct Matrix *mat1, struct Matrix *mat2); 

void uav_matrix_add(struct Matrix *mat1, struct Matrix *mat2, struct Matrix *res); 

void uav_matrix_sub_to(struct Matrix *mat1, struct Matrix *mat2); 

void uav_matrix_sub(struct Matrix *mat1, struct Matrix *mat2, struct Matrix *res); 

void uav_matrix_multiply(struct Matrix *mat1, struct Matrix *mat2, struct Matrix *res); 

void uav_matrix_scalar_multiply(struct Matrix *mat1, float scalar); 

void uav_matrix_transpose(struct Matrix *mat, struct Matrix *res); 

void uav_matrix_transpose(struct Matrix *mat, struct Matrix *res); 

float uav_matrix_max(struct Matrix *mat); 

float uav_matrix_min(struct Matrix *mat); 

// VECTOR MATH

float uav_vec_magnitude(struct Matrix *vec); 

float uav_vec_dot(struct Matrix *vec1, struct Matrix *vec2); 

void uav_vec_cross(struct Matrix *vec1, struct Matrix *vec2, struct Matrix *res); 

float uav_vec_angle_between(struct Matrix *vec1, struct Matrix *vec2); 

void uav_vec_differential(struct Matrix *vec, struct Matrix *w, struct Matrix *res, float delta_t_sec); 

// ROTATIONS MATH

void uav_rotation_inertial_to_body(struct Matrix *vec, struct Matrix *euler_angles);

void uav_rotation_body_to_inertial(struct Matrix *vec, struct Matrix *euler_angles);

// ORIENTATION MATH

void uav_orient_euler_to_q(struct Matrix *euler_angles, struct Matrix *q);

void uav_orient_q_to_euler(struct Matrix *q, struct Matrix *euler_angles);

void uav_orient_q_dot(struct Matrix *q, struct Matrix *w, struct Matrix *res, float delta_t_sec);

#endif

