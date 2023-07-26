#include "uav_math.h"
#include "stdlib.h"
#include <alloca.h>
#include <math.h>


#define RAD_TO_DEG 57.295779513082320876798154814105

void uav_matrix_init(struct Matrix *mat, uint8_t M, uint8_t N) {
	mat->M = M;
	mat->N = N;
	mat->rows = (float**)malloc(M * sizeof(uint32_t*));

	for (uint8_t i = 0; i < M; ++i) {
		mat->rows[i] = (float*)malloc(N * sizeof(float));
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

void uav_matrix_copy(struct Matrix *src, struct Matrix *dest) {
	for (uint8_t i = 0; i < src->M; ++i) {
		for (uint8_t j = 0; j < src->N; ++j) 
			dest->rows[i][j] = src->rows[i][j];
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

// ROTATIONS MATH

void uav_rotation_body_to_inertial(struct Matrix *vec, struct Matrix *euler_angles) {
	struct Matrix Rx;
	uav_matrix_init(&Rx, 3, 3);
	struct Matrix Ry;
	uav_matrix_init(&Ry, 3, 3);
	struct Matrix Rz;
	uav_matrix_init(&Rz, 3, 3);
	struct Matrix Rxy;
	uav_matrix_init(&Rxy, 3, 3);
	struct Matrix R;
	uav_matrix_init(&R, 3, 3);
	struct Matrix Rt;
	uav_matrix_init(&Rt, 3, 3);
	struct Matrix temp_vec;
	uav_matrix_init(&temp_vec, 3, 1);

	float phi = euler_angles->rows[0][0];
	float theta = euler_angles->rows[1][0];
	float psi = euler_angles->rows[2][0];

	// Rx

	Rx.rows[0][0] = 1.0f;
	Rx.rows[0][1] = 0.0f;
	Rx.rows[0][2] = 0.0f;

	Rx.rows[1][0] = 0.0f;
	Rx.rows[1][1] = cosf(phi);
	Rx.rows[1][2] = -sinf(phi);

	Rx.rows[2][0] = 0.0f;
	Rx.rows[2][1] = sinf(phi);
	Rx.rows[2][2] = cosf(phi);

	// Ry

	Ry.rows[0][0] = cosf(theta);
	Ry.rows[0][1] = 0.0f;
	Ry.rows[0][2] = sinf(theta);

	Ry.rows[1][0] = 0.0f;
	Ry.rows[1][1] = 1.0f;
	Ry.rows[1][2] = 0.0f;

	Ry.rows[2][0] = -sinf(theta);
	Ry.rows[2][1] = 0.0f;
	Ry.rows[2][2] = cosf(theta);

	// Rz

	Rz.rows[0][0] = cosf(psi);
	Rz.rows[0][1] = -sinf(psi);
	Rz.rows[0][2] = 0.0f;

	Rz.rows[1][0] = sinf(psi);
	Rz.rows[1][1] = cosf(psi);
	Rz.rows[1][2] = 0.0f;

	Rz.rows[2][0] = 0.0f;
	Rz.rows[2][1] = 0.0f;
	Rz.rows[2][2] = 1.0f;

	uav_matrix_multiply(&Ry, &Rx, &Rxy);


	uav_matrix_multiply(&Rz, &Rxy, &R);

	uav_matrix_transpose(&R, &Rt);


	uav_matrix_multiply(&Rt, vec, &temp_vec);
	uav_matrix_copy(&temp_vec, vec);

	uav_matrix_destroy(&Rx);
	uav_matrix_destroy(&Ry);
	uav_matrix_destroy(&Rz);
	uav_matrix_destroy(&Rxy);
	uav_matrix_destroy(&R);
	uav_matrix_destroy(&Rt);
	uav_matrix_destroy(&temp_vec);

}

void uav_rotation_inertial_to_body(struct Matrix *vec, struct Matrix *euler_angles) {
	struct Matrix Rx;
	uav_matrix_init(&Rx, 3, 3);
	struct Matrix Ry;
	uav_matrix_init(&Ry, 3, 3);
	struct Matrix Rz;
	uav_matrix_init(&Rz, 3, 3);
	struct Matrix Rxy;
	uav_matrix_init(&Rxy, 3, 3);
	struct Matrix R;
	uav_matrix_init(&R, 3, 3);
	struct Matrix temp_vec;
	uav_matrix_init(&temp_vec, 3, 1);

	float phi = euler_angles->rows[0][0];
	float theta = euler_angles->rows[1][0];
	float psi = euler_angles->rows[2][0];

	// Rx

	Rx.rows[0][0] = 1.0f;
	Rx.rows[0][1] = 0.0f;
	Rx.rows[0][2] = 0.0f;

	Rx.rows[1][0] = 0.0f;
	Rx.rows[1][1] = cosf(phi);
	Rx.rows[1][2] = -sinf(phi);

	Rx.rows[2][0] = 0.0f;
	Rx.rows[2][1] = sinf(phi);
	Rx.rows[2][2] = cosf(phi);

	// Ry

	Ry.rows[0][0] = cosf(theta);
	Ry.rows[0][1] = 0.0f;
	Ry.rows[0][2] = sinf(theta);

	Ry.rows[1][0] = 0.0f;
	Ry.rows[1][1] = 1.0f;
	Ry.rows[1][2] = 0.0f;

	Ry.rows[2][0] = -sinf(theta);
	Ry.rows[2][1] = 0.0f;
	Ry.rows[2][2] = cosf(theta);

	// Rz

	Rz.rows[0][0] = cosf(psi);
	Rz.rows[0][1] = -sinf(psi);
	Rz.rows[0][2] = 0.0f;

	Rz.rows[1][0] = sinf(psi);
	Rz.rows[1][1] = cosf(psi);
	Rz.rows[1][2] = 0.0f;

	Rz.rows[2][0] = 0.0f;
	Rz.rows[2][1] = 0.0f;
	Rz.rows[2][2] = 1.0f;


	uav_matrix_multiply(&Ry, &Rx, &Rxy);

	uav_matrix_multiply(&Rz, &Rxy, &R);

	uav_matrix_multiply(&R, vec, &temp_vec);
	uav_matrix_copy(&temp_vec, vec);

	uav_matrix_destroy(&Rx);
	uav_matrix_destroy(&Ry);
	uav_matrix_destroy(&Rz);
	uav_matrix_destroy(&Rxy);
	uav_matrix_destroy(&R);
	uav_matrix_destroy(&temp_vec);
}

// ORIENTATION MATH

void uav_orient_euler_to_q(struct Matrix *euler_angles, struct Matrix *q) {
	float phi = euler_angles->rows[0][0];
	float theta = euler_angles->rows[1][0];
	float psi = euler_angles->rows[2][0];

	float C11 = cosf(theta) * cosf(psi);
	float C12 = cosf(theta) * sinf(psi);
	float C13 = -sinf(theta);

	float C21 = sinf(phi) * sinf(theta) * cosf(psi) - cosf(phi) * sinf(psi);
	float C22 = cosf(phi) * cosf(psi) + sinf(phi) * sinf(theta) * sinf(psi);
	float C23 = sinf(phi) * cosf(theta);

	float C31 = cosf(phi) * sinf(theta) * cosf(psi) + sinf(phi) * sinf(psi);
	float C32 = cosf(phi) * sinf(theta) * sinf(psi) - sinf(phi) * cosf(psi);
	float C33 = cosf(phi) * cosf(theta);

	struct Matrix q_h;
	uav_matrix_init(&q_h, 4, 1);

	q_h.rows[0][0] = sqrtf(fabsf(0.25f * (1 + C11 + C22 + C33)));
	q_h.rows[1][0] = sqrtf(fabsf(0.25f * (1 + C11 - C22 - C33)));
	q_h.rows[2][0] = sqrtf(fabsf(0.25f * (1 - C11 + C22 - C33)));
	q_h.rows[3][0] = sqrtf(fabsf(0.25f * (1 - C11 - C22 + C33)));

	float q_h_max = uav_matrix_max(&q_h);

	if (q_h_max == q_h.rows[0][0]) {
		q->rows[0][0] = q_h_max;
		q->rows[1][0] = (C23 - C32) / (4.0f * q_h_max); 
		q->rows[2][0] = (C31 - C13) / (4.0f * q_h_max);
		q->rows[3][0] = (C12 - C21) / (4.0f * q_h_max);
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
	float C33 = qs * qs - qx * qx - qy * qy + qz * qz; 
	float C13 = 2.0f * (qx * qz - qy * qs); 
	float C12 = 2.0f * (qx * qy + qz * qs); 
	float C11 = qs * qs + qx * qx - qy * qy - qz * qz;

	if (C13 >= 1.0f) C13 = 1.0f;
	else if (C13 <= -1.0f) C13 = -1.0f;

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

// COORDINATE SYSTEM TRANSFORMATIONS

#define DEGREES_TO_RADIANS 0.0174532925
#define WGS84_A 6378137.0f
#define WGS84_E 0.0818191908426f

void uav_trans_geodetic_to_ECEF(float lat, float lon, float alt, float *x, float *y, float *z) {
    double clat = cos(lat * DEGREES_TO_RADIANS);
    double slat = sin(lat * DEGREES_TO_RADIANS);
    double clon = cos(lon * DEGREES_TO_RADIANS);
    double slon = sin(lon * DEGREES_TO_RADIANS);

    double N = WGS84_A / sqrt(1.0 - WGS84_E * WGS84_E * slat * slat);

    *x = (N + alt) * clat * clon;
    *y = (N + alt) * clat * slon;
    *z = (N * (1.0 - WGS84_E * WGS84_E) + alt) * slat;
}

void uav_trans_ECEF_to_geodetic(float x, float y, float z, float *lat, float *lon, float *alt) {
    // WGS84 constants
    double a = 6378137.0;
    double f = 1.0 / 298.257223563;

    // Derived constants
    double b = a - f * a;
    double e = sqrt(pow(a, 2.0) - pow(b, 2.0)) / a;
    double clambda = atan2(y, x);
    double p = sqrt(pow(x, 2.0) + pow(y, 2.0));
    double h_old = 0.0;

    // First guess with h=0 meters
    double theta = atan2(z, p * (1.0 - pow(e, 2.0)));
    double cs = cos(theta);
    double sn = sin(theta);
    double N = pow(a, 2.0) / sqrt(pow(a * cs, 2.0) + pow(b * sn, 2.0));
    double h = p / cs - N;

    while (fabs(h - h_old) > 1.0e-6) {
        h_old = h;
        theta = atan2(z, p * (1.0 - pow(e, 2.0) * N / (N + h)));
        cs = cos(theta);
        sn = sin(theta);
        N = pow(a, 2.0) / sqrt(pow(a * cs, 2.0) + pow(b * sn, 2.0));
        h = p / cs - N;
    }

    // Convert radians to degrees
    *lon = clambda * RAD_TO_DEG;
    *lat = theta * RAD_TO_DEG;
    *alt = h;
}

void uav_trans_ECEF_to_ENU(float x, float y, float z, float lat_r, float lon_r, float x_r, float y_r, float z_r, float *e, float *n, float *u) {
    double clat = cos(lat_r * DEGREES_TO_RADIANS);
    double slat = sin(lat_r * DEGREES_TO_RADIANS);
    double clon = cos(lon_r * DEGREES_TO_RADIANS);
    double slon = sin(lon_r * DEGREES_TO_RADIANS);
    double dx = x - x_r;
    double dy = y - y_r;
    double dz = z - z_r;

    *e = -slon*dx  + clon*dy;
    *n = -slat*clon*dx - slat*slon*dy + clat*dz;
    *u = clat*clon*dx + clat*slon*dy + slat*dz;
}

void uav_trans_ENU_to_ECEF(float e, float n, float u, float lat_r, float lon_r, float x_r, float y_r, float z_r, float *x, float *y, float *z) {
	 double a = 6378137.0;        
    double f = 1 / 298.257223563;
    double b = a * (1 - f);     

    lat_r *= 1 / RAD_TO_DEG;
    lon_r *= 1 / RAD_TO_DEG;

    double sinRefLat = sin(lat_r);
    double cosRefLat = cos(lat_r);
    double sinRefLon = sin(lon_r);
    double cosRefLon = cos(lon_r);

    double dX = -sinRefLon * e - sinRefLat * cosRefLon * n + cosRefLat * cosRefLon * u;
    double dY = cosRefLon * e - sinRefLat * sinRefLon * n + cosRefLat * sinRefLon * u;
    double dZ = cosRefLat * n + sinRefLat * u;

    *x = x_r + dX;
    *y = y_r + dY;
    *z = z_r + dZ;
}
