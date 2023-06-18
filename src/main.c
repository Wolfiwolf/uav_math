#include <math.h>
#include <stdio.h>
#include "uav_math/uav_math.h"

int main() {
	float lat = 46.397551f;
	float lon = 15.633963f;
	float alt = 270.0f;

	float lat0 = 46.404747f;
	float lon0 = 15.636317f;
	float alt0 = 270.0f;

	float x, y, z;
	float x_r, y_r, z_r;
	float e, n, u;

	uav_trans_geodetic_to_ECEF(lat, lon, alt, &x, &y, &z);
	uav_trans_geodetic_to_ECEF(lat0, lon0, alt0, &x_r, &y_r, &z_r);
	uav_trans_ECEF_to_ENU(x, y, z, lat0, lon0, x_r, y_r, z_r, &e, &n, &u);

	printf("%.3f %.3f %.3f\n", lat, lon, alt);
	printf("%.3f %.3f %.3f\n", x, y, z);
	printf("%.3f %.3f %.3f\n", e, n, u);

	return 0;
}

