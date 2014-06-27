/*
 * GammaEval.cpp
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of SPHDEMvariational.
 *
 * SPHDEMvariational is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SPHDEMvariational is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with SPHDEMvariational.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 27 Jun 2014
 *      Author: robinsonm
 */

#include "GammaEval.h"
#include "sph_common.h"
#include "gauss_2d_sphere.h"


double f(double x, double y, void* data) {
	const double r = sqrt(x*x+y*y);
	const double h = (double)data[0];
	const double q = r/h;
	return W(r/h,r);
}



void GammaEval::reset_limits(const double _h_min, const double _h_max,
		const double _r_min, const double _r_max,
		const double _n_h,const double _n_r) {

	const double dh = (h_max-h_min)/n_h;
	const double dr = (r_max-r_min)/n_r;

	h_min = _h_min;
	h_max = _h_max;
	r_min = _r_min;
	r_max = _r_max;
	n_h = _n_h;
	n_r = _n_r;

	gamma.resize(n_h+3,n_r+3);
	gamma_h.resize(n_h+3,n_r+3);
	gamma_r.resize(n_h+3,n_r+3);

	for (int i = 0; i <= n_h+2; ++i) {
		const double h = h_min + dh*(i-1);
		for (int j = 0; j <= n_r+2; ++j) {
			const double r = r_min + dr*(j-1);
			gamma(i,j) = gauss_product_2D_sphere(25,f,(void *)&h,d/2.0,r,0.0);
		}
	}

	for (int i = 0; i <= n_h+2; ++i) {
		const double h = h_min + dh*(i-1);
		for (int j = 0; j <= n_r+2; ++j) {
			const double r = r_min + dr*(j-1);
			if (i < 2) {
				gamma_h(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i+1,j) - (3.0)*gamma(i+2,j) + (4.0/3.0)*gamma(i+3,j) - (1.0/4.0)*gamma(i+4,j);
			} else if (i < 1) {
				gamma_h(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i+1,j) - (3.0)*gamma(i+2,j) + (4.0/3.0)*gamma(i+3,j) - (1.0/4.0)*gamma(i+4,j);
			} else if (i > n_h) {
				gamma_h(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i-1,j) + (3.0)*gamma(i-2,j) - (4.0/3.0)*gamma(i-3,j) + (1.0/4.0)*gamma(i-4,j);
			} else if (i > n_h+1) {
				gamma_h(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i-1,j) + (3.0)*gamma(i-2,j) - (4.0/3.0)*gamma(i-3,j) + (1.0/4.0)*gamma(i-4,j);
			} else {
				gamma_h(i,j) = (1.0/12.0)*gamma(i-2,j) - (2.0/3.0)*gamma(i-1,j) + (2.0/3.0)*gamma(i+1,j) - (1.0/12.0)*gamma(i+2,j);
			}

			if (j < 2) {
				gamma_r(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i,j+1) - (3.0)*gamma(i,j+2) + (4.0/3.0)*gamma(i,j+3) - (1.0/4.0)*gamma(i,j+4);
			} else if (j < 1) {
				gamma_r(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i,j+1) - (3.0)*gamma(i,j+2) + (4.0/3.0)*gamma(i,j+3) - (1.0/4.0)*gamma(i,j+4);
			} else if (j > n_r) {
				gamma_r(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i,j-1) + (3.0)*gamma(i,j-2) - (4.0/3.0)*gamma(i,j-3) + (1.0/4.0)*gamma(i,j-4);
			} else if (j > n_r+1) {
				gamma_r(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i,j-1) + (3.0)*gamma(i,j-2) - (4.0/3.0)*gamma(i,j-3) + (1.0/4.0)*gamma(i,j-4);
			} else {
				gamma_r(i,j) = (1.0/12.0)*gamma(i,j-2) - (2.0/3.0)*gamma(i,j-1) + (2.0/3.0)*gamma(i,j+1) - (1.0/12.0)*gamma(i,j+2);
			}
		}
	}

}

double bilinearInterpolate (double p[2][2], const double x, const double y) {
	return (1-x)*(1-y)*p[0][0] + x*(1-y)*p[1][0] + (1-x)*y*p[0][1] + x*y*p[1][1];
}

double cubicInterpolate (double p[4], const double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate (double p[4][4], const double x, const double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}

#define LINEAR
double GammaEval::interpolate(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &data, const double h, const double r) {
	const double dh = (h_max-h_min)/n_h;
	const double dr = (r_max-r_min)/n_r;
	const int i0 = (int)floor((h-h_min)/dh);
	const int j0 = (int)floor((r-r_min)/dr);
	const double x = (h-h_min)/dh-i0;
	const double y = (r-r_min)/dr-j0;

#ifndef LINEAR
	double p[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			p[i][j] = data(i0+i-1,j0+j-1);
		}
	}
	return  bicubicInterpolate(p,x,y);
#else
	double p[2][2];
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			p[i][j] = data(i0+i,j0+j);
		}
	}
	return  bilinearInterpolate(p,x,y);
#endif
}

double GammaEval::get_gamma(const double h, const double r) {
	return interpolate(gamma,h,r);
}

double GammaEval::get_gamma_h(const double h, const double r) {
	return interpolate(gamma_h,h,r);
}

double GammaEval::get_gamma_r(const double h, const double r) {
	return interpolate(gamma_r,h,r);
}


