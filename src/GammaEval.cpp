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



void GammaEval::reset_limits(const double h_min, const double h_max,
		const double r_min, const double r_max,
		const double n_h,const double n_r) {

	gamma.resize(n_h+1,n_r+1);
	gamma_h.resize(n_h+1,n_r+1);
	gamma_r.resize(n_h+1,n_r+1);

	const double dh = (h_max-h_min)/n_h;
	const double dr = (r_max-r_min)/n_r;

	for (int i = 0; i <= n_h; ++i) {
		const double h = h_min + dh*i;
		for (int j = 0; j <= n_r; ++j) {
			const double r = r_min + dr*i;
			gamma(i,j) = gauss_product_2D_sphere(25,f,(void *)&h,d/2.0,r,0.0);
		}
	}

	for (int i = 0; i <= n_h; ++i) {
		const double h = h_min + dh*i;
		for (int j = 0; j <= n_r; ++j) {
			const double r = r_min + dr*i;
			if (i < 2) {
				gamma_h(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i+1,j) - (3.0)*gamma(i+2,j) + (4.0/3.0)*gamma(i+3,j) - (1.0/4.0)*gamma(i+4,j);
			} else if (i < 1) {
				gamma_h(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i+1,j) - (3.0)*gamma(i+2,j) + (4.0/3.0)*gamma(i+3,j) - (1.0/4.0)*gamma(i+4,j);
			} else if (i > n_h-2) {
				gamma_h(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i-1,j) + (3.0)*gamma(i-2,j) - (4.0/3.0)*gamma(i-3,j) + (1.0/4.0)*gamma(i-4,j);
			} else if (i > n_h-1) {
				gamma_h(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i-1,j) + (3.0)*gamma(i-2,j) - (4.0/3.0)*gamma(i-3,j) + (1.0/4.0)*gamma(i-4,j);
			} else {
				gamma_h(i,j) = (1.0/12.0)*gamma(i-2,j) - (2.0/3.0)*gamma(i-1,j) + (2.0/3.0)*gamma(i+1,j) - (1.0/12.0)*gamma(i+2,j);
			}

			if (j < 2) {
				gamma_r(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i,j+1) - (3.0)*gamma(i,j+2) + (4.0/3.0)*gamma(i,j+3) - (1.0/4.0)*gamma(i,j+4);
			} else if (j < 1) {
				gamma_r(i,j) = -(25.0/12.0)*gamma(i,j) + (4.0)*gamma(i,j+1) - (3.0)*gamma(i,j+2) + (4.0/3.0)*gamma(i,j+3) - (1.0/4.0)*gamma(i,j+4);
			} else if (j > n_h-2) {
				gamma_r(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i,j-1) + (3.0)*gamma(i,j-2) - (4.0/3.0)*gamma(i,j-3) + (1.0/4.0)*gamma(i,j-4);
			} else if (j > n_h-1) {
				gamma_r(i,j) = +(25.0/12.0)*gamma(i,j) - (4.0)*gamma(i,j-1) + (3.0)*gamma(i,j-2) - (4.0/3.0)*gamma(i,j-3) + (1.0/4.0)*gamma(i,j-4);
			} else {
				gamma_r(i,j) = (1.0/12.0)*gamma(i,j-2) - (2.0/3.0)*gamma(i,j-1) + (2.0/3.0)*gamma(i,j+1) - (1.0/12.0)*gamma(i,j+2);
			}
		}
	}

}

double GammaEval::get_gamma(const double h, const double r) {
	const double dh = (h_max-h_min)/n_h;
	const double dr = (r_max-r_min)/n_r;
	const int i0 = (int)floor((h-h_min)/dh);
	const double h0 = i0*dh;
	const double h1 = (i0+1)*dh;

	const int j0 = (int)floor((r-r_min)/dr);
	const double r0 = j0*dr;
	const double r1 = (j0+1)*dr;


	return  gamma(i0,j0) + gamma(i0+1,j0) + gamma(i0,j0) + gamma(i0,j0+1) + gamma(i0+1,j0+1);
}

double GammaEval::get_gamma_h(const double h, const double r) {
}

double GammaEval::get_gamma_r(const double h, const double r) {
}


