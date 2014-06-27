/*
 * GammaEval.h
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

#ifndef GAMMAEVAL_H_
#define GAMMAEVAL_H_

#include "Eigen/Dense"

class GammaEval {
public:
	GammaEval(const double d):d(d),h_min(0),h_max(-1),r_min(0),r_max(-1) {

	}
	virtual ~GammaEval();
	void reset_limits(const double h_min,const double h_max,const double r_min,const double r_max,const double n_h,const double n_r);
	double get_gamma(const double h, const double r);
	double get_gamma_h(const double h, const double r);
	double get_gamma_r(const double h, const double r);
private:
	const double d;
	const double h_min,h_max,r_min,r_max;
	Eigen::Matrix<double, Dynamic, Dynamic> gamma,gamma_h,gamma_r;

};

#endif /* GAMMAEVAL_H_ */
