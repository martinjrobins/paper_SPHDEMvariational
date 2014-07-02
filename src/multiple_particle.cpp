/*
 * sphdem.cpp
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 7 Feb 2014
 *      Author: robinsonm
 */

#include "sphdem.h"
#include "Visualisation.h"

#include <vtkFloatArray.h>


int main(int argc, char **argv) {
	auto dem = DemType::New();
	auto sph = SphType::New();
	auto params = ptr<Params>(new Params());

	const int timesteps = 20000;
	const int nout = 100;
	const int timesteps_per_out = timesteps/nout;
	const double L = 0.004;
	const int nx = 10;


	 /* dem parameters
	 */
	const double porosity = 0.9;
	params->dem_diameter = 0.0001;
	params->dem_gamma = 0.0;
	params->dem_k = 10;
	params->dem_vol = (1.0/6.0)*PI*pow(params->dem_diameter,3);
	const double goal_dem_psep = pow(params->dem_vol/(1-porosity),1.0/NDIM);
	const double dem_dens = 2500.0;
	params->dem_mass = params->dem_vol*dem_dens;
	const double dem_min_reduced_mass = 0.5*params->dem_mass;
	params->dem_dt = (1.0/50.0)*PI/sqrt(params->dem_k/dem_min_reduced_mass-pow(0.5*params->dem_gamma/dem_min_reduced_mass,2));

	/*

	/*
	 * sph parameters
	 */
	params->sph_hfac = 1.5;
	params->sph_visc = 8.9e-07;
	params->sph_refd = 1000.0;
	params->sph_dens = 1000.0;
	params->sph_gamma = 7;
	const double VMAX = 2.0*sqrt(2*9.81*L);
	const double CSFAC = 10.0;
	params->sph_spsound = CSFAC*VMAX;
	params->sph_prb = pow(params->sph_refd/params->sph_dens,params->sph_gamma-1.0)*pow(params->sph_spsound,2)*params->sph_refd/params->sph_gamma;
	const double psep = 2*params->dem_diameter/params->sph_hfac;
	params->sph_dt = std::min(0.25*params->sph_hfac*psep/params->sph_spsound,0.125*pow(params->sph_hfac*psep,2)/params->sph_visc);
	params->sph_mass = params->sph_dens*pow(psep,NDIM);

	std::cout << "h = "<<params->sph_hfac*psep<<" vmax = "<<VMAX<<std::endl;

	params->sph_time_damping = params->sph_dt*500;
	params->dem_time_drop = params->sph_time_damping*2.0;

	params->time = 0;
	params->sph_maxh = params->sph_hfac*psep;


	/*
	 * define domain / geometry
	 */
	auto dem_geometry = [params](DemType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.8;
		REGISTER_DEM_PARTICLE(i);
		const double dem_diameter = params->dem_diameter;
		const double dem_k = params->dem_k;
		const double dem_gamma = params->dem_gamma;
		const double dem_mass = params->dem_mass;
		const double dem_vol = params->dem_vol;

		const double overlap = dem_diameter/2.0-r[2];
		if (overlap>0) {
			const double overlap_dot = -v[2];
			const Vect3d normal(0,0,1);
			acceleration += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
		}
		
		return acceleration;
	};

	auto sph_geometry = [](SphType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.8;
		return acceleration;
	};

	const double dem_width = L/4;
	const Vect3d min(0,0,0);
	const Vect3d max(L,L,L);
	const Vect3d dem_min(0,0,0.5*(L-dem_width));
	const Vect3d dem_max(L,L,0.5*(L+dem_width));
	const Vect3b periodic(true,true,false);
	const Vect3i sph_n = ((max-min)/psep).cast<int>();
	const Vect3d sph_psep = (max-min).cwiseQuotient(sph_n.cast<double>());
	const Vect3d min_domain(0,0,-3*sph_psep[2]);
	const Vect3d max_domain(max[0],max[1],2.0*L);

	const Vect3i dem_n = ((dem_max-dem_min)/goal_dem_psep).cast<int>();
	const Vect3d dem_psep = (dem_max-dem_min).cwiseQuotient(dem_n.cast<double>());


	/*
	 * create sph and dem particles
	 */
	dem->create_particles_grid(dem_min,dem_max,dem_n,[](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);
		v << 0,0,0;
		v0 << 0,0,0;
		f << 0,0,0;
		f0 << 0,0,0;
	});

	sph->create_particles_grid(min_domain,max,sph_n+Vect3i(0,0,3),[psep,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		h = params->sph_maxh;
		omega = 1.0;
		kappa = 0.0;
		v << 0,0,0;
		v0 << 0,0,0;
		dddt = 0;
		f0 << 0,0,0;
		e = 1;
		rho = params->sph_dens;
		f << 0,0,0;
		fdrag << 0,0,0;
		f0 << 0,0,0;
		if (r[2]<0) {
			fixed = true;
		} else {
			fixed = false;
		}
	});

	/*
	 * setup gamma calculation
	 */
	ptr<GammaEval> gamma = ptr<GammaEval>(new GammaEval(params->dem_diameter));
	gamma->reset_limits(params->sph_maxh*0.5,params->sph_maxh*1.5,0,2*params->sph_maxh+params->dem_diameter/2.0,100,100);


	/*
	 * setup output stuff
	 */
	auto sph_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto dem_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	sph->copy_to_vtk_grid(sph_grid);
	dem->copy_to_vtk_grid(dem_grid);
	Visualisation::vtkWriteGrid("vis/at_start_sph",0,sph_grid);
	Visualisation::vtkWriteGrid("vis/at_start_dem",0,dem_grid);


	std::cout << "init porosity and rho...."<<std::endl;
	sph->init_neighbour_search(min_domain,max_domain,2*params->sph_maxh,periodic);
	dem->init_neighbour_search(min_domain,max_domain,2*params->sph_maxh,periodic);

	/*
	 * init porosity and rho
	 */
	//std::cout << "calculate omega"<<std::endl;
	std::for_each(sph->begin(),sph->end(),[sph,dem,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		e = 1;
		//bool found = false;
		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			const double r = sqrt(r2);
			const double q = r/h;
			const double Wab = W(q,h);
			e -= params->dem_vol*Wab;
			//found = true;
		}
		rho = params->sph_dens/e;
	});


	/*
	 * Simulate!!!!
	 */
	std::cout << "starting...."<<std::endl;
	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			sphdem(sph,dem,params,sph_geometry,dem_geometry,gamma);
		}
		std::cout <<"iteration "<<i<<std::endl;
		
		sph->copy_to_vtk_grid(sph_grid);
		dem->copy_to_vtk_grid(dem_grid);
		Visualisation::vtkWriteGrid("vis/sph",i,sph_grid);
		Visualisation::vtkWriteGrid("vis/dem",i,dem_grid);
	}
	
	
}
