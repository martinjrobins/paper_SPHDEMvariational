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
#include <string>
#include <iostream>     // std::cout, std::right, std::endl
#include <stdio.h>


#include <vtkFloatArray.h>


int main(int argc, char **argv) {
	auto dem = DemType::New();
	auto sph = SphType::New();
	auto params = ptr<Params>(new Params());

	const int timesteps = 2000;
	const int nout = 100;
	const int timesteps_per_out = timesteps/nout;
	const double L = 0.004;
	const int nx = 20;

	const double Re = 1.0;


	if (argc < 2) {
		params->dem_diameter = 0.0006;
	} else {
		params->dem_diameter = 0.0006*atof(argv[1]);
	}
	std::cout << "using diameter = "<<params->dem_diameter<<std::endl;
	char dir[100];
	sprintf(dir, "results_diameter_%s", argv[1]);

	 /* dem parameters
	 */
	params->dem_gamma = 0.0;
	params->dem_k = 10;
	params->dem_vol = (1.0/6.0)*PI*pow(params->dem_diameter,3);
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
	params->sph_inletv = Re*params->sph_visc/params->dem_diameter;
	params->sph_refd = 1000.0;
	params->sph_dens = 1000.0;
	params->sph_gamma = 7;
	//const double VMAX = 2.0*sqrt(2*9.81*L);
	const double VMAX = params->sph_inletv;

	const double CSFAC = 10.0;
	params->sph_spsound = CSFAC*VMAX;
	params->sph_prb = pow(params->sph_refd/params->sph_dens,params->sph_gamma-1.0)*pow(params->sph_spsound,2)*params->sph_refd/params->sph_gamma;
	const double psep = L/nx;
	params->sph_dt = std::min(0.25*params->sph_hfac*psep/params->sph_spsound,0.125*pow(params->sph_hfac*psep,2)/params->sph_visc);
	params->sph_mass = params->sph_dens*pow(psep,NDIM);

	std::cout << "h = "<<params->sph_hfac*psep<<" vmax = "<<VMAX<<" psep = "<<psep<<std::endl;
	std::cout << "sph_dt = "<<params->sph_dt<<" 1/20 m/b = "<<(1.0/20.0)*params->dem_mass/(3.0*PI*params->sph_dens*params->sph_visc*params->dem_diameter)<<std::endl;
	std::cout << "dem_dt = "<<params->dem_dt<<std::endl;
	std::cout << "2*params->sph_maxh + params->dem_diameter/2.0 = "<<2*params->sph_maxh + params->dem_diameter/2.0<<std::endl;
	std::cout << " Re = "<<params->dem_diameter*params->sph_inletv/params->sph_visc<<std::endl;

	params->dem_time_drop = params->sph_dt*timesteps;
	params->time = 0;
	params->sph_maxh = params->sph_hfac*psep;


	/*
	 * define domain / geometry
	 */
	auto dem_geometry = [params](DemType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,0;
		return acceleration;
	};

	auto sph_geometry = [](SphType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,0;
		return acceleration;
	};

	const Vect3d min(0,0,-3.0*psep);
	const Vect3d max(L,L,L);
	const Vect3b periodic(true,true,false);

	dem->create_particles(1,[L](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);
		v << 0,0,0;
		v0 << 0,0,0;
		f << 0,0,0;
		f0 << 0,0,0;
		return Vect3d(L/2,L/2,L/2);
	});

	/*
	 * setup gamma calculation
	 */
	ptr<GammaEval> gamma = ptr<GammaEval>(new GammaEval(params->dem_diameter));
	gamma->reset_limits(params->sph_maxh*0.5,params->sph_maxh*1.5,0.0,3.0*params->sph_maxh+params->dem_diameter/2.0,100,100);


	std::cout << "starting...."<<std::endl;
	sph->init_neighbour_search(min,max,2*params->sph_maxh,periodic);
	dem->init_neighbour_search(min,max,2*params->sph_maxh + params->dem_diameter/2.0,periodic);

	/*
	 * create sph and dem particles
	 */
	sph->create_particles_grid(min,max,Vect3i(nx,nx,nx+3),[dem,psep,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		h = params->sph_maxh;
		omega = 1.0;
		kappa = 0.0;
		v << 0,0,params->sph_inletv;
		v0 << 0,0,params->sph_inletv;
		dddt = 0;
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
	
		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 < pow(params->dem_diameter/2.0,2)) {
				std::cout << "removing sph at r = "<<r<<" |r| = "<<sqrt(r2)<<"dist = "<<2*h+params->dem_diameter/2.0<<" dx = "<<dx<<std::endl;
				i.mark_for_deletion();
			}
		}
	});
	sph->delete_particles();


	/*
	 * init porosity and rho
	 */
	//std::cout << "calculate omega"<<std::endl;
	std::for_each(sph->begin(),sph->end(),[gamma,sph,dem,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		e = 1;
		//bool found = false;
		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > pow(2*h + params->dem_diameter/2.0,2)) continue;
			const double absr = sqrt(r2);
			e -= gamma->get_gamma(h,absr);
			//std::cout <<" particle at r = "<<r<<" absr = "<<absr<<" dx = "<<dx<<" e = "<<e<<std::endl;

			//found = true;
		}
		//h = pow(params->sph_mass/(e*rho),1.0/NDIM);

		//rho = params->sph_dens/e;
	});

	/*
	 * setup output stuff
	 */
	auto sph_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto dem_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	sph->copy_to_vtk_grid(sph_grid);
	dem->copy_to_vtk_grid(dem_grid);

	char sph_name[100];
	sprintf(sph_name, "%s/at_start_sph", dir);
	char dem_name[100];
	sprintf(dem_name, "%s/at_start_dem", dir);
	Visualisation::vtkWriteGrid(sph_name,0,sph_grid);
	Visualisation::vtkWriteGrid(dem_name,0,dem_grid);

	double next_layer_time = 0.5*psep/params->sph_inletv;
	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			//std::this_thread::sleep_for(std::chrono::seconds(1));
			sphdem(sph,dem,params,sph_geometry,dem_geometry,gamma);

			/*
			 * convert fixed to nonfixed sph particles
			 */
			std::for_each(sph->begin(),sph->end(),[L,gamma,sph,dem,params](SphType::Value& i) {
					REGISTER_SPH_PARTICLE(i);
					if (r[2]<0) {
						fixed = true;
					} else {
						fixed = false;
					}
					if (r[2]>L) {
						i.mark_for_deletion();
					}
			});
			sph->delete_particles();


			/*
			 * add new layer if neccessary
			 */
			if (params->time > next_layer_time) {
				sph->create_particles_grid(min,Vect3d(max[0],max[1],min[2]),Vect3i(nx,nx,1),[dem,psep,params](SphType::Value& i) {
					REGISTER_SPH_PARTICLE(i);
					h = params->sph_maxh;
					omega = 1.0;
					kappa = 0.0;
					v << 0,0,params->sph_inletv;
					v0 << 0,0,params->sph_inletv;
					dddt = 0;
					e = 1;
					rho = params->sph_dens;
					f << 0,0,0;
					fdrag << 0,0,0;
					f0 << 0,0,0;
					fixed = true;
				});
				next_layer_time += psep / params->sph_inletv;
			}
		}
		std::cout <<"iteration "<<i<<std::endl;

		sph->copy_to_vtk_grid(sph_grid);
		dem->copy_to_vtk_grid(dem_grid);
		sprintf(sph_name, "%s/sph", dir);
		sprintf(dem_name, "%s/dem", dir);
		Visualisation::vtkWriteGrid(sph_name,i,sph_grid);
		Visualisation::vtkWriteGrid(dem_name,i,dem_grid);
	}
	
	
}
