/*
 * VLM.hpp
 *
 *  Created on: 30 Jan 2013
 *      Author: rjs10
 */

#ifndef VLM_HPP_
#define VLM_HPP_
#include <datatypesx.hpp>

void cpp_solver_vlm(const double* Zeta, const double* ZetaDot, \
				    const double* Uext, \
				    double* ZetaStar, \
				    VMopts VMOPTS, \
				    double* Forces, \
				    double* Gamma_Vec, double* GammaStar_Vec, \
				    double* AIC_Vec,\
				    double* BIC_Vec);

#endif /* VLM_HPP_ */
