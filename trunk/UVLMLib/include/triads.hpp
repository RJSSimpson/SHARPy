/*
 * triads.hpp
 *
 *  Created on: 31 Jan 2013
 *      Author: rjs10
 */

#ifndef TRIADS_HPP_
#define TRIADS_HPP_

void AddTriad(const double* x1, const double* x2, double* xOut);

void SubTriad(const double* x1, const double* x2, double* xOut);

void MulTriad(const double* x1, const double Factor, double* xOut);

double DotTriad(const double* x1, const double* x2);

double NormTriad(const double* x1);

void NormaliseTriad(const double* x1, double* xOut);

void CrossTriad(const double* x1, const double* x2, double* xOut);

void BilinearMapTriad(const double* p1, const double* p2, \
					  const double* p3, const double* p4, \
					  double* pOut);

void CopyTriad(double* pTarget, double* pSrc);

void PrintTriad(double* x);


#endif /* TRIADS_HPP_ */
