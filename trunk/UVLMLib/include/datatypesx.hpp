/*
 * datatypesx.h
 *
 *  Created on: 30 Jan 2013
 *      Author: rjs10
 */

#ifndef DATATYPESX_H_
#define DATATYPESX_H_

#include <vector>
using std::vector;

class VMopts {
public:
	unsigned int M;
	unsigned int N;
	bool ImageMethod;
};

class Array3D {
public:
	Array3D(const unsigned int M, const unsigned int N, const unsigned int O) {
		/** initialise*/
		Elem.resize(M);
			for (unsigned int i = 0; i < M; i++) {
				Elem[i].resize(N);
				for (unsigned int j = 0; j < N; j++) {
					Elem[i][j].resize(O);
				}
			}
	}
	vector<vector<vector<double> > > Elem;
};


#endif /* DATATYPESX_H_ */
