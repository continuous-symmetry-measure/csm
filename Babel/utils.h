#ifndef UTILS_H
#define UTILS_H

#include <math.h>

class TransformationMatrix {
protected:
	double matrix[3][3];

public:
	void transform(double *vec, double *out_vec) {
		for (int i = 0; i < 3; i++) {
			out_vec[i] = 0.0;
			for (int j = 0; j < 3; j++) {
				out_vec[i] += matrix[i][j] * vec[j];
			}
		}
	}
};

class RotationMatrix : public TransformationMatrix {
public:
	RotationMatrix(double *axis, double angle) {
		double ct = cos(angle);
		double st = sin(angle);
		matrix[0][0] = ct + axis[0] * axis[0] * (1-ct);
		matrix[0][1] = axis[0] * axis[1] * (1-ct) - axis[2] * st;
		matrix[0][2] = axis[0] * axis[2] * (1-ct) + axis[1] * st;
		matrix[1][0] = axis[0] * axis[1] * (1-ct) + axis[2] * st; 
		matrix[1][1] = ct + axis[1] * axis[1] * (1-ct);
		matrix[1][2] = axis[1] * axis[2] * (1-ct) - axis[0] * st;
		matrix[2][0] = axis[2] * axis[0] * (1-ct) - axis[1] * st;
		matrix[2][1] = axis[1] * axis[2] * (1-ct) + axis[0] * st; 
		matrix[2][2] = ct + axis[2] * axis[2] * (1-ct);
	}
};

class ReflectionMatrix : public TransformationMatrix {
public:
	/* axis is the axis perpendicular to the plane */ 
	ReflectionMatrix(double *axis) {
		for (int i = 0; i < 3;i++) {
			for (int j = 0; j < 3; j++) {
				matrix[i][j] = -2 * axis[i] * axis[j];
				if (i==j) matrix[i][j] += 1;
			}
		}
	}
};

#endif