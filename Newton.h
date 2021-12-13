#ifndef NEWTON_METHOD_NEWTON_H
#define NEWTON_METHOD_NEWTON_H

#include <Eigen/Dense>
#include <iostream>

class Newton {
public:
    Newton() {}

    double solve(double u0, double (*f)(double), double (*df)(double), double epsilon) {
        double u = u0 - f(u0) / df(u0);
        while (fabs(u - u0) > epsilon) {
            u0 = u;
            u = u0 - f(u0) / df(u0);
        }
        return u;
    }

    double solve(double u0, double u1, double (*f)(double), double epsilon) {
        double u = u1 - f(u1) * (u1 - u0) / (f(u1) - f(u0));
        while (fabs(u - u1) > epsilon) {
            u0 = u1;
            u1 = u;
            u = u1 - f(u1) * (u1 - u0) / (f(u1) - f(u0));
        }
        return u;
    }
};


#endif //NEWTON_METHOD_NEWTON_H
