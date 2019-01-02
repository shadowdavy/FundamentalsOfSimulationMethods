#include <iostream>
#include <vector>
#include <math.h>

static const double G  = 6.67384e-11;
static const double ms = 1.99e30;
static const double me = 1;
static const double mm = 1;

typedef struct {
    double mass;
    double pos[3];
    double vel[3];
    double vel_old[3];
    double acc[3];
} object;

void calcAcc(object* o) {
    double g;
    for(int i = 0; i < 2; i++) {
        o[i].vel_old[0] = o[i].vel[0];
        o[i].vel_old[1] = o[i].vel[1];
        o[i].vel_old[2] = o[i].vel[2];
        o[i + 1].vel_old[0] = o[i + 1].vel[0];
        o[i + 1].vel_old[1] = o[i + 1].vel[1];
        o[i + 1].vel_old[2] = o[i + 1].vel[2];

        g = G / (std::pow(o[i].pos[0] - o[i + 1].pos[0], 2) + std::pow(o[i].pos[1] - o[i + 1].pos[1], 2) + std::pow(o[i].pos[2] - o[i + 1].pos[2], 2));
        o[i].vel[0] = - g * o[i].mass * (o[i].pos[0] - o[i + 1].pos[0]);
        o[i].vel[1] = - g * o[i].mass * (o[i].pos[1] - o[i + 1].pos[1]);
        o[i].vel[2] = - g * o[i].mass * (o[i].pos[2] - o[i + 1].pos[2]);

        o[i + 1].vel[0] = g * o[i + 1].mass * (o[i].pos[0] - o[i + 1].pos[0]);
        o[i + 1].vel[1] = g * o[i + 1].mass * (o[i].pos[1] - o[i + 1].pos[1]);
        o[i + 1].vel[2] = g * o[i + 1].mass * (o[i].pos[2] - o[i + 1].pos[2]);

        o[i].pos[0] = o[i].vel_old[0];
        o[i].pos[1] = o[i].vel_old[1];
        o[i].pos[2] = o[i].vel_old[2];

        o[i + 1].pos[0] = o[i + 1].vel_old[0];
        o[i + 1].pos[1] = o[i + 1].vel_old[1];
        o[i + 1].pos[2] = o[i + 1].vel_old[2];
    }
}

std::vector<double> leapfrog(function f, std::vector<double> y, double dt) {
    std::vector<double> y_firstHalf;
    std::vector<double> y_secondHalf;
    std::vector<double> x_half;
    std::vector<double> x;
    std::vector<double> v;
    int vecSize = y.size();

    for(int i = 0; i < vecSize / 2; i++) {
        y_firstHalf.push_back(y[i]);
    }
    for(int i = vecSize / 2; i < vecSize; i++) {
        y_secondHalf.push_back(y[i]);
    }
    
    for(int i = 0; i < vecSize / 2; i++) {
        x_half.push_back(y_secondHalf[i] + (dt / 2) * y_firstHalf[i]);
    }
    for(int i = 0; i < vecSize / 2; i++) {
        
    }
    
    return -1;
} 