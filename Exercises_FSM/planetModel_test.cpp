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

void calcAcc(object* o, int nbody) {
    double g, temp_vel1, temp_vel2;
    for(int i = 0; i < nbody - 1; i++) {
        /*o[i].vel_old[0] = o[i].vel[0];
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
        o[i + 1].pos[2] = o[i + 1].vel_old[2];*/

        for(int j = 1; j < nbody; j++) {
            if(i == j) continue;
            g = G / (std::pow(o[i].pos[0] - o[j].pos[0], 2) + std::pow(o[i].pos[1] - o[j].pos[1], 2) + std::pow(o[i].pos[2] - o[j].pos[2], 2));
            for(int k = 0; k < 3; k++) {
                temp_vel1 = o[i].vel[k];
                temp_vel2 = o[j].vel[k];

                o[i].vel[k] = - g * o[i].mass * (o[i].pos[k] - o[j].pos[k]);
                o[j].vel[k] = g * o[j].mass * (o[i].pos[k] - o[j].pos[k]);
            
                o[i].pos[k] = temp_vel1;
                o[j].pos[k] = temp_vel2;
            }
            
        }
    }
}


void calcVelAcc(object* o, int num, int k, int nbody) {
    double g, temp_vel1, temp_vel2;
    for(int i = num + 1; i < nbody; i++) {
        g = G / (std::pow(o[num].pos[0] - o[i].pos[0], 2) + std::pow(o[num].pos[1] - o[i].pos[1], 2) + std::pow(o[num].pos[2] - o[i].pos[2], 2));
                
        temp_vel1 = o[num].vel[k];
        temp_vel2 = o[i].vel[k];

        o[num].vel[k] = - g * o[num].mass * (o[num].pos[k] - o[i].pos[k]);
        o[i].vel[k] = g * o[i].mass * (o[num].pos[k] - o[i].pos[k]);
            
        o[num].pos[k] = temp_vel1;
        o[i].pos[k] = temp_vel2;
    }
}


void leapfrog(object* o, int nbody, double dt) {
    double x_half, v, a;
    for(int i = 0; i < nbody; i++) {
        for(int k = 0; k < 3; k++) {
            x_half = o[i].pos[k] + (dt / 2) * o[i].vel[k];
            v = o[i].vel[k];
            calcVelAcc(o, i, k, nbody);
            o[i].vel[k] = v + dt * o[i].vel[k];
            o[i].pos[k] = x_half + (dt / 2) * o[i].vel[k];
        }   
    }
} 

int main() {
    

    return 0;
}