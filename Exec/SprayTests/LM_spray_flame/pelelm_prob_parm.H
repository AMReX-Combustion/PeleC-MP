#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
    amrex::Real P_mean   = 101325.0_rt;
    amrex::Real standoff = 0.0_rt;
    amrex::Real pertmag  = 0.0004_rt;
    amrex::Real T_ox     = 900.;
    // Spray stuff
    amrex::Real jet_vel = 50.;
    amrex::Real jet_dia = 1.E-5;
    amrex::Real part_mean_dia = 1.E-6;
    amrex::Real part_stdev_dia = 0.;
    amrex::Real mass_flow_rate = 2.3E-3;
    amrex::Real part_temp = 300.;
    amrex::Real jet_start_time = 0.;
    amrex::Real jet_end_time = 10000.;
    amrex::Real spray_angle = 20.;
    amrex::Real jet_dx_mod = 50.;
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_jet = {{0.0}};
    unsigned int num_jets = 1;
    amrex::Gpu::HostVector<amrex::RealVect> jet_cents;
};
#endif
