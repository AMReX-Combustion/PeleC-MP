
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <Transport_F.H>
#include <drag_F.H>

using namespace amrex;


void
SprayParticleContainer::insertParticles (Real time, int nstep, int lev)
{

}


void
SprayParticleContainer::injectParticles (Real time, int nstep, int lev)
{

} 

void
SprayParticleContainer::InitParticlesUniform(const int& lev, const int& num_ppc)
{

  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    
    const Box& tile_box  = mfi.tilebox();

    Cuda::HostVector<ParticleType> host_particles;
    std::array<Cuda::HostVector<Real>, NSR_SPR> host_particles_rdata;

    for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
      for (int i_part=0; i_part<num_ppc;i_part++) {

	Real r = (rand()%100)/99.;
	Real x = plo[0] + (iv[0] + r)*dx[0];

	r = (rand()%100)/99.;
	Real y = plo[1] + (iv[1] + r)*dx[1];

	ParticleType p;
	p.id()  = ParticleType::NextID();
	p.cpu() = ParallelDescriptor::MyProc();
	p.pos(0) = x;
	p.pos(1) = y;
	if (AMREX_SPACEDIM>2) {
	  r = (rand()%100)/99.;
	  Real z = plo[2] + (iv[2] + r)*dx[2];
	  p.pos(2) = z;
	}

	host_particles.push_back(p);

	// fill in NSR_SPR items of particle real data
	host_particles_rdata[0].push_back( 0. );//u-vel
	host_particles_rdata[1].push_back( 0. );//v-vel
	if (AMREX_SPACEDIM>2) {
	  host_particles_rdata[2].push_back( 0. );//w-vel
	}    
	host_particles_rdata[AMREX_SPACEDIM].push_back( 293. ); // temperature
	host_particles_rdata[AMREX_SPACEDIM+1].push_back( 0.0200 ); // diameter
	host_particles_rdata[AMREX_SPACEDIM+2].push_back( 0.68141 ); // fuel density


      }   
    }

    auto& particles = GetParticles(lev);
    auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
    auto old_size = particle_tile.GetArrayOfStructs().size();
    auto new_size = old_size + host_particles.size();
    particle_tile.resize(new_size);

    //Copy the ArrayOfStructs part of host particles to the GPU
    Cuda::thrust_copy(host_particles.begin(),
		      host_particles.end(),
		      particle_tile.GetArrayOfStructs().begin() + old_size);

    //Copy the StructOfArrays part i.e. the rdata of host particles to GPU
    for (int kk = 0; kk < NSR_SPR; ++kk) {
      
      Cuda::thrust_copy(host_particles_rdata[kk].begin(),
			host_particles_rdata[kk].end(),
			particle_tile.GetStructOfArrays().GetRealData(kk).begin() + old_size);
      
    }

  }
    
}
