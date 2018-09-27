#include <vector>
#include "PhotosHepMCParticle.h"
#include "PhotosHepMCEvent.h"
#include "Log.h"
using namespace std;

namespace Photospp
{

PhotosHepMCEvent::PhotosHepMCEvent(HepMC::GenEvent * event)
{
	m_event=event;
	HepMC::GenEvent::particle_const_iterator part_itr = m_event->particles_begin();
	for( ; part_itr!=m_event->particles_end(); part_itr++)
	{
		PhotosParticle *particle = new PhotosHepMCParticle(*part_itr);
		particles.push_back(particle);
	}
    
	switch(m_event->momentum_unit()) {
		case HepMC::Units::GEV:
			Photos::setMomentumUnit(Photos::GEV);
			break;
		case HepMC::Units::MEV:
			Photos::setMomentumUnit(Photos::MEV);
			break;
		default:
			Log::Error()<<"PhotosHepMCEvent: undefined unit, important for pair emission only"<<endl;
			Photos::setMomentumUnit(Photos::DEFAULT_MOMENTUM);
			break;
	};
}

PhotosHepMCEvent::~PhotosHepMCEvent()
{
	while(particles.size())
	{
		PhotosParticle *p = particles.back();
		particles.pop_back();
		if(p) delete p;
	}
}

HepMC::GenEvent * PhotosHepMCEvent::getEvent()
{
	return m_event;
}

void PhotosHepMCEvent::print()
{
	if(!m_event) return;
	m_event->print();
}

vector<PhotosParticle*> PhotosHepMCEvent::getParticleList()
{
	return particles;
}

} // namespace Photospp
