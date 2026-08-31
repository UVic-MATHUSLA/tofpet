#include "CoincidenceGrouper.hpp"
#include <algorithm>
#include <math.h>
#include <vector>

using namespace PETSYS;

CoincidenceGrouper::CoincidenceGrouper(SystemConfig *systemConfig, EventSink<Coincidence> *sink):
	systemConfig(systemConfig), UnorderedEventHandler<GammaPhoton, Coincidence>(sink) {
	resetCounters();
}

CoincidenceGrouper::~CoincidenceGrouper() {}

EventBuffer<Coincidence> *CoincidenceGrouper::handleEvents(EventBuffer<GammaPhoton> *inBuffer) {
	double cWindow = systemConfig->sw_trigger_coincidence_time_window;   // Get coincidence time window from config
	double Tps = 5000;                                                   // harcoded clock time in ps!!!!!!!!!!
	long long tMin = inBuffer->getTMin() * (long long) Tps;   // Get minimum time from input buffer and convert to ps
	unsigned N = inBuffer->getSize();                         // Get number of events in input buffer
	EventBuffer<Coincidence> *outBuffer = new EventBuffer<Coincidence>(N, inBuffer);   // Output buffer with same size as input buffer
	u_int64_t lPrompts = 0;          // Number of prompt coincidences found
	u_int64_t lPromptsRegion = 0;    // Number of invalid region prompt coincidences found
	u_int64_t lPromptsTime = 0;      // Number of invalid time prompt coincidences found
	u_int64_t lCoincPhotopeak = 0;   // Number of coincidences that are photopeak events
	

	// Loop through each photon in input buffer to find coincidences with other photons based on time and region
	for(unsigned i = 0; i < N; i++) {
		GammaPhoton &photon1 = inBuffer->get(i);   // Get current photon from input buffer
		// Loop through remaining photons in input buffer to find coincidences with current photon
		for(unsigned j = i + 1; j < N; j++) {
			GammaPhoton &photon2 = inBuffer->get(j);   // Get next photon from input buffer to compare with current photon
			// If time diff between photons is greater than coincidence window + maximum unorder threshold break loop
			if((photon2.time - photon1.time) > (cWindow + MAX_UNORDER)) { 
				lPromptsTime++;
				break;
			}
			// If regions of the two photons are not allowed to coincide based on system config, skip to next iteration
			if(!systemConfig->isCoincidenceAllowed(photon1.region, photon2.region)) { 
				lPromptsRegion++;
				continue; 
			}

			// If time diff between photons is within coincidence window
			if(fabs(photon1.time - photon2.time) <= cWindow) {
				Coincidence &c = outBuffer->getWriteSlot();      // Get write slot in output buffer for coinc
				c.nPhotons = 2;                                  // Set number of photons in coinc to 2
				bool first1 = photon1.region > photon2.region;   // Order of photons based on region for consistency
				c.photons[0] = first1 ? &photon1 : &photon2;     // Assign first photon pointer in coinc
				c.photons[1] = first1 ? &photon2 : &photon1;     // Assign second photon pointer in coinc
				c.valid = true;                                  // Mark coinc as valid
				outBuffer->pushWriteSlot();   // Push write slot to output buffer to finalize coinc
				lPrompts++;                   // Increment local counter for number of coincs found
			}
		}
	}
	atomicAdd(nPrompts, lPrompts);
	atomicAdd(nPromptsRegion, lPromptsRegion);
	atomicAdd(nPromptsTime, lPromptsTime);
	atomicAdd(nCoincPhotopeak, lCoincPhotopeak);
	return outBuffer;
}

void CoincidenceGrouper::resetCounters() {
	nPrompts = 0;
	UnorderedEventHandler<GammaPhoton, Coincidence>::report();
}

void CoincidenceGrouper::report() {
	fprintf(stderr, ">> CoincidenceGrouper report\n");
	fprintf(stderr, "   events passed:\n");
	fprintf(stderr, "   %13lu \n", nPrompts);
	fprintf(stderr, "   %13lu (%4.1f%%) with photopeak events\n", nCoincPhotopeak, 100.0 * nCoincPhotopeak / nPrompts);
	fprintf(stderr, "   events rejected:\n");
	fprintf(stderr, "   %13lu (%4.1f%%) with invalid regions\n", nPromptsRegion, 100.0 * nPromptsRegion / nPrompts);
	fprintf(stderr, "   %13lu (%4.1f%%) with invalid times\n", nPromptsTime, 100.0 * nPromptsTime / nPrompts);
	UnorderedEventHandler<GammaPhoton, Coincidence>::report();
}
