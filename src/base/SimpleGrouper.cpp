#include "SimpleGrouper.hpp"
#include <algorithm>
#include <math.h>
#include <vector>

using namespace PETSYS;
using namespace std;

SimpleGrouper::SimpleGrouper(SystemConfig *systemConfig, EventSink<GammaPhoton> *sink):
	systemConfig(systemConfig), UnorderedEventHandler<Hit, GammaPhoton>(sink) {
	resetCounters();
}

SimpleGrouper::~SimpleGrouper() {}

auto comp = [](Hit *a, Hit *b) { return a->energy > b->energy; };   // Define a custom comparator to sort by energy in
																	// descending order

EventBuffer<GammaPhoton> *SimpleGrouper::handleEvents(EventBuffer<Hit> *inBuffer) {
	double timeWindow1 = systemConfig->sw_trigger_group_time_window;   // Get the time window for grouping hits into a photon
	// Get the square of the maximum distance for grouping hits into a photon
	float radius2 = (systemConfig->sw_trigger_group_max_distance) * (systemConfig->sw_trigger_group_max_distance);
	float minEnergy = systemConfig->sw_trigger_group_min_energy;   // Get the minimum energy threshold for a photon
	float maxEnergy = systemConfig->sw_trigger_group_max_energy;   // Get the maximum energy threshold for a photon
	int maxHits = systemConfig->sw_trigger_group_max_hits;         // Get the maximum number of hits allowed in a photon
	if(maxHits > GammaPhoton::maxHits)
		maxHits = maxHits;   // Ensure that maxHits does not exceed the maximum allowed hits in a GammaPhoton
	int minHits = systemConfig->sw_trigger_group_min_hits;   // Get the minimum number of hits required for a photon

	u_int64_t lPhotonsHits[maxHits];   // Local array to count the number of photons with a specific number of hits
	for(int i = 0; i < maxHits; i++) { lPhotonsHits[i] = 0; }   // Initialize the local photon hit count array to zero

	u_int64_t lHitsReceived = 0;           // Number of hits received
	u_int64_t lHitsReceivedValid = 0;      // Number of valid hits received
	u_int64_t lPhotonsFound = 0;           // Number of photons found
	u_int64_t lPhotonsHitsOverflow = 0;    // Number of photons with hits exceeding the maximum allowed
	u_int64_t lPhotonsHitsUnderflow = 0;   // Number of photons with hits below the minimum required
	u_int64_t lPhotonsLowEnergy = 0;       // Number of photons with energy below the minimum threshold
	u_int64_t lPhotonsHighEnergy = 0;      // Number of photons with energy above the maximum threshold
	u_int64_t lPhotonsPassed = 0;          // Number of photons that passed all criteria and were accepted

	unsigned N = inBuffer->getSize();
	EventBuffer<GammaPhoton> *outBuffer = new EventBuffer<GammaPhoton>(N, inBuffer);
	vector<bool> taken(N, false);
	Hit *hits[maxHits];

	// Loop through each hit in input buffer to group them into photons based on time and spatial proximity
	for(unsigned i = 0; i < N; i++) {
		Hit &hit = inBuffer->get(i);   // Get current hit from input buffer
		lHitsReceived += 1;            // Increment local counter for number of hits received

		if(!hit.valid) continue;   // If hit is not valid, skip to next iteration
		lHitsReceivedValid += 1;   // Increment local counter for number of valid hits received

		if(taken[i]) continue;   // If hit has already been grouped into a photon (taken), skip to next iteration
		taken[i] = true;         // Mark current hit as taken to avoid processing it again

		uint8_t eventFlags = 0x0;   // Initialize event flags to zero for current photon
		hits[0] = &hit;             // Store pointer to the current hit in hits array for current photon
		int nHits = 1;              // Initialize number of hits for current photon to 1 (the current hit)

		// Loop through remaining hits in input buffer to find hits that can be grouped with current hit into a photon
		for(int j = i + 1; j < N; j++) {
			Hit &hit2 = inBuffer->get(j);   // Get next hit from input buffer to compare with current hit
			if(!hit2.valid) continue;       // If next hit is not valid, skip to next iteration
			if(taken[j])continue;           // If next hit has already been taken, skip to next iteration

			// If time difference between two hits is greater than time window, stop searching for more hits for this photon
			if((hit2.time - hit.time) > (timeWindow1 + MAX_UNORDER)) break;
			// If two hits aren't allowed to be grouped together, skip to the next iteration
			if(!systemConfig->isMultiHitAllowed(hit2.region, hit.region)) continue;
			// If time difference between two hits is greater than time window, skip to the next iteration
			if(fabs(hit.time - hit2.time) > timeWindow1) continue;

			float u = hit.x - hit2.x;
			float v = hit.y - hit2.y;
			float w = hit.z - hit2.z;
			float d2 = u * u + v * v + w * w;
			if(d2 > radius2) continue;

			taken[j] = true;
			if(nHits >= maxHits) { nHits++; }   // Increase hit count but don't actually add a hit
			else {
				hits[nHits] = &hit2;   // Add pointer to next hit to hits array for current photon
				nHits++;               // Increment number of hits for current photon
			}
		}

		if(nHits > maxHits) {
			eventFlags |= 0x1;   // Flag event as having excessive hits
			nHits = maxHits;     // Set number of hits to maximum hits, as code below depends on it
		}
		else if(nHits < minHits) { eventFlags |= 0x8; }

		bool sorted = false;
		std::sort(hits, hits + nHits, comp);   // Sort to put highest energy event first

		float totalEnergy = 0;
		// Calculate total energy and assemble the output structure (input = hit, output = photon)
		GammaPhoton &photon = outBuffer->getWriteSlot();
		for(int k = 0; k < nHits; k++) {
			photon.hits[k] = hits[k];
			totalEnergy += photon.hits[k]->energy;
		}

		photon.nHits = nHits;
		photon.region = photon.hits[0]->region;
		photon.time = photon.hits[0]->time;
		photon.x = photon.hits[0]->x;
		photon.y = photon.hits[0]->y;
		photon.z = photon.hits[0]->z;
		photon.energy = totalEnergy;

		if(photon.energy < minEnergy) eventFlags |= 0x2;
		if(photon.energy > maxEnergy) eventFlags |= 0x4;

		// Count photons
		lPhotonsFound += 1;
		if((eventFlags & 0x1) == 0) { lPhotonsHits[photon.nHits - 1] += 1; }
		else { lPhotonsHitsOverflow += 1; }

		if((eventFlags & 0x8) != 0) lPhotonsHitsUnderflow += 1;
		if((eventFlags & 0x2) != 0) lPhotonsLowEnergy += 1;
		if((eventFlags & 0x4) != 0) lPhotonsHighEnergy += 1;

		if(eventFlags == 0) {
			lPhotonsPassed += 1;
			photon.valid = true;
			outBuffer->pushWriteSlot();
		}
	}

	for(int i = 0; i < maxHits; i++) atomicAdd(nPhotonsHits[i], lPhotonsHits[i]);

	atomicAdd(nHitsReceived, lHitsReceived);
	atomicAdd(nHitsReceivedValid, lHitsReceivedValid);
	atomicAdd(nPhotonsFound, lPhotonsFound);
	atomicAdd(nPhotonsHitsOverflow, lPhotonsHitsOverflow);
	atomicAdd(nPhotonsHitsUnderflow, lPhotonsHitsUnderflow);
	atomicAdd(nPhotonsLowEnergy, lPhotonsLowEnergy);
	atomicAdd(nPhotonsHighEnergy, lPhotonsHighEnergy);
	atomicAdd(nPhotonsPassed, lPhotonsPassed);

	return outBuffer;
}

void SimpleGrouper::resetCounters() {
	for(int i = 0; i < GammaPhoton::maxHits; i++) nPhotonsHits[i] = 0;
	nHitsReceived = 0;
	nHitsReceivedValid = 0;
	nPhotonsFound = 0;
	nPhotonsHitsOverflow = 0;
	nPhotonsHitsUnderflow = 0;
	nPhotonsLowEnergy = 0;
	nPhotonsHighEnergy = 0;
	nPhotonsPassed = 0;
	UnorderedEventHandler<Hit, GammaPhoton>::resetCounters();
}

void SimpleGrouper::report() {
	int maxHits = systemConfig->sw_trigger_group_max_hits;
	if(maxHits > GammaPhoton::maxHits) maxHits = maxHits;
	int minHits = systemConfig->sw_trigger_group_min_hits;

	fprintf(stderr, ">> SimpleGrouper report\n");
	fprintf(stderr, "   hits received:\n");
	fprintf(stderr, "   %13lu total\n", nHitsReceived);
	fprintf(stderr, "   %13lu (%4.1f%%) invalid\n", nHitsReceived - nHitsReceivedValid, 100.0 * (nHitsReceived - nHitsReceivedValid) / nHitsReceived);
	fprintf(stderr, "   photons found:\n");
	fprintf(stderr, "   %13lu total\n", nPhotonsFound);
	for(int i = 0; i < maxHits; i++) {
		float fraction = nPhotonsHits[i] / ((float) nPhotonsFound);
		if(fraction > 0.05) { fprintf(stderr, "   %13lu (%4.1f%%) with %d hits\n", nPhotonsHits[i], 100.0 * fraction, i + 1); }
	}
	fprintf(stderr, "   %13.1f hits/photon\n", float(nHitsReceived) / nPhotonsFound);
	fprintf(stderr, "   photons rejected:\n");
	fprintf(stderr, "   %13lu (%4.1f%%) with more than %d hits\n", nPhotonsHitsOverflow, 100.0 * nPhotonsHitsOverflow / nPhotonsFound, maxHits);
	fprintf(stderr, "   %13lu (%4.1f%%) with less than %d hits\n", nPhotonsHitsUnderflow, 100.0 * nPhotonsHitsUnderflow / nPhotonsFound, minHits);
	fprintf(stderr, "   %13lu (%4.1f%%) failed minimum energy\n", nPhotonsLowEnergy, 100.0 * nPhotonsLowEnergy / nPhotonsFound);
	fprintf(stderr, "   %13lu (%4.1f%%) failed maximim energy\n", nPhotonsHighEnergy, 100.0 * nPhotonsHighEnergy / nPhotonsFound);
	fprintf(stderr, "   photons passed:\n");
	fprintf(stderr, "   %13lu (%4.1f%%) passed\n", nPhotonsPassed, 100.0 * nPhotonsPassed / nPhotonsFound);
	UnorderedEventHandler<Hit, GammaPhoton>::report();
}

