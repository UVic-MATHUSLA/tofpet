#include "CoarseSorter.hpp"
#include <algorithm>
#include <functional>
#include <vector>

using namespace std;
using namespace PETSYS;

CoarseSorter::CoarseSorter(EventSink<RawHit> *sink): UnorderedEventHandler<RawHit, RawHit>(sink) { resetCounters(); }

struct SortEntry {
	long long time;
	RawHit *p;
};

static bool operator < (SortEntry lhs, SortEntry rhs) { return lhs.time < rhs.time; }

EventBuffer<RawHit> *CoarseSorter::handleEvents(EventBuffer<RawHit> *inBuffer) {
	unsigned N = inBuffer->getSize();										// get number of events in input buffer
	EventBuffer<RawHit> *outBuffer = new EventBuffer<RawHit>(N, inBuffer);	// allocate output buffer with same size as input buffer
	u_int64_t lSingleRead = 0;												// local counter for number of events passed through this handler

	vector<SortEntry> sortList;	// create a vector to hold the events and their times for sorting
	sortList.reserve(N);		// reserve space in the vector to avoid reallocations during push_back

	auto pi = inBuffer->getPtr();						// get pointer to the first event in the input buffer
	auto pe = pi + N;									// get pointer to the end of the input buffer (one past the last event)
	for(; pi < pe; pi++) {								// populate the sortList with events and their times
		SortEntry entry = {.time = pi->time, .p = pi};	// create a SortEntry with the event's time and a pointer to the event
		sortList.push_back(entry);						// add the SortEntry to the sortList
	}

	sort(sortList.begin(), sortList.end());	// sort the sortList based on the event times using the overloaded < operator defined above

	auto po = outBuffer->getPtr();	// get pointer to the first event in the output buffer
	for(auto iter = sortList.begin(); iter != sortList.end(); iter++) {
		auto p = (*iter).p;			// get the pointer to the RawHit from the sorted list
		*po = *p;					// copy the RawHit to the output buffer
		po++;						// increment the output buffer pointer to the next slot
		lSingleRead++;				// increment the local counter for the number of events passed through this handler
	}

	atomicAdd(nSingleRead, lSingleRead);// atomically add the local counter to the class member counter for the total number of events passed through this handler
	outBuffer->setUsed(lSingleRead);	// set the number of used slots in the output buffer to the number of events passed through this handler
	return outBuffer;
}

void CoarseSorter::resetCounters() {
	nSingleRead = 0;
	UnorderedEventHandler<RawHit, RawHit>::resetCounters();
}

void CoarseSorter::report() {
	u_int64_t nTotal = nSingleRead;
	fprintf(stderr, ">> CoarseSorter report\n");
	fprintf(stderr, "   events passed\n");
	fprintf(stderr, "  %13lu\n", nSingleRead);
	UnorderedEventHandler<RawHit, RawHit>::report();
}
