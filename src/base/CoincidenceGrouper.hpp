#ifndef __PETSYS__COINCIDENCEGROUPER_HPP__DEFINED__
#define __PETSYS__COINCIDENCEGROUPER_HPP__DEFINED__

#include <Event.hpp>
#include <Instrumentation.hpp>
#include <SystemConfig.hpp>
#include <UnorderedEventHandler.hpp>

namespace PETSYS {
	class CoincidenceGrouper: public UnorderedEventHandler<GammaPhoton, Coincidence> {
	  public:
		CoincidenceGrouper(SystemConfig *systemConfig, EventSink<Coincidence> *sink);
		~CoincidenceGrouper();
		virtual void report();
		virtual void resetCounters();
	  private:
		virtual EventBuffer<Coincidence> *handleEvents(EventBuffer<GammaPhoton> *inBuffer);

		SystemConfig *systemConfig;
		u_int64_t nPrompts;
		u_int64_t nCoincPhotopeak;
		u_int64_t nListModeControl;
	};
}   // namespace PETSYS
#endif   // __PETSYS__COINCIDENCEGROUPER_HPP__DEFINED__
