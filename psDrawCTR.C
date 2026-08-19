#include <stdint.h>
#include <stdio.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TTree.h>

struct Event {
	uint8_t mh_n1;
	uint8_t mh_j1;
	long long time1;
	float e1;
	int id1;
	uint8_t mh_n2;
	uint8_t mh_j2;
	long long time2;
	float e2;
	int id2;
} __attribute__((__packed__));

int psDrawCTR(char const *filePrefix, Int_t channelA = -1, Int_t channelB = -1) {
	char fn[1024];
	sprintf(fn, "%s.ldat", filePrefix);   // get binary data file name
	FILE *dataFile = fopen(fn, "rb");     // open binary data file
	sprintf(fn, "%s.lidx", filePrefix);   // get index file name
	FILE *indexFile = fopen(fn, "rb");    // open index file

	Double_t tBinWidth = 30;         // time bin width in ps
	Double_t wMax = 5000;            // maximum time window in ps
	const Int_t N_CHANNELS = 1024;   // number of channels in the data

	// [channelA, channelB]
    //  [0,14],   [1,7],  [2,15],  [3,12],  [4,27],  [5,10],  [6,11],  [8,13],
    //  [9,50], [16,22], [17,20], [18,23], [19,21], [24,31], [25,32], [26,30],
    // [28,29], [33,41], [34,40], [35,37], [36,39], [38,57], [42,48], [43,49],
    // [44,47], [45,46], [51,60], [52,63], [53,56], [54,55], [58,61], [59,62]]
	// channels with the same index in the two arrays are coincident pairs to be analyzed for CTR
	int chsA[32] = {0, 1, 2, 3, 4, 5, 6, 8, 9,16,17,18,19,24,25,26,28,33,34,35,36,38,42,43,44,45,51,52,53,54,58,59};
	int chsB[32] = {14,7,15,12,27,10,11,13,50,22,20,23,21,31,32,30,29,41,40,37,39,57,48,49,47,46,60,63,56,55,61,62};

	// histogram for coincidence time difference
	TH1F *hDelta = new TH1F("hDelta", "Coincidence time difference", 2 * wMax / tBinWidth, -wMax, wMax);
	// histogram for 2D channel correlation
	TH2F *hC2 = new TH2F("hC2", "C2", N_CHANNELS, 0, N_CHANNELS, N_CHANNELS, 0, N_CHANNELS);

	int minToT = 0;                                                              // minimum ToT to consider for analysis
	int maxToT = 15;                                                             // maximum ToT to consider for analysis
	TH1F *hE1 = new TH1F("hE1", "E 1", (maxToT - minToT) * 5, minToT, maxToT);   // histogram for charge spectrum of channel 1
	TH1F *hE2 = new TH1F("hE2", "E 2", (maxToT - minToT) * 5, minToT, maxToT);   // histogram for charge spectrum of channel 2

	struct CoincEvent {
		Long64_t time1;   // time of event 1
		Long64_t time2;   // time of event 2
		Float_t tot1;     // ToT of event 1
		Float_t tot2;     // ToT of event 2
		Int_t id1;        // channel ID of event 1
		Int_t id2;        // channel ID of event 2
	};

	const Int_t MAX_EVENTS = 1000000;
	const Int_t MAX_EVENTS_TOTAL = 250000000;
	CoincEvent *eventBuffer = new CoincEvent[MAX_EVENTS];   // buffer to hold events for analysis

	long long stepBegin, stepEnd;
	float step1, step2;

	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = NULL;

	sprintf(fn, "%s_stats.txt", filePrefix);
	FILE *stats = fopen(fn, "w");

	int stepIndex = 0;
	while(fscanf(indexFile, "%lld\t%lld\t%f\%f\n", &stepBegin, &stepEnd, &step1, &step2) == 4) {
		stepIndex += 1;
		hDelta->Reset();
		hC2->Reset();
		hE1->Reset();
		hE2->Reset();

		fseek(dataFile, stepBegin, SEEK_SET);
		Int_t count = 0;
		Event ei;
		int nRead = 0;

		// iterate through indeces of channelA and channelB to find the maximum bin in the 2D histogram hC2
		printf("Coincident channel pairs:\n");
		for(int i = 0; i < 32; i++) {
			printf("[%d,%d] ", chsA[i], chsB[i]);
		}
		printf("\n");

		if(channelA == -1 && channelB == -1) {
			while(ftell(dataFile) < stepEnd && fread(&ei, sizeof(ei), 1, dataFile) == 1 && count < MAX_EVENTS_TOTAL) {
				if(ei.e1 < minToT || ei.e2 < minToT) continue;
				hC2->Fill(ei.id1, ei.id2);
				count++;
				if(count % 1000000 == 0) cout << count << endl;
			}
			fseek(dataFile, stepBegin, SEEK_SET);
			Int_t binMax = hC2->GetMaximumBin();
			Int_t binx, biny, binz;
			hC2->GetBinXYZ(binMax, binx, biny, binz);
			channelA = binx - 1;
			channelB = biny - 1;
			printf("ChA=%d ChB=%d\n", channelA, channelB);
		}
		count = 0;
		char hName[128];

		while(ftell(dataFile) < stepEnd && fread(&ei, sizeof(ei), 1, dataFile) == 1 && count < MAX_EVENTS) {
			nRead++;
			if(nRead % 10000 == 0) { fflush(stdout); }
			if(ei.e1 < minToT || ei.e2 < minToT) continue;
			hC2->Fill(ei.id1, ei.id2);

			bool selectedPair = ((ei.id1 == channelA) || (ei.id1 == channelB)) && ((ei.id2 == channelA) || (ei.id2 == channelB));
			if(!selectedPair) { continue; }

			CoincEvent &e = eventBuffer[count];
			e.time1 = ei.time1;
			e.time2 = ei.time2;
			e.tot1 = ei.e1;
			e.tot2 = ei.e2;
			e.id1 = ei.id1;
			e.id2 = ei.id2;
			count++;

			hE1->Fill(ei.e1);
			hE2->Fill(ei.e2);
		}
		printf("\n");

		// hE1->Smooth();
		// hE2->Smooth();

		TSpectrum *spectrum = new TSpectrum(1, 3);
		spectrum->Search(hE1, 3, " ", 0.1);
		Int_t nPeaks = spectrum->GetNPeaks();
		if(nPeaks == 0) {
			printf("No peaks in hE1!!!\n");
			continue;
		}

	#if ROOT_VERSION_CODE > ROOT_VERSION(6, 0, 0)
		Double_t *xPositions1 = spectrum->GetPositionX();
	#else
		Float_t *xPositions1 = spectrum->GetPositionX();
	#endif
		Float_t x1_psc = 0;
		Float_t x1_pe = 0;

		for(Int_t i = 0; i < nPeaks; i++) {
			if(xPositions1[i] > x1_pe) {
				if(x1_pe > x1_psc) x1_psc = x1_pe;
				x1_pe = xPositions1[i];
			}
		}
		spectrum = new TSpectrum(1, 3);
		spectrum->Search(hE2, 3, " ", 0.1);
		nPeaks = spectrum->GetNPeaks();
		if(nPeaks == 0) {
			printf("No peaks in hE2!!!\n");
			continue;
		}

	#if ROOT_VERSION_CODE > ROOT_VERSION(6, 0, 0)
		Double_t *xPositions2 = spectrum->GetPositionX();
	#else
		Float_t *xPositions2 = spectrum->GetPositionX();
	#endif
		Float_t x2_psc = 0;
		Float_t x2_pe = 0;

		for(Int_t i = 0; i < nPeaks; i++) {
			if(xPositions2[i] > x2_pe) {
				if(x2_pe > x2_psc) x2_psc = x2_pe;
				x2_pe = xPositions2[i];
			}
		}
		Float_t x1 = x1_pe;
		Float_t x2 = x2_pe;

		hE1->Fit("gaus", "", "", x1 - 0.7, x1 + 0.7);
		if(hE1->GetFunction("gaus") == NULL) { return -1; }
		x1 = hE1->GetFunction("gaus")->GetParameter(1);
		Float_t sigma1 = hE1->GetFunction("gaus")->GetParameter(2);

		hE2->Fit("gaus", "", "", x2 - 0.7, x2 + 0.7);
		if(hE2->GetFunction("gaus") == NULL) { return -1; }
		x2 = hE2->GetFunction("gaus")->GetParameter(1);
		Float_t sigma2 = hE2->GetFunction("gaus")->GetParameter(2);

		// Redo, now with sN
		Float_t sN = 2.0;
		hE1->Fit("gaus", "", "", x1 - sN * sigma1, x1 + sN * sigma1);
		if(hE1->GetFunction("gaus") == NULL) { return -1; }
		x1 = hE1->GetFunction("gaus")->GetParameter(1);
		sigma1 = hE1->GetFunction("gaus")->GetParameter(2);

		hE2->Fit("gaus", "", "", x2 - sN * sigma2, x2 + sN * sigma2);
		if(hE2->GetFunction("gaus") == NULL) { return -1; }
		x2 = hE2->GetFunction("gaus")->GetParameter(1);
		sigma2 = hE2->GetFunction("gaus")->GetParameter(2);

		sN = 1;
		for(Int_t i = 0; i < count; i++) {
			CoincEvent &e = eventBuffer[i];
			if((e.tot1 < (x1 - sN * sigma1)) || (e.tot1 > (x1 + sN * sigma1))) continue;
			if((e.tot2 < (x2 - sN * sigma2)) || (e.tot2 > (x2 + sN * sigma2))) continue;
			Float_t delta = e.time1 - e.time2;
			hDelta->Fill(delta);
		}
		int binmax = hDelta->GetMaximumBin();
		Float_t hmax = hDelta->GetXaxis()->GetBinCenter(binmax);
		hDelta->Fit("gaus", "", "", hmax - 350, hmax + 350);

		char title[1024];

		gStyle->SetPalette(1);
		gStyle->SetOptFit(1);

		c2 = c1;
		c2->Clear();
		c2->Divide(3, 1);

		c2->cd(1);
		hE1->Draw();
		sprintf(title, "Charge Spectrum (Channel %d)", channelB);
		hE1->SetTitle(title);
		hE1->GetXaxis()->SetTitle("Charge [a.u.]");

		c2->cd(2);
		hDelta->Draw();
		hDelta->GetXaxis()->SetTitle("time1-time2 [ps]");

		c2->cd(3);
		hE2->Draw();
		sprintf(title, "Charge Spectrum (Channel %d)", channelA);
		hE2->SetTitle(title);
		hE2->GetXaxis()->SetTitle("Charge [a.u.]");

		c2->Modified();
		char pdfName[1024];
		sprintf(pdfName, "%s_%3.2f_%3.2f.pdf", filePrefix, step1, step2);
		c2->SaveAs(pdfName);

		float sigma = hDelta->GetRMS();
		float mean = hDelta->GetMean();
		float sigmaError = 0.1 * sigma;
		float meanError = 0.1 * mean;
		TF1 *g = hDelta->GetFunction("gaus");
		if(g != NULL) {
			mean = g->GetParameter(1);
			meanError = g->GetParError(1);
			sigma = g->GetParameter(2);
			sigmaError = g->GetParError(2);
		}
		int counts = hDelta->GetEntries();

		printf("---- step1 = %3.2f step2 = %3.2f ----\n", step1, step2);
		printf("CTR = %f sigma (%f FWHM)\n", sigma, 2.354 * sigma);
		printf("E1 = %f\n", x1);
		printf("E2 = %f\n", x2);
		fprintf(stats, "step1\tstep2\tsigma\tsigmaError\tmean\tmeanError\tx1\tx2\tcounts\n");
		fprintf(stats, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\n", step1, step2, sigma, sigmaError, mean, meanError, x1, x2, counts);
	}
	delete[] eventBuffer;
	return 0;
}
